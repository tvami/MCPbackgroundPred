import glob
import re
import ROOT, os, warnings, pandas, math, time
from PIL import Image
from collections import OrderedDict
from TwoDAlphabet.helpers import set_hist_maximums, execute_cmd, cd, hist2array
from TwoDAlphabet.binning import stitch_hists_in_x, convert_to_events_per_unit, get_min_bin_width
from TwoDAlphabet.ext import tdrstyle, CMS_lumi
from TwoDAlphabet.plotstyle import * # dictionaries containing plotting styles for mplhep
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import mplhep as hep

class Plotter(object):
    '''Class to manage output distributions, manipulate them, and provide access to plotting
    standard groups of distributions.

    Attributes:
        fit_tag (str): Either 's' or 'b'.
        fit_results (RooFitResult): The RooFitResult corresponding to the fit_tag.
        signal_strength (float): The post-fit signal strength.
        twoD (TwoDAlphabet): TwoDAlphabet object storing various meta information needed for access.
        yaxis1D_title (str): Title for "counts" axis of 1D plots. Defaults to 'Events / bin' but can change if plotting events per unit.
        df (pandas.DataFrame): DataFrame organizing the post-fit and pre-fit plots and their 1D projections.
        dir (str): Directory path to save final images.
        slices (dict): Stores edges to slice "x" and "y" axes. 
        root_out (ROOT.TFile): File storing all histograms that are made.
    '''
    def __init__(self,ledger,twoD,fittag,loadExisting=False):
        '''Constructor.

        Args:
            twoD (TwoDAlphabet): Object with meta information about the run.
            fittag (str): Either 's' or 'b'.
            loadExisting (bool, optional): Flag to load existing projections instead of remaking everything. Defaults to False.
        '''

        self.fittag = fittag
        self.twoD = twoD
        self.ledger = ledger
        self.yaxis1D_title = 'Events / bin'
        self.df = pandas.DataFrame(columns=['process','region','process_type','title'])
        self.dir = 'plots_fit_{f}'.format(f=self.fittag)
        self.slices = {'x': {}, 'y': {}}
        self.root_out = None
        self.df_path = None

        if not loadExisting:
            self._make()
        else:
            self._load()

    def __del__(self):
        '''On deletion, save DataFrame to csv and close ROOT file.'''
        if self.df_path is not None:
            self.df.to_csv(self.df_path)
        if self.root_out is not None and hasattr(self.root_out, 'Close'):
            self.root_out.Close()

    def _load(self):
        '''Open pickled DataFrame and output ROOT file
        and reference with `self.df` and `self.root_out` attributes.'''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name)
        
        df_csv_path = os.path.join(self.dir, 'df.csv')
        self.df = pandas.read_csv(df_csv_path)
        self.df_path = df_csv_path

    def _format_1Dhist(self, hslice, title, xtitle, ytitle, color, proc_type):
        '''Perform some basic formatting of a 1D histogram so that the ROOT.TH1
        already has some of the meta information set like line and fill colors.
        Also renormalizes bins to the minimum bin width if the y-axis is requested
        to be plotted as events per unit rather than events per bin.

        Args:
            hslice (ROOT.TH1): Histogram.
            title (str): Title for the output histogram.
            xtitle (str): X-axis title for the output histogram.
            ytitle (str): Y-axis title for the output histogram. Will be modified if requesting to plot events/unit.
            color (int): ROOT color code.
            proc_type (str): Process type. Either 'BKG', 'SIGNAL', or 'DATA'.

        Raises:
            NameError: If proc_type is not 'BKG', 'SIGNAL', or 'DATA'.

        Returns:
            TH1: Formatted histogram.
        '''
        ytitle = self.yaxis1D_title
        if self.twoD.options.plotEvtsPerUnit:
            hslice = convert_to_events_per_unit(hslice)
            ytitle = 'Events / %s GeV' % get_min_bin_width(hslice)
        hslice.SetMinimum(0)
        hslice.SetTitle(title)
        hslice.GetXaxis().SetTitle(xtitle)
        hslice.GetYaxis().SetTitle(ytitle)

        # Obtain the ROOT color code from the dictionary defined in TwoDAlphabet.plotstyle
        if color not in mpl_to_root_colors.keys():
            available_colors = '", "'.join(mpl_to_root_colors.keys())
            raise ValueError(f'Color "{color}" not defined. Please add the ROOT TColor code to the "mpl_to_root_colors" dictionary defined in TwoDAlphabet.plotstyle. Available default colors are: "{available_colors}"')
        else:
            color = int(mpl_to_root_colors[color]) #ROOT call in C++ sometimes cannot convert it to int

        if proc_type == 'BKG':
            hslice.SetFillColor(color)
            hslice.SetLineColorAlpha(0,0)
        elif proc_type == 'SIGNAL' or proc_type == 'TOTAL':
            hslice.SetLineColor(color)
        elif proc_type == 'DATA':
            hslice.SetLineColor(color)
            hslice.SetMarkerColor(color)
        else:
            raise NameError('Process type "%s" is not supported.'%proc_type)

        return hslice

    def _make(self):
        '''Make the DataFrame and output ROOT file from scratch
        and reference with `self.df` and `self.root_out` attributes.
        
        Loops over all regions and processes from the pre-fit and post-fit shapes
        and tracks/constructs the 2D histograms and six projections (three each for "x" and "y"). 
        '''
        root_out_name = '%s/all_plots.root'%self.dir
        self.root_out = ROOT.TFile.Open(root_out_name,'RECREATE')

        shapes_file = ROOT.TFile.Open('postfitshapes_%s.root'%self.fittag)
        loc_base = '{r}_{c}_{t}/{p}'

        proc_reg_pairs = self.ledger.GetProcRegPairs()+[('TotalBkg', r) for r in self.ledger.GetRegions()]

        for region in self.ledger.GetRegions():
            binning,_ = self.twoD.GetBinningFor(region)
            if len(self.twoD.options.blindedPlots) > 0:
                if region in self.twoD.options.blindedPlots:
                    blinding = [1,2]
                else: blinding = []
            else: blinding = []

            self.slices['x'][region] = {'vals': binning.xSlices,'idxs':binning.xSliceIdx}
            self.slices['y'][region] = {'vals': binning.ySlices,'idxs':binning.ySliceIdx}
            
            for process in self.ledger.GetProcesses()+['TotalBkg']:
                # Skip processes not in this region
                if process not in [pair[0] for pair in proc_reg_pairs if pair[1] == region]:
                    continue
                
                if process != 'TotalBkg':
                    color = self.ledger.GetProcessColor(process)
                    proc_type = self.ledger.GetProcessType(process)
                    proc_title = self.ledger.GetProcessTitle(process)
                else:
                    #color = ROOT.kBlack
                    color = 'black'
                    proc_type = 'TOTAL'
                    proc_title = 'TotalBkg'

                self.df = pandas.concat([self.df,pandas.DataFrame([{'process':process,
                                          'region':region,
                                          'process_type': proc_type,
                                          'title': proc_title}])], ignore_index=True)

                for time in ['prefit','postfit']:
                    # 2D distributions first
                    out2d_name = '%s_%s_%s_2D'%(process,region,time)

                    low_name = loc_base.format(r=region, c='LOW', t=time, p=process)
                    sig_name = loc_base.format(r=region, c='SIG', t=time, p=process)
                    high_name = loc_base.format(r=region, c='HIGH',t=time, p=process)
                    
                    low  = shapes_file.Get(low_name)
                    sig  = shapes_file.Get(sig_name)
                    high = shapes_file.Get(high_name)

                    if low == None: raise IOError('Could not find histogram %s in postfitshapes_%s.root'%(low_name, self.fittag))
                    if sig == None: raise IOError('Could not find histogram %s in postfitshapes_%s.root'%(sig_name, self.fittag))
                    if high == None: raise IOError('Could not find histogram %s in postfitshapes_%s.root'%(high_name, self.fittag))

                    full = stitch_hists_in_x(out2d_name, binning, [low,sig,high], blinded=blinding if process == 'data_obs' else [])
                    full.SetMinimum(0)
                    full.SetTitle('%s, %s, %s'%(proc_title,region,time))

                    self.root_out.WriteTObject(full,full.GetName())

                    # Now do projections using the 2D
                    out_proj_name = '{p}_{r}_{t}_proj{x}{i}'
                    #for proj in ['X','Y']:
                    for proj in ['X']:
                        slices = self.slices['x' if proj == 'Y' else 'y'][region]

                        for islice in range(3):
                            hname = out_proj_name.format(p=process,r=region,t=time,x=proj.lower(),i=islice)
                            start,stop = _get_start_stop(islice,slices['idxs'])
                            
                            hslice = getattr(full,'Projection'+proj)(hname,start,stop,'e')
                            hslice_title = '%s, %s, %s, %s-%s'%(proc_title,region,time,slices['vals'][islice],slices['vals'][islice+1])
                            hslice = self._format_1Dhist(
                                hslice, hslice_title,
                                binning.xtitle if proj == 'X' else binning.ytitle,
                                self.yaxis1D_title,
                                color, proc_type)

                            self.root_out.WriteTObject(hslice,hslice.GetName())

        shapes_file.Close()
        self.root_out.Close()
        self.root_out = ROOT.TFile.Open(root_out_name)

    def Get(self,hname=None,row=None,hist_type=None):
        '''Get a histogram by name from the master ROOT file.
        Does quick check for if the histogram is saved.
        
        Args:
            hname (str): Histogram name.

        Raises:
            LookupError: If histogram cannot be found.
        '''
        if hname != None:
            if hname not in [k.GetName() for k in self.root_out.GetListOfKeys()]:
                raise LookupError('Histogram %s not found in %s'%(hname,self.root_out.GetName()))
            name = hname
        else:
            name = '_'.join([row.process,row.region,hist_type])
        return self.root_out.Get(name)

    def _order_df_on_proc_list(self,df,proc_type,proclist=[],alphaBottom=True):
        '''Re-order input dataframe (`df`) based on the ordered list of process names (`proclist`).
        Useful for pre-ordering process before trying to construct a THStack.

        Args:
            df (pandas.DataFrame): Input DataFrame to manipulate.
            proc_type (str): Process type.  Either 'BKG', 'SIGNAL', or 'DATA'.
            proclist (list, optional): Ordered list of processes to order by. Defaults to [] in which case order is determined based on alphaBottom.
            alphaBottom (bool, optional): Only matters if proclist == []. Defaults to True in which case parametric Alphabet objects included first (thus, on the bottom of the THStack).

        Returns:
            pandas.DataFrame: Ordered DataFrame.
        '''
        if proclist == []:
            if alphaBottom:
                process_order = self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=True, includeConfig=False) + self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=False, includeConfig=True)
            else:
                process_order = self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=False, includeConfig=True) + self.ledger.GetProcesses(ptype=proc_type, includeNonConfig=True, includeConfig=False)
        else:
            process_order = proclist

        process_order_df = pandas.DataFrame({'process':process_order})
        return process_order_df.merge(df,on='process',how='inner')

    def plot_2D_distributions(self):
        '''Take the saved 2D distributions and plot them together on sub-pads
        based on process and region groupings.

        Plots are grouped based on process and then the regions and pre-fit/post-fit
        plots share the same canvas as sub-pads.

        Returns:
            None
        '''
        for pr, _ in self.df.groupby(['process','region']):
            process, region = pr[0], pr[1]
            out_file_name = '{d}/base_figs/{p}_{r}_%s_2D'.format(d=self.dir,p=process,r=region)
            make_pad_2D(outname=out_file_name%('prefit'), hist=self.Get('{p}_{r}_{t}'.format(p=process,r=region,t='prefit_2D')),
                            year=self.twoD.options.year, savePDF=True, savePNG=True)
            make_pad_2D(outname=out_file_name%('postfit'), hist=self.Get('{p}_{r}_{t}'.format(p=process,r=region,t='postfit_2D')),
                            year=self.twoD.options.year, savePDF=True, savePNG=True)

            make_can('{d}/{p}_{r}_2D'.format(d=self.dir,p=process,r=region), [out_file_name%('prefit')+'.png', out_file_name%('postfit')+'.png'])

    def plot_projections(self, lumiText=r'138 $fb^{-1}$ (13 TeV)', extraText='Internal', pf_slice_str={}, subtitles={}, units='GeV', regionsToGroup=[]):
        '''Plot comparisons of data and the post-fit background model and signal
        using the 1D projections. Canvases are grouped based on projection axis.
        The canvas rows are separate selection regions while the columns
        are the different slices of the un-plotted axis.

        Args:
            lumiText (string): LaTeX-formatted string containing luminosity information. Defaults to Run2 conditions.
            extraText (str): Additional text to place after experiment (CMS) text.
            pf_slice_str ({str:str}): Dict of raw strings corresponding to the pass/fail regions specified in the JSON to be placed as the axis slice text in top left corner of pad.
            subtitles ({str:str}, optional): Dict of raw strings corresponding to each region specified in the JSON to be placed underneath axis slice text in top left corner of pad. If multiple titles are desired, separate with semicolon character.
                Example: {"SR_fail": r"$ParticleNet_{TvsQCD}$ Pass\n$ParticleNetMD_{Xbb}$ Fail", "SR_pass": r"$ParticleNet_{TvsQCD}$ Pass;$ParticleNetMD_{Xbb}$ Pass"}
            units (str): Units of measurement for observable. Placed on x-axis and in slice string.

            regionsToGroup ([[str]]): 
                List of list of strings representing the desired regions to group. For example if the fit involved
                four regions: CR_fail, CR_pass, SR_fail, SR_pass then 2DAlphabet will try to plot all
                (4 regions) x (3 slices) = 12 plots on the same page, which is greater than the number that can be
                plotted. Instead, pass regionsToGroup = [['CR'],['SR']] to have the CR and SR plotted on separate
                2x3 canvases. 
                If you wanted to plot, e.g. SR_fail, SR_pass, and CR_pass, pass in the following list of lists:
                    [['CR'], ['SR'], ['SR_fail','SR_pass','CR_pass']]
                The order matters if the sub-list contains multiple strings - the regions are plotted in order
                with the first on top and last on bottom of the canvas.

        Returns:
            None
        '''
        axes = pandas.DataFrame() # Book a dataframe for creating a full group of plots
        for region, group in self.df.groupby('region'):
            binning, _ = self.twoD.GetBinningFor(region)

            # Make both regular and logarithmic y-axis plots
            for logyFlag in [False, True]:
                # Get reduced dataframes for all backgrounds and the signals
                ordered_bkgs = self._order_df_on_proc_list(
                                            group[group.process_type.eq('BKG')], proc_type='BKG',
                                            alphaBottom=(not logyFlag))
                signals = group[group.process_type.eq('SIGNAL')]

                for proj in ['postfit_projx']:
                #for proj in ['postfit_projx','postfit_projy']:
                    for islice in range(3):
                        projn     = f'{proj}{islice}'
                        if projn != 'postfit_projx2': continue
                        sig_projn = projn
                        if self.twoD.options.plotPrefitSigInFitB and self.fittag == 'b':
                            sig_projn = projn.replace('postfit','prefit') # Plot prefit signal in b-only plots

                        # TH1s representing the sum of all data and background histograms in the workspace for this projection and slice
                        this_data =      self.Get(row=group.loc[group.process_type.eq('DATA')].iloc[0], hist_type=projn)
                        this_totalbkg =  self.Get(row=group.loc[group.process_type.eq('TOTAL')].iloc[0], hist_type=projn)
                        # Lists of individual TH1s comprising all bkg and sig histograms for this projection and slice.
                        these_bkgs =    [self.Get(row=ordered_bkgs.iloc[irow], hist_type=projn) for irow in range(ordered_bkgs.shape[0])]
                        these_signals = [self.Get(row=signals.iloc[irow], hist_type=sig_projn) for irow in range(signals.shape[0])]

                        slice_edges = (
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice],
                            binning.xtitle if 'y' in proj else binning.ytitle,
                            self.slices['x' if 'y' in proj else 'y'][region]['vals'][islice+1],
                        )
                        slice_str = '%s < %s < %s'%slice_edges

                        if (region == "pass") :
                            slice_str = pf_slice_str[region]
                        elif (region == "fail") :
                            slice_str = pf_slice_str[region]

                        out_pad_name = f'{self.dir}/base_figs/{projn}_{region}{"" if not logyFlag else "_logy"}'

                        # If user requested sub-titles, obtain the ones for this specific region
                        if subtitles:
                            subtitle = subtitles[region]
                        else:
                            subtitle = {}

                        # Produce a matplotlib axis containing the full data + stacked bkg histogram for this projection and slice
                        make_ax_1D(
                            out_pad_name, 
                            binning, 
                            data=this_data, 
                            bkgs=these_bkgs, 
                            signals=these_signals, 
                            totalBkg=this_totalbkg, 
                            logyFlag=logyFlag, 
                            slicetitle=slice_str,   # opposite axis slices (top left, in pad)
                            subtitle=subtitle,       # (top left, under slice string)
                            lumiText=lumiText,      # (top right, above pad)
                            extraText=extraText,    # top left, above pad, next to "CMS"
                            units=units,
                            savePDF=True,
                            savePNG=True
                        )

                        # Append projection plot to master dataframe
                        axes = pandas.concat([axes, pandas.DataFrame([{'ax':out_pad_name+'.png', 'region':region, 'proj':projn, 'logy':logyFlag}])], ignore_index=True)

        # Create a canvas with the full set of projection plots
        for logy in ['', '_logy']:
            for proj in ['postfit_projx']:
            #for proj in ['postfit_projx','postfit_projy']:
                these_axes = axes.loc[axes.proj.str.contains(proj)]
                if logy == '':
                    these_axes = these_axes.loc[these_axes.logy.eq(False)]
                else:
                    these_axes = these_axes.loc[these_axes.logy.eq(True)]

                if (len(these_axes) > 9) and (len(regionsToGroup) == 0):
                    raise RuntimeError('histlist of size %s not currently supported. Instead, call plot_projections() with regionsToGroup list describing the regions you want to group together.'%len(these_axes))
                elif (len(these_axes) > 9) and (len(regionsToGroup) > 0):
                    validRegions = these_axes['region'].to_list()
                    for i in range(len(regionsToGroup)):
                        region = regionsToGroup[i]
                        if (len(region) > 1):	# e.g. ['SR_fail', 'SR_pass', 'ttbarCR_pass']
                            # Ensure that the regions are plotted in order
                            new_axes = []
                            for r in region:
                                if r not in validRegions:
                                    raise ValueError(f'Region "{r}" specified in regionsToGroup is not available in the 2DAlphabet workspace. Available regions:\n\t{validRegions}')
                                r_axes = these_axes[these_axes['region'].str.match(r)].sort_values(by=['region','proj'])['ax'].to_list()
                                new_axes += r_axes
                            '''
                            new_axes = []
                            # Ensure that the regions are plotted in order provided
                            for r in region:
                                for ax in these_axes:
                                    # As an edge case, assume r would be: CR_pass
                                    # pad would be: plots_fit_b/base_figs/postfit_projx2_CR_pass.png
                                    # but also: plots_fit_b/base_figs/postfit_projx2_ttbarCR_pass.png
                                    # so we have to append an underscore to ensure we get the right one 
                                    rNew = '_'+r
                                    if rNew in ax:
                                        new_axes.append(ax)
                            '''
                            out_can_name = '{d}/{reg}_{proj}{logy}'.format(d=self.dir,proj=proj,logy=logy,reg='_and_'.join(region))
                            make_can(out_can_name, new_axes)
                        else:	# e.g. ['SR']
                            new_axes = these_axes[these_axes['region'].str.match(f'{region[0]}_')].sort_values(by=['region','proj'])['ax'].to_list()
                            out_can_name = '{d}/{reg}_{proj}{logy}'.format(d=self.dir,proj=proj,logy=logy,reg=region[0])
                            make_can(out_can_name, new_axes)
                else:
                    these_axes = these_axes.sort_values(by=['region','proj'])['ax'].to_list()
                    out_can_name = '{d}/{proj}{logy}'.format(d=self.dir, proj=proj, logy=logy)
                    make_can(out_can_name, these_axes)


def _save_pad_generic(outname, pad, ROOTout, savePDF, savePNG):
    if isinstance(ROOTout, ROOT.TFile):
        ROOTout.WriteTObject(pad,outname.split('/')[-1])
    if savePDF:
        pad.Print(outname+'.pdf','pdf')
    if savePNG:
        pad.Print(outname+'.png','png')

def _make_pad_gen(name):
    tdrstyle.setTDRStyle()
    ROOT.gStyle.SetLegendFont(42)
    ROOT.gStyle.SetTitleBorderSize(0)
    ROOT.gStyle.SetTitleAlign(33)
    ROOT.gStyle.SetTitleX(.77)

    pad = ROOT.TCanvas(name, name, 800, 700)
    pad.cd(); pad.SetRightMargin(0.0); pad.SetTopMargin(0.0); pad.SetBottomMargin(0.0)
    return pad

def make_pad_2D(outname, hist, style='lego', logzFlag=False, ROOTout=None,
                savePDF=False, savePNG=False, year=1, extraText='Internal'):
    '''Make a pad holding a 2D plot with standardized formatting conventions.

    Args:
        outname (str): Output file path name.
        hist (TH2): Histogram to draw on the pad.
        style (str, optional): ROOT drawing style. Defaults to 'lego'.
        logzFlag (bool, optional): Make log z-axis. Defaults to False.
        saveROOT (bool, optional): Save to master ROOT file. Defaults to True.
        savePDF (bool, optional): Save to PDF. Defaults to False.
        savePNG (bool, optional): Save to PNG. Defaults to False.
        year (int, optional): Luminosity formatting. Options are 16, 17, 18, 1 (full Run 2), 2 (16+17+18). Defaults to 1.
        extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.

    Returns:
        [type]: [description]
    '''
    pad = _make_pad_gen(outname)
    pad.SetLeftMargin(0.15)
    pad.SetRightMargin(0.2)
    pad.SetBottomMargin(0.12)
    pad.SetTopMargin(0.1)
    if logzFlag: pad.SetLogz()

    hist.GetXaxis().SetTitleOffset(1.15); hist.GetXaxis().SetLabelSize(0.05); hist.GetXaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetTitleOffset(1.5);  hist.GetYaxis().SetLabelSize(0.05); hist.GetYaxis().SetTitleSize(0.05)
    hist.GetZaxis().SetTitleOffset(1.5);  hist.GetZaxis().SetLabelSize(0.05); hist.GetZaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetNdivisions(505)
    
    if 'lego' in style.lower():
        hist.GetZaxis().SetTitleOffset(1.4)

    hist.Draw(style)
    
    CMS_lumi.extraText = extraText
    CMS_lumi.CMS_lumi(pad, year, 11, sim=False if 'data' in hist.GetName().lower() else True)

    ROOT.SetOwnership(pad, False)
    _save_pad_generic(outname, pad, ROOTout, savePDF, savePNG)

    return pad


def make_ax_1D(outname, binning, data, bkgs=[], signals=[], title='', subtitle='', slicetitle='',
            totalBkg=None, logyFlag=True, ROOTout=None, savePDF=False, savePNG=False,
            dataOff=False, datastyle='pe X0', year=1, addSignals=True, 
            lumiText=r'$138 fb^{-1} (13 TeV)$', extraText='Internal', units='GeV', hspace=0.0):
    '''Create a matplotlib.axis.Axis object holding a 1D plot with standardized CMS formatting conventions
    Args:
        outname (str): Output file path + name.
        binning (Binning): TwoDAlphabet binning object for the given region from which to obtain binning info.
        data (TH1): Data histogram.
        bkgs ([TH1]): List of background histograms (to be stacked).
        signals ([TH1]): List of signal histograms.
        title (str, optional): Title of plot. Only applicable if bkgs is empty. Defaults to ''.
        subtitle (str, optional): Subtitle text for physics information (like slice ranges). Defaults to ''.
        totalBkg (TH1, optional): Total background estimate from fit. Used to get total background uncertianty. Defaults to None.
        logyFlag (bool, optional): Make log y-axis. Defaults to False (now to True //TAV).
        savePDF (bool, optional): Save to PDF. Defaults to True.
        savePNG (bool, optional): Save to PNG. Defaults to True.
        dataOff (bool, optional): Turn off the data from plotting. Defaults to False.
        extraText (str, optional): Prepended to the CMS subtext. Defaults to 'Preliminary'.
        units (str, optional): Units of measurement for x- and y-axes. Defaults to 'GeV'.
        hspace (float, optional): Spacing between the main plotting axis and the ratio subplot axis. Defaults to 0.0 (no spacing). Original 2DAlphabet spacing would correspond to hspace=0.7.
    Returns:
        ax (matplotlib.axis.Axis): Output axis object for further manipulation.
    '''
    # Convert all histograms to numpy arrays. First, determine the bin edges
    projn  = outname.split('/')[-1].split('_')[1].split('proj')[-1][0]
    islice = outname.split('/')[-1].split('_')[1].split('proj')[-1][1]; islice = int(islice)
    if projn == 'x':
        xbins = binning.xbinByCat
        edges = np.array(xbins['LOW'][:-1]+xbins['SIG'][:-1]+xbins['HIGH'])
    else:
        edges = np.array(binning.ybinList)
    widths = np.diff(edges)     # obtain bin widths

    # Convert all unique bkg and signal hists to arrays and group by process. They will come already ordered, so keep the ordering
    bkgDict = OrderedDict()
    sigDict = OrderedDict()
    bkgNames = list(dict.fromkeys([hist.GetTitle().split(',')[0] for hist in bkgs]))
    sigNames = list(dict.fromkeys([hist.GetTitle().split(',')[0] for hist in signals]))
    # Replace the ROOT latex "#" with standard latex "\" escape character for python rstring
    bkgNamesLatex = [r'${}$'.format(re.sub(r'^CMS_[A-Z]+\d*_', '', bkgName.replace("#","\\"))) for bkgName in bkgNames] # strip the CMS_<CADI>_ prefix
    sigNamesLatex = [r'${}$'.format(sigName.replace("#","\\")) for sigName in sigNames]
    # Sum the common backgrounds and signals
    for bkg in bkgNames:
        bkg_arrs = [hist2array(i, return_errors=True)[0] for i in bkgs if bkg in i.GetTitle()]
        bkg_errs = [hist2array(i, return_errors=True)[1] for i in bkgs if bkg in i.GetTitle()]
        bkgDict[f'{bkg}_arr'] = sum(bkg_arrs)
        bkgDict[f'{bkg}_err'] = sum(bkg_errs)
    for sig in sigNames:
        sigDict[f'{sig}_arr'] = sum([hist2array(i) for i in signals if sig in i.GetTitle()]) # ignore uncertainty
    # For the colors, we need to translate from ROOT TColor to matplotlib 
    bkgColors = list(dict.fromkeys([hist.GetFillColor() for hist in bkgs]))
    bkgColors = [root_to_matplotlib_color(TColor) for TColor in bkgColors]
    sigColors = list(dict.fromkeys([hist.GetLineColor() for hist in signals]))  # signals have no fill, only line color
    sigColors = [root_to_matplotlib_color(TColor) for TColor in sigColors]

    # Get the data array and total bkg
    data_arr = hist2array(data);
    for i, val in enumerate(data_arr):
        if np.isnan(val):
            data_arr[i] = 0.0
    data_arr = np.array([int(i) for i in data_arr]) # hist2array converts to floats which may leave very small differences when converting TH1 (int) -> array (float). This fixes it
    if totalBkg:
        totalBkg_arr, totalBkg_err = hist2array(totalBkg, return_errors=True)

    # Begin plotting
    plt.style.use([hep.style.CMS])
    fig, (ax, rax) = plt.subplots(2, 1, sharex=True, **ratio_fig_style) # ax is the stacked plot, rax contains the ratio plot
    fig.subplots_adjust(hspace=hspace)

    bkg_stack = np.vstack([arr for key, arr in bkgDict.items() if '_arr' in key]) # stack all unique background processes
    # depending on step option ('pre' or 'post'), the last bin needs be concatenated on one side, so that the edge bin is drawn (annoying)
    bkg_stack = np.hstack([bkg_stack, bkg_stack[:,-1:]])
    bkg_stack = np.hstack([bkg_stack])

    ax.stackplot(edges, bkg_stack, labels=bkgNamesLatex, colors=bkgColors, step='post', **stack_style)
    unc = totalBkg_err # the hist2array() call returns the sqrt of the sumw2 for the histogram
    unc = np.hstack([unc, unc[-1]]) # concatenation of last bin so that edge bin is drawn
    ax.fill_between(x=edges, y1=bkg_stack.sum(axis=0)-unc, y2=bkg_stack.sum(axis=0)+unc, label='Bkg. Unc.', step='post', **hatch_style)
    bin_centers = (edges[:-1] + edges[1:])/2
    # Determine if we have variable binning, and if so, add the CMS-required horizontal bars to denote bin width
    if len(np.unique(widths)) > 1: # Detected variable bin widths
        xerrs = (edges[1:]-edges[:-1])/2
        ax.set_ylabel('Events / bin')
    else:
        ax.set_ylabel(f'Events / {widths[0]} {units}')
        xerrs = None

    #Data plot with Garwood Poisson CI as error bars
    lower_errors, upper_errors = poisson_conf_interval(data_arr)
    yerr = [data_arr - lower_errors, upper_errors - data_arr]
    
    ax.errorbar(x=bin_centers, y=data_arr, yerr=yerr, xerr=xerrs, label='Data', **errorbar_style)
    #ax.errorbar(x=bin_centers, y=data_arr, yerr=np.sqrt(data_arr), xerr=xerrs, label='Data', **errorbar_style)

    # Plot signals
    for i, sig in enumerate(sigNames):
        sigarr = sigDict[f'{sig}_arr']
        ax.step(x=edges, y=np.hstack([sigarr, sigarr[-1]]), where='post', color=sigColors[i], label=sigNamesLatex[i])

    if logyFlag:
        ax.set_ylim(0.0001, totalBkg_arr.max()*1e5)
        ax.set_yscale('log')
    else:
        ax.set_ylim(0, totalBkg_arr.max()*1.38)

    # Make sure data and signal(s) come first 
    handles, labels = ax.get_legend_handles_labels()
    data_idx = labels.index('Data') # should always be last entry in legend
    unc_idx  = labels.index('Bkg. Unc.')
    sig_idxs = [i for i in range(unc_idx+1, data_idx)] # indices of all signals in legend
    # Resulting legend ordering will be Data, signal(s), bkgs, then bkg unc. The bkg ordering is already sent to this function properly.
    leg_order = [data_idx] + sig_idxs + [idx for idx in range(len(labels)) if idx not in [data_idx] + sig_idxs]
    ax.legend([handles[idx] for idx in leg_order],[labels[idx] for idx in leg_order])
    ax.autoscale(axis='x', tight=True)
    ax.margins(x=0) # remove white space at left and right margins of plot 

    hep.cms.label(loc=0, ax=ax, data = not dataOff, label=extraText, rlabel='') # CMS + label, where label is typically “Preliminary” “Supplementary”, “Private Work” or “Work in Progress”
    hep.cms.lumitext(lumiText, ax=ax)                       # Typically luminosity + sqrt(s)
    # Can't use the hep.cms.text() wrapper without "CMS" being added, so add the slice text manually
    if slicetitle:
        slicetitle = r'${}$'.format(slicetitle.replace('#','\\'))
        ax.text(0.3, 0.95, slicetitle, ha='center', va='top', fontsize='small', transform=ax.transAxes)
    if subtitle:
        # Check if user requested multi-line
        if (len(subtitle.split(';')) == 1):
            ax.text(0.3, 0.89, r'%s'%subtitle, ha='center', va='top', fontsize='small', transform=ax.transAxes)
        else:
            for i, title in enumerate(subtitle.split(';')):
                ax.text(0.3, 0.95-(0.06*(i+1)), r'%s'%title, ha='center', va='top', fontsize='small', transform=ax.transAxes)

    # pull
    # TAV: After unblinding one needs to come back to this
    #dataMinusBkg = data_arr - totalBkg_arr
    dataMinusBkg = np.where(data_arr == 0, 0, data_arr - totalBkg_arr)
    data_error = np.where(dataMinusBkg<0, upper_errors - data_arr, data_arr - lower_errors)
    sigmas = np.sqrt(data_error**2 + totalBkg_err**2) 
    #sigmas = np.sqrt(np.sqrt(data_arr)*np.sqrt(data_arr) + totalBkg_err*totalBkg_err)
    sigmas[sigmas==0.0] = 1e-5 # avoid division by zero 
    pulls =  dataMinusBkg/sigmas
    rax.bar(bin_centers,pulls, width=widths, color='gray')
    rax.axhline(0.0, color='black', linestyle='--', linewidth=0.8, zorder=3) # reference line at 0 pull
    rax.set_ylim(-3,3)
    rax.set_ylabel(r'$\frac{Data-Bkg}{\sigma}$')
    axisTitle = binning.xtitle if projn == 'x' else binning.ytitle
    axisTitle = axisTitle.replace("#","\\")
    rax.set_xlabel(r'${}$ [{}]'.format(axisTitle, units) if units else r'${}$'.format(axisTitle))  # no brackets without units
    rax.autoscale(axis='x', tight=True)
    rax.margins(x=0)

    if savePDF:
        plt.savefig(f'{outname}.pdf')
    if savePNG:
        plt.savefig(f'{outname}.png')
    if ((not savePDF) and (not savePNG)):
        print(f'WARNING: plot "{outname}" has not been saved.')

    print(f'Plotting {outname}')
    plt.close()
    # return ax, rax      ????????????

def make_can(outname, padnames, padx=0, pady=0):
    '''Combine multiple pads/canvases into one canvas for convenience of viewing.
    Input pad order matters.

    Args:
        outname (str): Output file path name.
        pads ([TCanvas,TPad]): List of canvases/pads to plot together on one canvas.
        saveROOT (bool, optional): Save to master ROOT file. Defaults to False.
        savePDF (bool, optional): Save to PDF. Defaults to True.
        savePNG (bool, optional): Save to PNG. Defaults to True.

    Raises:
        RuntimeError: If 10 or more subdivisions are requested.

    Returns:
        None
    '''
    if padx == 0 or pady == 0:
        if len(padnames) == 1:
            padx = 1; pady = 1
        elif len(padnames) == 2:
            padx = 2; pady = 1
        elif len(padnames) == 3:
            padx = 3; pady = 1
        elif len(padnames) == 4:
            padx = 2; pady = 2
        elif len(padnames) <= 6:
            padx = 3; pady = 2
        elif len(padnames) <= 9:
            padx = 3; pady = 3
        else:
            raise RuntimeError('histlist of size %s not currently supported'%len(padnames))
            #raise RuntimeError('histlist of size %s not currently supported: %s'%(len(padnames),[p.GetName() for p in padnames]))

    pads = [Image.open(os.path.abspath(pname)) for pname in padnames]
    w, h = pads[0].size
    grid = Image.new('RGB', size=(padx*w, pady*h))
    
    for i, pad in enumerate(pads):
        grid.paste(pad, box=(i%padx*w, i//padx*h))
    
    print ('Writing grid of images %s.pdf'%outname)
    grid.save(outname+'.pdf')
    grid.save(outname+'.png')

def _get_start_stop(i,slice_idxs):
    start = slice_idxs[i]+1
    stop  = slice_idxs[i+1]
    return start, stop

def gen_projections(ledger, twoD, fittag, loadExisting=False, lumiText=r'138 $fb^{-1}$ (13 TeV)', extraText='Internal', pf_slice_str={}, subtitles={}, units='GeV', regionsToGroup=[]):
    plotter = Plotter(ledger, twoD, fittag, loadExisting)
    plotter.plot_2D_distributions()
    plotter.plot_projections(lumiText, extraText, pf_slice_str, subtitles, units, regionsToGroup)

def make_systematic_plots(twoD):
    '''Make plots of the systematic shape variations of each process based on those
    processes and systematic shapes specified in the config. Shapes are presented
    as projections onto 1D axis where no selection has been made on the axis not
    being plotted. Plots are saved to UncertPlots/.
    '''
    for (p,r), _ in twoD.df.groupby(['process','region']):
        if p == 'data_obs': continue

        nominal_full = twoD.organizedHists.Get(process=p, region=r, systematic='')
        binning, _ = twoD.GetBinningFor(r)

        for axis in ['X','Y']:

            nominal_hist = getattr(nominal_full,'Projection'+axis)('%s_%s_%s_%s'%(p,r,'nom','proj'+axis))
            nominal = hist2array(nominal_hist)

            # Get the bin edges from the 2DAlphabet binning object. Avoid edge duplication in the case of X-axis stitching
            if axis == 'X':
                xbins = binning.xbinByCat
                edges = xbins['LOW'][:-1]+xbins['SIG'][:-1]+xbins['HIGH']
            else:
                edges = binning.ybinList

            # GetShapeSystematics() lists all shape systs, even if the given process is not subject to them.
            check_df = twoD.df.loc[(twoD.df.process == p) & (twoD.df.region == r)] # Get the entries for this proc and region
            proc_vars = check_df.variation.unique() # Get all unique variations for this process

            for s in twoD.ledger.GetShapeSystematics(drop_norms=True):
                if s not in proc_vars: continue

                up_hist = getattr(twoD.organizedHists.Get(process=p,region=r,systematic=s+'Up'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Up','proj'+axis))
                up = hist2array(up_hist)
                down_hist = getattr(twoD.organizedHists.Get(process=p,region=r,systematic=s+'Down'),'Projection'+axis)('%s_%s_%s_%s'%(p,r,s+'Down','proj'+axis))
                down = hist2array(down_hist)

                # Begin plotting
                labels = ['Nominal', r'$+1\sigma$', r'$-1\sigma$']
                colors = ['black', 'red', 'blue']
                styles = ['solid', 'dashed', 'dashed']
                histos = [nominal, up, down]

                plt.style.use([hep.style.CMS])
                fig, ax = plt.subplots(figsize=(12, 8), dpi=200)
                hep.cms.text("WiP",loc=1)
                ax.set_title(f'{p}, {s}', pad=12)
                hep.histplot(histos, edges, stack=False, ax=ax, label=labels, histtype='step', linestyle=styles, color=colors)
                handles, labelsproj = ax.get_legend_handles_labels()
                ax.set_xlabel(getattr(binning,axis.lower()+'title'))
                ax.set_ylabel('Events')
                plt.legend(loc='best')

                outname = f'{twoD.tag}/UncertPlots/Uncertainty_{p}_{r}_{s}_proj{axis}.png'
                print(f'[2DAlphabet.plot] INFO: Plotting histogram\n\t{outname}')

                plt.savefig(outname)
                plt.close() # free up memory

		# make the same plot, but add a ratio plot below in the same canvas
                # Calculate the ratio to nominal
                ratio_up = np.where(nominal == 0, 1, up / nominal) # avoid division by zero
                ratio_down = np.where(nominal == 0, 1, down / nominal) # avoid division by zero
                ratio_up = np.where(np.logical_and(ratio_up == 0, nominal == 0), 1, ratio_up) # avoid zero in ratio plot
                ratio_down = np.where(np.logical_and(ratio_down == 0, nominal == 0), 1, ratio_down) # avoid zero in ratio plot
                # Plot the ratio plot
                fig, (ax, rax) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]}, figsize=(12, 10), dpi=200)
                hep.cms.text("WiP",loc=0,ax=ax) # Add CMS text to the main plot
                ax.set_title(f'{p}, {s}', pad=12, fontsize=18)
                # Plot the nominal, up and down histograms
                hep.histplot(histos, edges, stack=False, ax=ax, label=labels, histtype='step', linestyle=styles, color=colors)
                handles, labelsproj = ax.get_legend_handles_labels()
                ax.set_xlabel("")
                ax.set_ylabel('Events')
                ax.legend([handles[idx] for idx in range(len(labels))], [labelsproj[idx] for idx in range(len(labels))])
                
                # Plot the ratio
                hep.histplot([ratio_up, ratio_down], edges,
                            stack=False, ax=rax,
                            histtype='step', linestyle=['dashed', 'dashed'], color=['red', 'blue'])

                max_ratio = max(np.max(ratio_up), np.max(ratio_down))
                if max_ratio == 1:
                    # If both ratios are 1, set limits to 0.8-1.2
                    rax.set_ylim(0.8, 1.2)
                else:
                    rax.set_ylim(0.5, max_ratio*1.1 if max_ratio > 1 else 1.5) # Set limits for ratio plot, ensure it doesn't go below 0.5 or above 1.2 for visual clarity
                rax.yaxis.set_major_locator(plt.MaxNLocator(5)) # Set major ticks to 5 for the ratio plot
                rax.tick_params(axis='y', which='major', labelsize=15) # Set major tick labels to size 10
                rax.axhline(1, color='black', linestyle='solid', linewidth=1) # Add a horizontal line at y=1
                rax.axhline(1.2, color='grey', linestyle='dashed', linewidth=1) # Add a horizontal line at y=1.2
                rax.axhline(0.8, color='grey', linestyle='dashed', linewidth=1) # Add a horizontal line at y=0.8
                rax.set_ylabel('Ratio')
                rax.set_xlabel(f"${getattr(binning, axis.lower()+'title')}$") # Set x-axis label for ratio plot
                
                outname_ratio = f'{twoD.tag}/UncertPlots/Uncertainty_{p}_{r}_{s}_proj{axis}_ratio.png'
                print(f'[2DAlphabet.plot] INFO: Plotting histogram with ratio\n\t{outname_ratio}')
                plt.savefig(outname_ratio)
                plt.close() # free up memory

def _make_pull_plot(data, bkg, preVsPost=False):
    pull = data.Clone(data.GetName()+"_pull")
    pull.Add(bkg,-1)
    
    sigma = 0.0
    for ibin in range(1,pull.GetNbinsX()+1):
        d = data.GetBinContent(ibin)
        b = bkg.GetBinContent(ibin)
        if d >= b:
            derr = data.GetBinErrorLow(ibin)
            berr = bkg.GetBinErrorUp(ibin)
        elif d < b:
            derr = data.GetBinErrorUp(ibin)
            berr = bkg.GetBinErrorLow(ibin)
        
        if d == 0:
            derr = 1

        sigma = math.sqrt(derr*derr + berr*berr)
        if sigma != 0:
            pull.SetBinContent(ibin, (pull.GetBinContent(ibin))/sigma)
        else:
            pull.SetBinContent(ibin, 0.0 )

    pull.SetFillColor(ROOT.kBlue)
    pull.SetTitle(";"+data.GetXaxis().GetTitle()+";({})/#sigma".format('Post-Pre' if preVsPost else 'Data-Bkg'))
    pull.SetStats(0)

    pull.GetYaxis().SetRangeUser(-2.9,2.9)
    pull.GetYaxis().SetTitleOffset(0.4)                             
    pull.GetYaxis().SetLabelSize(0.13)
    pull.GetYaxis().SetTitleSize(0.12)
    pull.GetYaxis().SetNdivisions(306)

    pull.GetXaxis().SetLabelSize(0.13)
    pull.GetXaxis().SetTitleSize(0.15)
    pull.GetYaxis().SetTitle("({})/#sigma".format('Post-Pre' if preVsPost else 'Data-Bkg'))
    return pull

def _get_good_fit_results(tfile):
    successful_fits = []
    for fittag in ['b','s']:
        if 'fit_'+fittag not in [k.GetName() for k in tfile.GetListOfKeys()]:
            warnings.warn('Unable to find result fit_%s...'%fittag,RuntimeWarning)
        else:
            successful_fits.append(fittag)
    return successful_fits

def nuis_pulls(vtol=0.3, stol=0.1, vtol2=2.0, stol2=0.5, regex='^(?!.*(_bin_|_par))'):
    diffnuis_cmd = f"python $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit/test/diffNuisances.py fitDiagnosticsTest.root --abs -g nuisance_pulls.root --vtol={vtol} --stol={stol} --vtol2={vtol2} --stol2={stol2} --all --regex='{regex}'"
    execute_cmd(diffnuis_cmd)
    # Make a PDF of the nuisance_pulls.root
    if os.path.exists('nuisance_pulls.root'):
        nuis_file = ROOT.TFile.Open('nuisance_pulls.root')
        nuis_can = nuis_file.Get('nuisances')
        nuis_can.Print('nuisance_pulls.pdf','pdf')
        nuis_file.Close()

def save_post_fit_parametric_vals():
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    goodFitTags = _get_good_fit_results(fit_result_file)
    for fittag in goodFitTags:
        coeffs_final = fit_result_file.Get('fit_'+fittag).floatParsFinal()
        all_par_names = [] # get names of anything matching *_par<N>
        for i in range(coeffs_final.getSize()):
            name = coeffs_final.at(i).GetName()
            if name.split('_')[-1].startswith('par'):
                all_par_names.append(name)
        all_par_names.sort()
        # Get unique prefixes to _par<N>
        all_obj_names = set(['_'.join(name.split('_')[:-1]) for name in all_par_names])

        for obj_name in all_obj_names:
            with open('rpf_params_%s_fit%s.txt'%(obj_name,fittag),'w') as param_out:
                for par_name in [p for p in all_par_names if p.startswith(obj_name)]:
                    var = coeffs_final.find(par_name)
                    param_out.write('%s: %s +/- %s\n'%(par_name, var.getValV(), var.getError()))

    fit_result_file.Close()

def gen_post_fit_shapes():
    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    goodFitTags = _get_good_fit_results(fit_result_file)
    for t in goodFitTags:
        if os.path.exists('higgsCombineTest.FitDiagnostics.mH120.root'):
            workspace_file = 'higgsCombineTest.FitDiagnostics.mH120.root'
        else:
            workspace_file = 'higgsCombineTest.FitDiagnostics.mH120.123456.root'
        shapes_cmd = 'PostFit2DShapesFromWorkspace -w {w} --output postfitshapes_{t}.root -f fitDiagnosticsTest.root:fit_{t} --postfit --samples 100 --print > PostFitShapes2D_stderr_{t}.txt'.format(t=t,w=workspace_file)
        execute_cmd(shapes_cmd)
    fit_result_file.Close()

def _reduced_corr_matrix(fit_result, varsToIgnore=[], varsOfInterest=[], threshold=0):
    if threshold < 0:
        raise ValueError('Threshold for correlation matrix values to plot must be a positive number.')

    ROOT.gStyle.SetOptStat(0)
    # ROOT.gStyle.SetPaintTextFormat('.3f')
    CM = fit_result.correlationMatrix()
    finalPars = fit_result.floatParsFinal()

    nParams = CM.GetNcols()
    finalParamsDict = {}
    for cm_index in range(nParams):
        if varsOfInterest == []:
            if finalPars.at(cm_index).GetName() not in varsToIgnore:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index
        else:
            if finalPars.at(cm_index).GetName() in varsOfInterest:
                finalParamsDict[finalPars.at(cm_index).GetName()] = cm_index

    nFinalParams = len(finalParamsDict.keys())
    out = ROOT.TH2D('correlation_matrix','correlation_matrix',nFinalParams,0,nFinalParams,nFinalParams,0,nFinalParams)
    out_txt = ''

    for out_x_index, paramXName in enumerate(sorted(finalParamsDict.keys())):
        cm_index_x = finalParamsDict[paramXName]

        if not any([abs(CM[cm_index_x][cm_index_y]) > threshold for cm_index_y in range(nParams) if cm_index_y != cm_index_x]):
            continue

        for out_y_index, paramYName in enumerate(sorted(finalParamsDict.keys())):
            cm_index_y = finalParamsDict[paramYName]
            if cm_index_x > cm_index_y:
                out_txt += '%s:%s = %s\n'%(paramXName,paramYName,CM[cm_index_x][cm_index_y])
            out.Fill(out_x_index+0.5,out_y_index+0.5,CM[cm_index_x][cm_index_y])

        out.GetXaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
        out.GetYaxis().SetBinLabel(out_x_index+1,finalPars.at(cm_index_x).GetName())
    out.SetMinimum(-1)
    out.SetMaximum(+1)

    return out, out_txt

def plot_correlation_matrix(varsToIgnore, threshold=0, corrText=False, cmsText='Internal', lumiText=''):
    '''Plot the post-fit correlation matrix for each successful fit (b-only, s+b).

    Reads fitDiagnosticsTest.root in the current directory and writes a CMS-styled
    correlation_matrix.{png,pdf,txt} into plots_fit_<fittag>/. The bin-by-bin nuisances
    are usually passed via varsToIgnore so only the TF parameters, the POI, and the
    named systematics remain.

    Args:
        varsToIgnore (list): Parameter names to drop from the matrix.
        threshold (float, optional): Skip rows whose off-diagonal correlations are all
            below this absolute value. Defaults to 0.
        corrText (bool, optional): Overlay the numerical correlation values. Defaults to False.
        cmsText (str, optional): Sublabel drawn after 'CMS'. Defaults to 'Internal'.
        lumiText (str, optional): Text drawn at the top right (e.g. luminosity/energy).
            Defaults to '' (nothing drawn).
    '''
    from array import array

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)
    ROOT.gStyle.SetOptTitle(False)
    # Symmetric blue-white-red diverging palette so zero correlation maps to white
    stops = array('d', [0.00, 0.50, 1.00])
    red   = array('d', [0.00, 1.00, 1.00])
    green = array('d', [0.00, 1.00, 0.00])
    blue  = array('d', [1.00, 1.00, 0.00])
    ROOT.TColor.CreateGradientColorTable(3, stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)
    ROOT.gStyle.SetPaintTextFormat('.2f')

    fit_result_file = ROOT.TFile.Open('fitDiagnosticsTest.root')
    for fittag in _get_good_fit_results(fit_result_file):
        fit_result = fit_result_file.Get("fit_"+fittag)
        if not hasattr(fit_result,'correlationMatrix'):
            warnings.warn('Not able to produce correlation matrix.',RuntimeWarning)
            continue

        corrMtrx, corrTxt = _reduced_corr_matrix(fit_result, varsToIgnore=varsToIgnore, threshold=threshold)

        # Shorten the axis labels for readability (drop the long common prefixes)
        for axis in (corrMtrx.GetXaxis(), corrMtrx.GetYaxis()):
            for b in range(1, axis.GetNbins()+1):
                lab = axis.GetBinLabel(b)
                lab = re.sub(r'^CMS_[A-Z]+\d*_(Background_)?', '', lab)  # strip the CMS_<CADI>_ prefix
                axis.SetBinLabel(b, lab)

        nbins = corrMtrx.GetXaxis().GetNbins()
        # Scale label size with the number of bins so labels stay legible but don't overlap
        labSize = min(0.03, 0.45/max(nbins,1))

        corrMtrxCan = ROOT.TCanvas('c_corr_'+fittag,'c_corr_'+fittag,1000,900)
        corrMtrxCan.SetBottomMargin(0.20)
        corrMtrxCan.SetLeftMargin(0.15)
        corrMtrxCan.SetRightMargin(0.20)
        corrMtrxCan.SetTopMargin(0.08)

        corrMtrx.SetTitle('')
        corrMtrx.SetMinimum(-1)
        corrMtrx.SetMaximum(+1)
        corrMtrx.GetXaxis().SetLabelSize(labSize)
        corrMtrx.GetYaxis().SetLabelSize(labSize)
        corrMtrx.GetXaxis().LabelsOption('v')  # vertical x labels so long names don't overlap
        corrMtrx.GetZaxis().SetTitle('Correlation')
        corrMtrx.GetZaxis().SetTitleOffset(1.10)
        corrMtrx.GetZaxis().SetTitleSize(0.04)
        corrMtrx.GetZaxis().SetLabelSize(0.035)
        corrMtrx.SetMarkerColor(ROOT.kBlack)
        corrMtrx.SetMarkerSize(0.9)
        corrMtrx.Draw('colz text' if corrText else 'colz')

        # CMS Internal label (matches the GoF / TF / F-test plots)
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAlign(11)
        latex.SetTextFont(62)
        latex.SetTextSize(0.05)
        latex.DrawLatex(0.15,0.93,"CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.32,0.93,cmsText)
        if lumiText:
            latex.SetTextFont(42)
            latex.SetTextAlign(31)
            latex.SetTextSize(0.04)
            latex.DrawLatex(0.80,0.93,lumiText)

        corrMtrxCan.RedrawAxis()
        corrMtrxCan.Print('plots_fit_%s/correlation_matrix.png'%fittag,'png')
        corrMtrxCan.Print('plots_fit_%s/correlation_matrix.pdf'%fittag,'pdf')

        with open('plots_fit_%s/correlation_matrix.txt'%fittag,'w') as corrTxtFile:
            corrTxtFile.write(corrTxt)

    fit_result_file.Close()

def plot_gof(tag, subtag, seed=123456, condor=False):
    with cd(tag+'/'+subtag):
        if condor:
            tmpdir = 'notneeded/tmp/'
            execute_cmd('mkdir '+tmpdir) 
            execute_cmd('cat %s_%s_gof_toys_output_*.tgz | tar zxvf - -i --strip-components 2 -C %s'%(tag,subtag,tmpdir))
            toy_limit_tree = ROOT.TChain('limit')
            if len(glob.glob(tmpdir+'higgsCombine_gof_toys.GoodnessOfFit.mH120.*.root')) == 0:
                raise Exception('No files found')
            toy_limit_tree.Add(tmpdir+'higgsCombine_gof_toys.GoodnessOfFit.mH120.*.root') 
            
        else:
            toyOutput = ROOT.TFile.Open('higgsCombine_gof_toys.GoodnessOfFit.mH120.{seed}.root'.format(seed=seed))
            toy_limit_tree = toyOutput.Get('limit')

        # Now to analyze the output
        # Get observation
        ROOT.gROOT.SetBatch(True)
        ROOT.gStyle.SetOptStat(False)
        gof_data_file = ROOT.TFile.Open('higgsCombine_gof_data.GoodnessOfFit.mH120.root')
        gof_limit_tree = gof_data_file.Get('limit')
        gof_limit_tree.GetEntry(0)
        gof_data = gof_limit_tree.limit

        # Get toys. Bin ALL valid toys (TS>0) over a range that covers the toy distribution AND
        # the observed, so the drawn histogram integrates to the full toy count. The old
        # 'limit < 2*observed' cut hid most toys when the data fits very well (observed TS far
        # below the toy bulk -> e.g. a high p-value), making the plot look sparse / under-filled.
        print(gof_data)
        toy_max = toy_limit_tree.GetMaximum('limit')
        xhi = 1.15*max(toy_max, gof_data, 1.0)
        toy_limit_tree.Draw('limit>>hlimit(70,0,%f)'%xhi, 'limit>0')
        htoy_gof = ROOT.gDirectory.Get('hlimit')
        # p-value = fraction of *valid* toys with TS >= observed. The saturated GoF statistic is
        # chi2-like (skewed), so an empirical count beats a Gaussian fit. Toys whose fit failed
        # are recorded with TS <= 0; they are excluded from the denominator so a failed toy
        # doesn't bias the p-value low (the legend still reports the raw number thrown).
        n_thrown = int(toy_limit_tree.GetEntries())
        n_valid  = float(toy_limit_tree.GetEntries('limit>0'))
        n_above  = float(toy_limit_tree.GetEntries('limit>=%s'%gof_data))
        pvalue   = n_above/n_valid if n_valid > 0 else 0.0

        # Write out for reference
        with open('gof_results.txt','w') as out:
            out.write('Test statistic in data = %s\n'%gof_data)
            out.write('Toys thrown = %d\n'%n_thrown)
            out.write('Valid toys (TS>0) = %d\n'%int(n_valid))
            out.write('Toys >= observed = %d\n'%int(n_above))
            out.write('p-value = %s\n'%pvalue)


        # Canvas (CMS style)
        cout = ROOT.TCanvas('cout','cout',800,800)
        cout.SetLeftMargin(0.13)
        cout.SetRightMargin(0.05)
        cout.SetTopMargin(0.08)
        cout.SetBottomMargin(0.12)
        cout.SetTickx(1)
        cout.SetTicky(1)

        # Style the toy histogram
        htoy_gof.SetTitle('')
        htoy_gof.SetMarkerStyle(20)
        htoy_gof.SetMarkerSize(0.7)
        htoy_gof.SetMarkerColor(ROOT.kBlack)
        htoy_gof.SetLineColor(ROOT.kBlack)
        htoy_gof.GetXaxis().SetTitle('Saturated goodness-of-fit test statistic')
        htoy_gof.GetYaxis().SetTitle('Number of toys')
        htoy_gof.GetXaxis().SetTitleSize(0.045)
        htoy_gof.GetYaxis().SetTitleSize(0.045)
        htoy_gof.GetXaxis().SetTitleOffset(1.1)
        htoy_gof.GetYaxis().SetTitleOffset(1.35)
        htoy_gof.GetXaxis().SetLabelSize(0.04)
        htoy_gof.GetYaxis().SetLabelSize(0.04)
        htoy_gof.SetMaximum(1.25*htoy_gof.GetMaximum())

        # Gamma fit, overlaid for visualization. The saturated GoF statistic is chi2-like; a
        # gamma (chi2 with a free scale) hugs the peaked toy distribution much better than a
        # Gaussian or fixed-scale chi2. The reported p-value stays empirical.
        binwidth = htoy_gof.GetXaxis().GetBinWidth(1)
        mean = max(1.0, htoy_gof.GetMean())
        chi2fit = ROOT.TF1('gammafit','[0]*ROOT::Math::gamma_pdf(x,[1],[2],0)',
                           htoy_gof.GetXaxis().GetXmin(), htoy_gof.GetXaxis().GetXmax())
        chi2fit.SetParameters(htoy_gof.Integral()*binwidth, 2.0, mean/2.0)  # norm, shape, scale
        chi2fit.SetParName(0,'norm')
        chi2fit.SetParName(1,'shape')
        chi2fit.SetParName(2,'scale')
        # Wide shape range: as shape->inf the gamma -> Gaussian, so a symmetric/Gaussian-like
        # toy distribution is still fit well (gamma contains the Gaussian as a limiting case).
        chi2fit.SetParLimits(1, 0.1, 1000)
        chi2fit.SetParLimits(2, 0.001, 50)
        chi2fit.SetLineColor(ROOT.kAzure+1)
        chi2fit.SetLineWidth(2)
        htoy_gof.Fit(chi2fit,'QR')

        # Arrow for observed
        arrow = ROOT.TArrow(gof_data,0.25*htoy_gof.GetMaximum(),gof_data,0,0.02,'|>')
        arrow.SetLineWidth(3)
        arrow.SetLineColor(ROOT.kRed+1)
        arrow.SetFillColor(ROOT.kRed+1)

        # Legend
        leg = ROOT.TLegend(0.46,0.70,0.88,0.88)
        leg.SetLineWidth(0)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.035)
        # Total toys thrown (the drawn histogram may show fewer, since a display
        # cut of limit < 2*observed can crop far-tail toys).
        n_toys = int(toy_limit_tree.GetEntries())
        leg.AddEntry(htoy_gof,"Toy data (%d thrown)"%n_toys,"lep")
        leg.AddEntry(chi2fit,"Gamma fit","l")
        leg.AddEntry(arrow,"Observed = %.1f"%gof_data,"l")
        leg.AddEntry(0,"#it{p}-value = %.2f"%(pvalue),"")

        # Draw
        htoy_gof.Draw('pe')
        chi2fit.Draw('same')
        arrow.Draw()
        leg.Draw()

        # CMS Internal label
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAlign(11)
        latex.SetTextFont(62)
        latex.SetTextSize(0.05)
        latex.DrawLatex(0.13,0.93,"CMS")
        latex.SetTextFont(52)
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.23,0.93,"Internal")

        cout.RedrawAxis()
        cout.Print('gof_plot.pdf','pdf')
        cout.Print('gof_plot.png','png')

    if condor:
            execute_cmd('rm -r '+tmpdir)

def plot_transfer_funcs(tag, subtag, projection='projx2', savename=None, cmsText='Internal', xtitle='p_{T} [GeV]'):
    '''Plot the transfer factor (postfit pass/fail ratio) vs pT for one fit area.

    Reads the b-only postfit projections from <tag>/<subtag>/plots_fit_b/all_plots.root
    and divides the pass projection by the fail projection to form the transfer factor.
    Integrated (and CMS-restyled) from the former standalone helper_scripts/plot_TF.py.

    Args:
        tag (str): working area (e.g. 'rpfmult_Binningv9_Inputv25_SR_BkgMC').
        subtag (str): fit area (e.g. 'Signal_M3000GeV_e4_SR-2x0_area').
        projection (str, optional): which postfit x-projection to use. Defaults to 'projx2'.
        savename (str, optional): output path stem (no extension). Defaults to
            <tag>/<subtag>/transfer_func.
        cmsText (str, optional): sublabel drawn after 'CMS'. Defaults to 'Internal'.
        xtitle (str, optional): x-axis title. Defaults to 'p_{T} [GeV]'.
    '''
    rundir = tag+'/'+subtag
    fpath = rundir+'/plots_fit_b/all_plots.root'
    if not os.path.exists(fpath):
        raise Exception('Postfit shapes not found: %s (run the fit + StdPlots first)'%fpath)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(False)

    f = ROOT.TFile.Open(fpath)
    h_pass = f.Get('TotalBkg_pass_postfit_'+projection)
    h_fail = f.Get('TotalBkg_fail_postfit_'+projection)
    if not h_pass or not h_fail:
        f.Close()
        raise Exception('Missing TotalBkg_{pass,fail}_postfit_%s in %s'%(projection,fpath))
    h_ratio = h_pass.Clone('h_TF_'+subtag)
    h_ratio.Divide(h_fail)
    h_ratio.SetDirectory(0)  # detach so it survives the file close
    f.Close()

    # Canvas (CMS style, matches the GoF / F-test plots)
    c = ROOT.TCanvas('cTF','cTF',800,800)
    c.SetLeftMargin(0.19)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)

    h_ratio.SetTitle('')
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetLineWidth(2)
    h_ratio.SetMarkerStyle(0)
    h_ratio.GetXaxis().SetTitle(xtitle)
    h_ratio.GetYaxis().SetTitle('Transfer Factor')
    h_ratio.GetXaxis().SetTitleSize(0.045)
    h_ratio.GetYaxis().SetTitleSize(0.045)
    h_ratio.GetXaxis().SetTitleOffset(1.1)
    # 1.35 overlapped the tick labels whenever the transfer factor is small enough that
    # ROOT writes them out in full (e.g. "0.005" in the closure fit) rather than pulling
    # out an exponent.
    h_ratio.GetYaxis().SetTitleOffset(1.9)
    h_ratio.GetXaxis().SetLabelSize(0.04)
    h_ratio.GetYaxis().SetLabelSize(0.04)
    # Give the curve headroom instead of letting ROOT put the frame top exactly on the
    # maximum bin. With SetMinimum(0) and no maximum, a CONSTANT transfer factor (par1
    # railed at its lower bound, as in the SR_BkgMC closure fit) draws as a single line
    # lying on the top frame edge over an otherwise empty panel, which reads as a failed
    # plot rather than as a flat TF. A curved TF hid the problem because its maximum sits
    # at one pT only.
    _hi = h_ratio.GetMaximum()
    h_ratio.SetMinimum(0)
    h_ratio.SetMaximum(1.3*_hi if _hi > 0 else 1.0)
    h_ratio.Draw('hist ][ C')

    # CMS Internal label
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.19,0.93,"CMS")
    latex.SetTextFont(52)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.30,0.93,cmsText)

    c.RedrawAxis()
    if savename is None:
        savename = rundir+'/transfer_func'
    c.Print(savename+'.pdf','pdf')
    c.Print(savename+'.png','png')

def plot_ftest(base_fstat, nRpfs1, nRpfs2, nBins, poly1, poly2, savedir, cmsText='Internal'):
    '''Plot the Fisher F-test comparing two transfer-function orders (no toys).

    Draws the analytic F-distribution with ndf=(p2-p1, nBins-p2), marks the observed
    F-statistic, and reports the p-value as its upper-tail integral. Saves
    ftest_<poly1>_vs_<poly2>_notoys.{png,pdf} under savedir.

    Args:
        base_fstat (list): output of TwoDAlphabet.ftest.FstatCalc; base_fstat[0] is the observed F.
        nRpfs1, nRpfs2 (int): number of RPF parameters of the two fits.
        nBins (int): number of bins entering the fit.
        poly1, poly2 (str): TF order labels (e.g. '1x0', '2x0').
        savedir (str): directory to write the output plots into (the working area).
        cmsText (str, optional): sublabel after 'CMS'. Defaults to 'Internal'.
    '''
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    if len(base_fstat) == 0:
        base_fstat = [0.0]

    ftest_p1    = min(nRpfs1,nRpfs2)
    ftest_p2    = max(nRpfs1,nRpfs2)
    xmax        = max(10,1.3*base_fstat[0])
    fdist       = ROOT.TF1("fDist", "[0]*TMath::FDist(x, [1], [2])", 0, xmax)
    fdist.SetParameter(0,1)
    fdist.SetParameter(1,ftest_p2-ftest_p1)
    fdist.SetParameter(2,nBins-ftest_p2)

    pval = fdist.Integral(0.0,base_fstat[0])
    print('F-test p-value: %s'%(1-pval))

    # Canvas (CMS style, matches the GoF plots)
    c = ROOT.TCanvas('cFtest','cFtest',800,800)
    c.SetLeftMargin(0.13)
    c.SetRightMargin(0.05)
    c.SetTopMargin(0.08)
    c.SetBottomMargin(0.12)
    c.SetTickx(1)
    c.SetTicky(1)

    # The F-dist PDF diverges at x->0 when ndf1==1, so fdist.GetMaximum() is unusable.
    # Cap the y-axis where the bulk of the distribution lives.
    ymax = 1.1

    frame = ROOT.TH1F("Fhist","",100,0,xmax)
    frame.GetXaxis().SetTitle("F statistic")
    frame.GetYaxis().SetTitle("Probability density")
    frame.GetXaxis().SetTitleSize(0.045)
    frame.GetYaxis().SetTitleSize(0.045)
    frame.GetXaxis().SetTitleOffset(1.1)
    frame.GetYaxis().SetTitleOffset(1.35)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.SetMaximum(ymax)
    frame.SetMinimum(0)
    frame.Draw()

    fdist.SetLineColor(ROOT.kAzure+1)
    fdist.SetLineWidth(2)
    fdist.Draw('same')

    # Observed arrow (red, matches GoF observed marker)
    ftestobs = ROOT.TArrow(base_fstat[0],0.25*ymax,base_fstat[0],0,0.02,'|>')
    ftestobs.SetLineColor(ROOT.kRed+1)
    ftestobs.SetFillColor(ROOT.kRed+1)
    ftestobs.SetLineWidth(3)
    ftestobs.Draw()

    tLeg = ROOT.TLegend(0.46,0.70,0.88,0.88)
    tLeg.SetLineWidth(0)
    tLeg.SetBorderSize(0)
    tLeg.SetFillStyle(0)
    tLeg.SetTextFont(42)
    tLeg.SetTextSize(0.035)
    tLeg.AddEntry(fdist,"F-dist, ndf = (%.0f, %.0f)"%(fdist.GetParameter(1),fdist.GetParameter(2)),"l")
    tLeg.AddEntry(ftestobs,"Observed = %.1f"%base_fstat[0],"l")
    tLeg.AddEntry(0,"%s vs. %s"%(poly1,poly2),"")
    tLeg.AddEntry(0,"#it{p}-value = %.2f"%(1-pval),"")
    tLeg.Draw("same")

    # CMS Internal label
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.SetTextSize(0.05)
    latex.DrawLatex(0.13,0.93,"CMS")
    latex.SetTextFont(52)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.23,0.93,cmsText)

    c.RedrawAxis()
    c.SaveAs(savedir+'/ftest_{0}_vs_{1}_notoys.png'.format(poly1,poly2))
    c.SaveAs(savedir+'/ftest_{0}_vs_{1}_notoys.pdf'.format(poly1,poly2))

def plot_signalInjection(tag, subtag, injectedAmount, seed=123456, stats=True, condor=False):
    # if injectedAmount is not an integer, need to look for different file
    # see: https://github.com/lcorcodilos/2DAlphabet/blob/e089ed1da63172770726b3e6f406c11e611e057d/TwoDAlphabet/twoDalphabet.py#L488
    injectedName = str(injectedAmount).replace('.','p')
    with cd(tag+'/'+subtag):
        if condor:
            tmpdir = 'notneeded/tmp/'
            execute_cmd('mkdir '+tmpdir) 
            execute_cmd('cat %s_%s_sigInj_r%s_output_*.tgz | tar zxvf - -i --strip-components 2 -C %s'%(tag,subtag,injectedName,tmpdir))
            tree_fit_sb = ROOT.TChain('tree_fit_sb')
            tree_fit_sb.Add(tmpdir+'fitDiagnostics_sigInj_r{rinj}*.root'.format(rinj=injectedName)) 
            
        else:
            toyOutput = ROOT.TFile.Open('fitDiagnostics_sigInj_r{rinj}_{seed}.root'.format(rinj=injectedName,seed=seed))
            tree_fit_sb = toyOutput.Get('tree_fit_sb')

        ROOT.gROOT.SetBatch(True)
        if stats:
            ROOT.gStyle.SetOptStat(True)
        else:
            ROOT.gStyle.SetOptStat(False)
        # Final plotting
        result_can = ROOT.TCanvas('sigpull_can','sigpull_can',800,700)

        tree_fit_sb.Draw("(r-{rinj})/(rHiErr*(r<{rinj})+rLoErr*(r>{rinj}))>>sigpull(20,-5,5)".format(rinj=injectedAmount),"fit_status>=0")
        hsigpull = ROOT.gDirectory.Get('sigpull')
        tree_fit_sb.Draw("(r-{rinj})>>sigstrength(20,-1,1)".format(rinj=injectedAmount),"fit_status>=0")
        hsignstrength = ROOT.gDirectory.Get('sigstrength')

        hsigpull.Fit("gaus","L")
        hsigpull.SetTitle('')
        hsigpull.GetXaxis().SetTitle('(r-%s)/rErr'%injectedAmount)
        result_can.cd()
        hsigpull.Draw('pe')
        result_can.Print('signalInjection_r%s_pull.png'%(str(injectedAmount).replace('.','p')),'png')

        hsignstrength.Fit("gaus","L")
        hsignstrength.SetTitle('')
        hsignstrength.GetXaxis().SetTitle('r-%s'%injectedAmount)
        result_can.cd()
        hsignstrength.Draw('pe')
        result_can.Print('signalInjection_r%s.png'%(str(injectedAmount).replace('.','p')),'png')

    if condor:
            execute_cmd('rm -r '+tmpdir)


'''
def plotRPF(postfitShapesFile,odir,qcdTag,passTag="M",failTag="F",xRange=[60,200],yTitle="$R_{M/F}$"):
    hPass2D = get2DPostfitPlot(postfitShapesFile,qcdTag,passTag)
    hFail2D = get2DPostfitPlot(postfitShapesFile,"qcd",failTag)
    nBins   = hPass2D.GetNbinsY()
    foutName = "R{0}{1}".format(passTag,failTag)
    maxRpf  = 0.
    histos  = []
    if(nBins==2):
        legendFlag  = True
        labels      = ["300<$p_T$<400 GeV","$p_T$>400 GeV"]
    else:
        labels      = ["_dummy"]
        legendFlag  = False

    for i in range(1,nBins+1):
        hPass = hPass2D.ProjectionX("hPass_temp",i,i)
        hFail = hFail2D.ProjectionX("hFail_temp",i,i)
        hRpf  = hPass.Clone("hRPF")
        hRpf.Divide(hFail)
        hRpf.Scale(1000)
        maxRpf = max(maxRpf,hRpf.GetMaximum())
        hRpf, edges = hist2array(hRpf,return_edges=True)
        histos.append(hRpf)

    plt.style.use([hep.style.CMS])
    f, ax = plt.subplots()
    hep.histplot(histos, bins=edges[0],label=labels)
    if legendFlag:
        plt.legend()

    hep.cms.text("Internal",loc=0)
    
    ax.set_xlim(xRange)
    ax.set_ylim([0.,maxRpf*1.3])

    plt.xlabel("$M_{PNet}$ (GeV)",horizontalalignment='right', x=1.0)
    plt.ylabel(yTitle+" x $10^{3}$",horizontalalignment='right', y=1.0)
    #ax.yaxis.set_tick_params(which='minor', left=False)    
    #ax.yaxis.set_tick_params(which='minor', right=False)    

    print("Saving "+odir+"{0}.pdf".format(foutName))
    plt.savefig(odir+"/{0}.pdf".format(foutName), bbox_inches='tight')
    plt.savefig(odir+"/{0}.png".format(foutName), bbox_inches='tight')
    plt.cla()
    plt.clf()
'''

def poisson_conf_interval(k):
    """
    Calculate Poisson (Garwood) confidence intervals using ROOT's TH1 with kPoisson error option.
    
    Parameters:
    k (array): The number of counts (events) per bin.

    Returns:
    lower (array): Bin count - lower error.
    upper (array): Bin count + upper error.
    """
    lower = np.zeros_like(k, dtype=float)
    upper = np.zeros_like(k, dtype=float)
    #Temp hist to exploit ROOT's built-in CI calculating function
    hist = ROOT.TH1F("hist_delete", "", 1, 0, 1)
    hist.SetBinErrorOption(ROOT.TH1.kPoisson)
    hist.Sumw2()

    for i, count in enumerate(k):
        hist.SetBinContent(1, count)
        
        lower[i] = hist.GetBinContent(1) - hist.GetBinErrorLow(1)
        upper[i] = hist.GetBinContent(1) + hist.GetBinErrorUp(1)
        
    hist.Delete()
    
    return lower, upper
