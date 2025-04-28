#!/bin/env python

import math
import os
import json
import shutil
import argparse
import logging
import colorlog

LOG_LEVEL = logging.DEBUG
LOGFORMAT = "  %(log_color)s%(levelname)-8s%(reset)s | %(log_color)s%(message)s%(reset)s"
logging.root.setLevel(LOG_LEVEL)
formatter = colorlog.ColoredFormatter(LOGFORMAT)
stream = logging.StreamHandler()
stream.setLevel(LOG_LEVEL)
stream.setFormatter(formatter)
logger = logging.getLogger("Calculators42HDM")
logger.setLevel(LOG_LEVEL)
logger.addHandler(stream)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib as mpl
mpl.use('Agg')

from packaging import version
if version.parse(mpl.__version__) >= version.parse('2.0.0'):
    # Override some of matplotlib 2 new style
    mpl.rcParams['grid.color'] = 'k'
    mpl.rcParams['grid.linestyle'] = 'dotted'
    mpl.rcParams['grid.linewidth'] = 0.5
    mpl.rcParams['figure.figsize'] = [8.0, 6.0]
    mpl.rcParams['figure.dpi'] = 80
    mpl.rcParams['savefig.dpi'] = 100
    #mpl.rcParams['font.size'] = 12
    mpl.rcParams['legend.fontsize'] = 35 #'large'
    mpl.rcParams['figure.titlesize'] = 'medium'
    mpl.rcParams['lines.linewidth'] = 1.0
    mpl.rcParams['lines.dashed_pattern'] = [6, 6]
    mpl.rcParams['lines.dashdot_pattern'] = [3, 5, 1, 5]
    mpl.rcParams['lines.dotted_pattern'] = [1, 3]
    mpl.rcParams['lines.scale_dashes'] = False

from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import interpolate
from scipy.interpolate import make_interp_spline, BSpline
from matplotlib.gridspec import GridSpec

#
import CMSStyle

from cp3_llbb.Calculators42HDM.Calc2HDM import *
from cp3_llbb.Calculators42HDM.labellines import *



def float_to_str(x, digits=2):
    tmp = ':.{:d}f'.format(digits)
    tmp = ('{' + tmp + '}').format(x)
    return tmp.replace('.', 'p')

def mass_to_string(m):
    r = '{:.2f}'.format(m)
    r = r.replace('.', 'p')
    return r

def ZABR(type2hdm= 2, func_of_cba=False, func_of_tb=False, func_of_mA=False, func_of_mH=False):
    sqrts = 13000
    type = type2hdm
    mh = 125.
    mZ = 91.1876
    if func_of_tb:
        cba = 0.13
        sba = math.sqrt(1 - pow(cba, 2))
        dicvars = {'tb': [] } 
        mA = 465
        mH = 660
    elif func_of_cba:
        tb = 1.5
        dicvars = {'cba': [] } 
        mA = 465
        mH = 660
    elif func_of_mA:
        cba = 0.13
        sba = math.sqrt(1 - pow(cba, 2))
        tb= 5.
        dicvars = {'mA': [] } 
    elif func_of_mH:
        cba = 0.13
        sba = math.sqrt(1 - pow(cba, 2))
        tb= 1.5
        dicvars = {'mH': [] } 
  
    outputpath = '/home/ucl/cp3/kjaffel/scratch/Calculators42HDM_outputs'
    cwd = os.getcwd()
    if not os.path.exists(outputpath):
        os.makedirs(outputpath)
    results = {
        'BRHtoZA': [],
        'BRHtoZga': [],
        'BRHtobb': [],
        'BRHtoss': [],
        'BRHtocc': [],
        'BRHtoee': [],
        'BRHtomumu': [],
        'BRHtohh': [],
        'BRHtoWW': [],
        'BRHtoZZ': [],
        'BRHtott': [],
        'BRHtotautau': [],
        'BRHtogg': [],
        'BRHtogluglu': [],
        
        'BRAtobb': [],
        'BRAtoee': [],
        'BRAtomumu': [],
        'BRAtoss': [],
        'BRAtocc': [],
        'BRAtotautau': [],
        'BRAtogluglu': [],
        'BRAtoZh': [],
        'BRAtoZga': [],
        'BRAtogg' : []
            }
    results.update( dicvars )
    print (results)
    xaxis = ( np.arange(0.13, math.pi, 0.01) if func_of_cba 
              else (np.arange(0.01, 20.1, 0.001) if func_of_tb 
                  else(np.arange(30.,1100.,10.) if func_of_mA 
                      else (np.arange(120.,1190.,10.)))))
    for x in xaxis : 
        if func_of_tb:
            tb = x
            results['tb'].append(tb)
        elif func_of_cba:
            cba = math.cos(x)
            if x >= 0 and x <= math.pi/2:
                sba = math.sin(x)
            elif x > math.pi/2 and x <= math.pi:
                sba = math.sin(x-math.pi)
            results['cba'].append(cba)
        elif func_of_mA:
            mA = x
            mH = mA + mZ
            results['mA'].append(mA)
        elif func_of_mH:
            mH= x
            mA= mH - mZ
            results['mH'].append(mH)

        if mA > mH:
            logger.info("MA_{} > MH_{} switching to A->ZH mode!".format(mA, mH))
            mode = 'A'
        elif mH >= mA and mH> 125.:
            logger.info("MA_{} =< MH_{} switching to H->ZA mode!".format(mA, mH))
            mode = 'H'
        elif mH >= mA and mH <= 125.:
            logger.info("MA_{} >= MH_{} && H <= 125. GeV switching to h->ZH mode!".format(mA, mH))
            mode ='h'
        
        outputFile = 'out.dat'
        if options.debug:    
            outputFile = 'out_mH-{}_mA-{}_tb-{}_cba-{}.dat'.format(
                    mass_to_string(mH), mass_to_string(mA), mass_to_string(tb), mass_to_string(cba))
        
        mhc = max(mH, mA)
        m12 = math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
        res = Calc2HDM(mode  = mode,
                sqrts = sqrts,
                type  = type,
                tb    = tb,
                m12   = m12,
                mh    = mh,
                mH    = mH,
                mA    = mA,
                mhc   = mhc,
                sba   = sba,
                outputFile = outputFile,
                muR4ggh = muR4ggh,
                muF4ggh = muF4ggh,
                muR4bbh = muR4bbh,
                muF4bbh = muF4bbh
                )
 
        res.setpdf('NNPDF31_nnlo_as_0118_mc_hessian_pdfas')
        res.computeBR()

        results['BRHtoss'].append(res.HtossBR)
        results['BRHtocc'].append(res.HtoccBR)
        results['BRHtobb'].append(res.HtobbBR)
        results['BRHtott'].append(res.HtottBR)
        results['BRHtoee'].append(res.HtoeeBR)
        results['BRHtomumu'].append(res.HtomumuBR)
        results['BRHtotautau'].append(res.HtotautauBR)
        results['BRHtogg'].append(res.HtoggBR) # gamma-gamma
        results['BRHtoZZ'].append(res.HtoZZBR)
        results['BRHtoWW'].append(res.HtoWWBR)
        results['BRHtoZga'].append(res.HtoZgaBR)
        results['BRHtogluglu'].append(res.HtoglugluBR)
        #results['BRHtohh'].append(res.HtohhBR)
        results['BRHtoZA'].append(res.HtoZABR)
        results['BRAtoss'].append(res.AtossBR)
        results['BRAtocc'].append(res.AtoccBR)
        results['BRAtobb'].append(res.AtobbBR)
        results['BRAtoee'].append(res.AtoeeBR)
        results['BRAtomumu'].append(res.AtomumuBR)
        results['BRAtotautau'].append(res.AtotautauBR)
        results['BRAtogg'].append(res.AtoggBR)
        results['BRAtoZga'].append(res.AtoZgaBR)
        results['BRAtogluglu'].append(res.AtoglugluBR)
        results['BRAtoZh'].append(res.AtoZhBR)
        
        if options.debug:
            shutil.move(os.path.join(cmssw, outputFile), os.path.join( outputpath, outputFile ))
            shutil.move(os.path.join(cmssw, outputFile.replace('.dat', '.log')), outputpath)
    
    var = (cba if func_of_tb else(tb))
    if func_of_tb or func_of_cba:
        jsonf = 'BR_mH-%s_mA-%s_2hdmtype-%s_%s-%s_function_of_%s.json' % (mass_to_string(mH), mass_to_string(mA), 
                                                                            type2hdm,
                                                                        ('cba' if func_of_tb else('tb')),
                                                                        float_to_str(var, 2), 
                                                                        ('tb' if func_of_tb else('cba')))
    else:
        jsonf = 'BR_2hdmtype-%s_tb-%s_cba-%s_function_of_%s.json'% ( #('mA' if func_of_mH else ('mH')),
                                                                     #      mass_to_string(mA) if func_of_mH else( mass_to_string(mH)),
                                                                           type2hdm,
                                                                           float_to_str(tb, 2),
                                                                           float_to_str(cba, 2),
                                                                           ('mA' if func_of_mA else('mH')))
    with open(jsonf, 'w+') as f:
        json.dump(results, f)
    return jsonf


def MakeBRPlots( plots=[] ): 
    plt.rc('axes', labelsize=23, titlesize=23)
    fig, axs = plt.subplots(2, 2, figsize=(23, 15), dpi=300, gridspec_kw={'height_ratios': [1.3, 1]})
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.27, hspace=None)
    
    for plot in plots:
        x = []
        y_HtoZA = []
        y_Htobb = []
        y_Htohh = []
        y_HtoWW = []
        y_HtoZZ = []
        y_Htott = []
        y_Htotautau = []
        y_Htogg = []
        y_Htogluglu = []
        y_Atobb = []
        y_Atotautau = []
        y_Atogluglu = []
        y_Atogg = []
        y_AtoZh = []

        if options.lazy:
            input_files = {'func_of_cba': 'data/BR_mH-300p00_mA-200p00_2hdmtype-2_tb-1p50_function_of_cba.json', 
                            'func_of_tb': 'data/BR_mH-300p00_mA-200p00_2hdmtype-2_cba-0p01_function_of_tb.json',
                            'func_of_mA': 'BR_2hdmtype-2_tb-20p00_cba-0p01_function_of_mA.json',
                            'func_of_mH': 'BR_2hdmtype-2_tb-1p50_cba-0p01_function_of_mH.json'
                        }
        else:
            input_files = {'func_of_cba': ZABR(type2hdm=options.type, func_of_cba=True, func_of_tb=False, func_of_mA=False, func_of_mH=False), 
                           'func_of_tb': ZABR(type2hdm=options.type, func_of_cba=False, func_of_tb=True, func_of_mA=False, func_of_mH=False),
                            #'func_of_mA': ZABR(type2hdm=options.type, func_of_cba=False, func_of_tb=False, func_of_mA=True, func_of_mH=False),
                            #'func_of_mH': ZABR(type2hdm=options.type, func_of_cba=False, func_of_tb=False, func_of_mA=False, func_of_mH=True)
                        }

        with open(input_files[plot]) as f:
            theory = json.load(f)

        x = (theory['tb'] if plot =='func_of_tb' else ( 
             theory['cba'] if plot =='func_of_cba' else ( 
             theory['mA'] if plot=='func_of_mA' else ( 
             theory['mH']))))
            
        y_HtoZA = theory['BRHtoZA']
        y_HtoZga= theory['BRHtoZga']
        y_Htobb = theory['BRHtobb']
        y_Htocc = theory['BRHtocc']
        y_Htoss = theory['BRHtoss']
        y_Htoee = theory['BRHtoee']
        y_Htohh = theory['BRHtohh']
        y_HtoWW = theory['BRHtoWW']
        y_HtoZZ = theory['BRHtoZZ']
        y_Htott = theory['BRHtott']
        y_Htogg = theory['BRHtogg']
        y_Htogluglu = theory['BRHtogluglu']
        y_Htotautau = theory['BRHtotautau']
        y_Htomumu   = theory['BRHtomumu']
        y_Atobb = theory['BRAtobb']
        y_Atocc = theory['BRAtocc']
        y_Atoss = theory['BRAtoss']
        y_Atoee = theory['BRAtoee']
        y_AtoZga = theory['BRAtoZga']
        y_Atomumu = theory['BRAtomumu']
        y_Atotautau = theory['BRAtotautau']
        y_Atogluglu = theory['BRAtogluglu']
        y_Atogg = theory['BRAtogg']
        y_AtoZh = theory['BRAtoZh']
        x = np.asarray(x)
        #xnew = np.linspace(x.min(),x.max(),len(x))
        y_HtoZA = np.asarray(y_HtoZA)
        y_HtoZga = np.asarray(y_HtoZga)
        y_Htobb = np.asarray(y_Htobb) 
        y_Htocc = np.asarray(y_Htocc) 
        y_Htoss = np.asarray(y_Htoss) 
        y_Htoee = np.asarray(y_Htoee) 
        y_Htomumu = np.asarray(y_Htomumu) 
        y_Htohh = np.asarray(y_Htohh)
        y_HtoWW = np.asarray(y_HtoWW)
        y_HtoZZ = np.asarray(y_HtoZZ)
        y_Htott = np.asarray(y_Htott)
        #spl = make_interp_spline(xnew, y_Htott, k=1)
        #y_Htott_smooth = spl(xnew)
        y_Htotautau = np.asarray(y_Htotautau) 
        y_Htogg = np.asarray(y_Htogg)
        y_Htogluglu = np.asarray(y_Htogluglu)
        y_Atobb = np.asarray(y_Atobb)
        y_Atocc = np.asarray(y_Atocc)
        y_Atoss = np.asarray(y_Atoss)
        y_Atoee = np.asarray(y_Atoee)
        y_Atomumu = np.asarray(y_Atomumu)
        y_AtoZga = np.asarray(y_AtoZga)
        y_Atotautau = np.asarray(y_Atotautau)
        y_Atogluglu = np.asarray(y_Atogluglu)
        y_Atogg = np.asarray(y_Atogg)
        y_AtoZh = np.asarray(y_AtoZh)

        xmin = (0.01 if plot == 'func_of_tb' else (-1. if plot == 'func_of_cba' else (30. if  plot == 'func_of_mA' else(120.))))
        xmax = (20.1 if plot == 'func_of_tb' else (1. if plot == 'func_of_cba' else (1100. if  plot == 'func_of_mA' else(1190.))))
        print ( plot, xmin, xmax )
        
        # BR
        ymin = 0.000000001
        ymax = 10.
    
        #CMSStyle.changeFont()
        npoints = len(x)
        if plot == "func_of_cba" or plot == 'func_of_mA':
            rowH = 0
            colH = 0
            rowA = 1
            colA = 0
        else:
            rowH = 0
            colH = 1
            rowA = 1
            colA = 1
        
        #First subplot
        print (rowH,colH,rowA,colA)
        axs[rowH,colH].plot(x, y_HtoZA, 'red', label='ZA', linewidth=2.5)
        axs[rowH,colH].plot(x, y_HtoZga, 'violet', label=r'Z$\gamma$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htobb, 'purple', label=r'b$\bar{b}$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htocc, 'blueviolet', label=r'c$\bar{c}$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htoss, 'deeppink', label=r's$\bar{s}$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htoee, 'darkred', label=r'$e^+e^-$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htomumu, 'navy', label=r'$\mu^+\mu^-$', linewidth=2.5)
        #axs[rowH,colH].plot(x, y_Htohh, 'darkgoldenrod',label='hh', linewidth=2.5)
        axs[rowH,colH].plot(x, y_HtoWW, 'green', label='WW', linewidth=2.5)
        axs[rowH,colH].plot(x, y_HtoZZ, 'lightskyblue', label='ZZ', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htott, 'mediumblue', label=r't$\bar{t}$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htotautau, 'palevioletred', label=r'$\tau \tau$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htogg, 'lime', label=r'$\gamma \gamma$', linewidth=2.5)
        axs[rowH,colH].plot(x, y_Htogluglu, 'aqua', label='gg', linewidth=2.5)
        axs[rowH,colH].xaxis.set_tick_params(which='both', labelbottom=False, labeltop=False)
        axs[rowH,colH].xaxis.offsetText.set_visible(False)

        print (xmin, xmax, ymin, ymax )
        axs[rowH,colH].axis([xmin, xmax, ymin, ymax])
        if options.logy:
            axs[rowH,colH].set_yscale('log')
        if plot in ["func_of_tb", "func_of_mH", "func_of_mA"] and options.logy:
            axs[rowH,colH].set_xscale('log')
        
        #if plot == "func_of_tb":
        #    axs[rowH,colH].legend(loc='center right', bbox_to_anchor=(1.33, 0.5), prop={'size': 30}, frameon=False)
        if plot == "func_of_tb":
            axs[rowH,colH].text(0.02,  1.02, 
                    r"2HDM-Type{}: $m_H=300$ GeV, $m_A=200$ GeV, $cos(\beta-\alpha) = 0.01$".format('II' if options.type==2 else('I')), 
                    fontsize=18, color='black', transform=axs[rowH,colH].transAxes)
        elif plot == "func_of_cba":
            axs[rowH,colH].text(0.02,  1.02, 
                    r"2HDM-Type{}: $m_H=300$ GeV, $m_A=200$ GeV, $tan\beta = 1.5$".format('II' if options.type==2 else('I')), 
                    fontsize=18, color='black', transform=axs[rowH,colH].transAxes)
        elif plot == "func_of_mA":
            axs[rowH,colH].text(0.02,  1.02, 
                    r"2HDM-Type{}: $m_H= m_Z + m_A$, $cos(\beta-\alpha) = 0.01$, $tan\beta = 20.$".format('II' if options.type==2 else('I')), 
                    fontsize=18, color='black', transform=axs[rowH,colH].transAxes)
        elif plot == "func_of_mH":
            axs[rowH,colH].text(0.02,  1.02, 
                    r"2HDM-Type{}: $m_A= m_H - m_Z$, $cos(\beta-\alpha) = 0.01$, $tan\beta = 1.5$".format('II' if options.type==2 else('I')), 
                    fontsize=18, color='black', transform=axs[rowH,colH].transAxes)
        
        #Second subplot
        axs[rowA,colA].plot(x, y_AtoZh, 'darkslategray', label='Zh',               linewidth=2.5)
        axs[rowA,colA].plot(x, y_AtoZga, 'violet', label=r'Z$\gamma$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atobb, 'purple', label=r'b$\bar{b}$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atocc, 'blueviolet', label=r'c$\bar{c}$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atoss, 'deeppink', label=r's$\bar{s}$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atoee, 'darkred', label=r'$e^+e^-$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atomumu, 'navy', label=r'$\mu^+\mu^-$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atotautau, 'palevioletred', label=r'$\tau \tau$', linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atogg, 'lime', label=r'$\gamma \gamma$',          linewidth=2.5)
        axs[rowA,colA].plot(x, y_Atogluglu, 'aqua', label='gg',                    linewidth=2.5)
        axs[rowA,colA].axis([xmin, xmax, 0.000000001, ymax])
    
        if options.logy:
            axs[rowA,colA].set_yscale('log')
        if plot in ["func_of_tb"] and options.logy:
            axs[rowA,colA].set_xscale('log')
        
        if plot == "func_of_cba":
            majorLocator   = MultipleLocator(0.5)
            majorFormatter = FormatStrFormatter('%.1f')
            minorLocator   = MultipleLocator(0.25)
    
            axs[rowA,colA].xaxis.set_major_locator(majorLocator)
            axs[rowA,colA].xaxis.set_major_formatter(majorFormatter)
            axs[rowA,colA].xaxis.set_minor_locator(minorLocator)
            #ticks = np.arange(-1, 1.01, 0.5)
            #axs[rowA,colA].set_xticks(ticks, minor=False)
        
        #if plot == "func_of_tb":
        #    axs[rowA,colA].legend(loc='center right', bbox_to_anchor=(1.29, 0.5), prop={'size': 30}, frameon=False)
        plt.setp(axs[rowH,colH].get_xticklabels(), fontsize=24)
        plt.setp(axs[rowH,colH].get_yticklabels(), fontsize=24)
        labelLines(axs[rowH,colH].get_lines(),align=False,fontsize=18)
        plt.setp(axs[rowA,colA].get_xticklabels(), fontsize=24)
        plt.setp(axs[rowA,colA].get_yticklabels(), fontsize=24)
        labelLines(axs[rowA,colA].get_lines(),align=False,fontsize=18)
    
        if plot == "func_of_cba":
            axs[rowA,colA].set_xlabel(r'cos($\beta-\alpha$)', fontsize=26, horizontalalignment='right', x=0.56)
            axs[rowH,colH].set_xlabel(r'cos($\beta-\alpha$)', fontsize=26, horizontalalignment='right', x=0.56)
        elif plot == "func_of_tb":
            axs[rowA,colA].set_xlabel(r'tan$\beta$', fontsize=26, horizontalalignment='right', x=0.56)
            axs[rowH,colH].set_xlabel(r'tan$\beta$', fontsize=26, horizontalalignment='right', x=0.56)
        elif plot == "func_of_mA":
            axs[rowA,colA].set_xlabel(r'$m_A$ [GeV]', fontsize=26, horizontalalignment='right', x=0.56)
            axs[rowH,colH].set_xlabel(r'$m_A$ [GeV]', fontsize=26, horizontalalignment='right', x=0.56)
        elif plot == "func_of_mH":
            axs[rowA,colA].set_xlabel(r'$m_H$ [GeV]', fontsize=26, horizontalalignment='right', x=0.56)
            axs[rowH,colH].set_xlabel(r'$m_H$ [GeV]', fontsize=26, horizontalalignment='right', x=0.56)
    
        axs[rowH,colH].set_ylabel(r'BR($H \to XX$)', fontsize=26, horizontalalignment='right', y=1.0)
        axs[rowH,colH].grid()
        axs[rowA,colA].set_ylabel(r'BR($A \to XX$)', fontsize=26, horizontalalignment='right', y=1.0)
        axs[rowA,colA].grid()
        
        #fig.tight_layout()
        #suffix= ('cosbeta-alpha_and_tb' if 'tb' in plot else ("MA_and_MH"))
    
    suffix ='new' 
    if fig:
        fig.savefig('2HDM{}_BRs_func_{}.pdf'.format('II' if options.type==2 else('I'), suffix), bbox_inches='tight', pad_inches=0)
        fig.savefig('2HDM{}_BRs_func_{}.png'.format('II' if options.type==2 else('I'), suffix), bbox_inches='tight', pad_inches=0)
        print ("plot saved in:",'forhh_2HDM{}_BRs_func_{}.png'.format('II' if options.type==2 else('I'), suffix)) 
        # clean the figure before next plot
        plt.gcf().clear()
    return fig


if __name__ == '__main__': 
    
    parser = argparse.ArgumentParser(description='Create Validation Plots for 2HDM BR')

    parser.add_argument('-lazy', '--lazy', action='store_true', required=False,
            help='This will allow you to do fast run relying on old jsons files with calling 2HDMC-1.8.0/CalcPhys/ Class')
    parser.add_argument('-s', '--scan', action='store_true', required=False,
            help='BR=f(mA) or BR=f(mH)')
    parser.add_argument('-t', '--type', action='store', type=int, dest='type', default=1, choices=[1,2],
            help='2HDM Type ')
    parser.add_argument('-log', '--logy', action='store_true', required=False,
            help='log scale')
    parser.add_argument('-d', '--debug', action='store_true', required=False,
            help='will not ovwerite out.dat and out.log but instead you will get in format out_mH-{}_mA-{}_tb-{}_cba-{}.dat')
    options = parser.parse_args()

    cmssw = os.path.join(os.environ["CMSSW_BASE"], "src", "cp3_llbb/Calculators42HDM")
    
    use_default_scales = True
    
    if use_default_scales:
        muR4ggh = 0.5
        muF4ggh = 0.5
        muR4bbh = 1.
        muF4bbh = 0.25
    else:
        muR4ggh = mH/2
        muF4ggh = muR4ggh
        muR4bbh = (mh + mh + mb + mb__tilde__ )
        muF4bbh = muR4bbh
    
    plots = [
        'func_of_cba',
        'func_of_tb'
        ]
    if options.scan:
        plots += ['func_of_mA', 
                  'func_of_mH']
    MakeBRPlots(plots=plots)
    
