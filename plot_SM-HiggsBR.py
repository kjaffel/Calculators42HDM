#!/usr/bin/python

import math
import argparse
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

from cp3_llbb.Calculators42HDM.Calc2HDM import *
from labellines import *


def get_options():
    parser = argparse.ArgumentParser(description='Compare gg fusion and b-associated production xsc as function of tan beta ')
    parser.add_argument('-lazy', '--lazy', action='store_true', required=False, 
            help='This will allow you to do fast test, with few files saved in data/ ')
    options = parser.parse_args()
    return options

def which_points(grid):
    grid['fullsim'] = [
        #( MH, MA )
         ( 800, 660),
        #( 200, 50), 
        #( 200, 100),
        #( 250, 50), ( 250, 100),
        #( 300, 50), ( 300, 100),
        #( 300, 200),
        #( 500, 50), 
        #( 500, 100), ( 500, 200), 
        #( 500, 300), 
        #( 500, 400),
        #( 650, 50),
        #( 800, 50), ( 800, 100), 
        #( 800, 200),              
        #( 800, 400),              
        #( 800, 700),
        #(1000, 50),              
        #(1000, 200),                           
        #(1000, 500),
        ]
    return grid 

def check_sm_higgs_BRin2hdm(H, A):
    mH = H
    mA = A
    mhc = max(mH, mA)
    sqrts = 13000
    type = 1
    tb = 0.5
    mh = 125.
    mZ = 91.1876
    cba = 0.2
    mode = 'H'
    tb_list = [] 
    results = {
                'BRhtoss': [],
                'BRhtocc': [],
                'BRhtobb': [],
                'BRhtoee': [],
                'BRhtomumu': [],
                'BRhtotautau': [],
                'BRhtogg': [],
                'BRhtoZZ': [],
                'BRhtoWW': [],
                'BRhtoZga': [],
                'BRhtogluglu': [],}
    
    while tb < 60.:
        m12 = math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
        tb_list.append(tb)
        beta=math.atan(tb)
        alpha=math.atan(tb)-math.acos(cba)
        sba = math.sin(math.atan(tb)-alpha)
        outputFile = "forhh_out_mH-{}_mA-{}_tb-{}_type-{}_mode-{}_cba-{}.dat".format(mH, mA, tb , type, mode, cba)
        
        test = Calc2HDM(mode = 'H', 
                        sqrts = sqrts, 
                        type = type, 
                        tb = tb, 
                        m12 = m12, 
                        mh = mh, 
                        mH = mH, 
                        mA = mA, 
                        mhc = mhc, 
                        sba = sba, 
                        outputFile = outputFile)
        
        test.setm12(m12)
        test.settb(tb)
        test.setsba(sba)
        test.setmA(mA)
        test.setmH(mH)
        test.computeBR()
        
        results['BRhtoss'].append(test.htossBR)
        results['BRhtocc'].append(test.htoccBR)
        results['BRhtobb'].append(test.htobbBR)
        results['BRhtoee'].append(test.htoeeBR)
        results['BRhtomumu'].append(test.htomumuBR)
        results['BRhtotautau'].append(test.htotautauBR)
        results['BRhtogg'].append(test.htoggBR) # gamma-gamma
        results['BRhtoZZ'].append(test.htoZZBR)
        results['BRhtoWW'].append(test.htoWWBR)
        results['BRhtoZga'].append(test.htoZgaBR)
        results['BRhtogluglu'].append(test.htoglugluBR)

        tb+=1.5
        
    return results, tb_list

options = get_options()
Colors=['forestgreen', 'pink', 'crimson', 'magenta', 'indigo', 'limegreen', 'blueviolet', 'plum', 'purple', 'hotpink', 'mediumseagreen', 'springgreen', 'aquamarine', 'turquoise', 'aqua', 'mediumslateblue', 'orchid', 'deeppink', 'darkturquoise', 'teal', 'mediumslateblue']
grid = {}
grid = which_points(grid)

for decayto in ['2fermions', '2gaugebosons']:
    idx=0
    for H, A in (grid['fullsim']):
        results,tb_list = check_sm_higgs_BRin2hdm(H, A)
        fig= plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        if decayto =='2fermions':
            plt.plot(tb_list,results['BRhtobb'], color=Colors[idx], marker='o', label=r'b$\bar{b}$')
            #plt.text( 5., 0.01, r'$H \rightarrow bb$', color= Colors[idx]) 
            plt.plot(tb_list,results['BRhtotautau'], color=Colors[idx+1], marker='o', label=r'$\tau \tau$')
            plt.plot(tb_list,results['BRhtomumu'], color=Colors[idx+2], marker='o', label=r'$\mu^+\mu^-$')
            plt.plot(tb_list,results['BRhtocc'], color=Colors[idx+3], marker='o', label=r'c$\bar{c}$')
            
            # from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR
            # H --> 2 fermions
            # M(GeV)       BR ->bb          THU             PU(mq)         PU(alphas)
            # 125.00    5.824E-01   +0.65   -0.65   +0.72   -0.74   +0.78   -0.80 
            ax.axhspan((5.824E-01- 0.65*5.824E-01/100.) , (5.824E-01+ 0.65*5.824E-01/100.) , color=Colors[idx])
            #            h -> tautau
            # 125.00    6.272E-02     +1.17   -1.16   +0.98   -0.99   +0.62   -0.62
            ax.axhspan((6.272E-02- 1.160*6.272E-02/100.) , (6.272E-02+ 1.17*6.272E-02/100.) , color=Colors[idx+1])
            #           h -> mumu 
            # 125.00    2.176E-04     +1.23   -1.23   +0.97   -0.99   +0.59   -0.64   
            ax.axhspan((2.176E-04- 1.23*2.176E-04/100.) , (2.176E-04+ 1.23*2.176E-04/100.) , color=Colors[idx+2])
            #           h -> cc    
            #           2.891E-02     +1.20   -1.20   +5.26   -0.98   +1.25   -1.25 
            ax.axhspan((2.891E-02- 1.20*2.891E-02/100.) , (2.891E-02+ 1.20*2.891E-02/100.) , color=Colors[idx+3])
            
            labelLines(plt.gca().get_lines(),align=True,fontsize=10)
        
        else:
            plt.plot(tb_list,results['BRhtogg'], color=Colors[idx], marker='o', label=r'$\gamma \gamma$')
            plt.plot(tb_list,results['BRhtoZZ'], color=Colors[idx+1], marker='o', label=r'$ZZ$')
            plt.plot(tb_list,results['BRhtoWW'], color=Colors[idx+2], marker='o', label=r'$WW$')
            plt.plot(tb_list,results['BRhtoZga'], color=Colors[idx+3], marker='o', label=r'Z$\gamma$')
            plt.plot(tb_list,results['BRhtogluglu'], color=Colors[idx+4], marker='o', label=r'$gg$')
        
            # H --> 2 gauge bosons
            # 2.270E-03     +1.73   -1.72   +0.93   -0.99   +0.61   -0.62
            ax.axhspan((2.270E-03- 1.73*2.270E-03/100.) , (2.270E-03+ 1.72*2.270E-03/100.) , color=Colors[idx])
            # 2.619E-02     +0.99   -0.99   +0.99   -0.98   +0.66   -0.63
            ax.axhspan((2.619E-02- 0.99*2.619E-02/100.) , (2.619E-02+ 0.99*2.619E-02/100.) , color=Colors[idx+1])
            # 2.137E-01   +0.99   -0.99   +0.99   -0.98   +0.66   -0.63
            ax.axhspan((2.137E-01- 0.99*2.137E-01/100.) , (2.137E-01+ 0.99*2.137E-01/100.) , color=Colors[idx+2])
            # 1.533E-03   +5.71   -5.71   +0.98   -1.01   +0.58   -0.65 
            ax.axhspan((1.533E-03- 5.71*1.533E-03/100.) , (1.533E-03+ 5.71*1.533E-03/100.) , color=Colors[idx+3])
            # 8.187E-02     +3.40   -3.41   +1.12   -1.13   +3.69   -3.61
            ax.axhspan((8.187E-02- 3.40*8.187E-02/100.) , (8.187E-02+ 3.41*8.187E-02/100.) , color=Colors[idx+4])
            
            labelLines(plt.gca().get_lines(),align=True,fontsize=10)

        idx +=1

        plt.ylabel(r'BR(SM-Higgs vs SM-like $\rightarrow$ XX)')
        plt.xlabel(r'$tan\beta$')
        plt.yscale('log')
        #plt.xscale('log')
        plt.title(r'2HDM-typeI: $M_{H^\pm}=M_{H}$, $cos(\beta-\alpha)$= 0.01, mh= 125. GeV, SM-like (o), SM-Higgs (-)', fontsize=10.)
        plt.xlim(min(tb_list), max(tb_list))
        plt.ylim(1E-7, 1)
        plt.legend()
        #plt.grid()
        fig.savefig('smhiggs_BR{}_comparetoBRfrom-2hdmt2_asfunc-tb_Benchmark{}_MH-{}_MA-{}.png'.format(decayto, idx, H, A))
        plt.gcf().clear()
