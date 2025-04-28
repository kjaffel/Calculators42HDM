#!/usr/bin/python

import math
import matplotlib.pyplot as plt

from cp3_llbb.Calculators42HDM.Calc2HDM import *

use_default_scales    = True 
split_to_subprocesses = True
doplot_crosssection = False

thdmtype = 1
mode = 'H' 
sqrts = 13600
tb = 10.
mh = 125.
mH = 600.
mA = 600.
mhc = max(mH, mA)
cba = 0.13
m12 = pow(mH,2)* pow(cba,2)/tb # or math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
beta=math.atan(tb)
alpha=math.atan(tb)-math.acos(cba)
sba = math.sin(math.atan(tb)-alpha) # or simply : sba = math.sqrt(1 - pow(cba, 2))
mb = 4.92
mb__tilde__ = 4.92

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

print( muR4ggh, muF4ggh, muR4bbh, muF4bbh )
outputFile = "forhh_out_mH-{}_mA-{}_tb-{}_type-{}.dat".format(mH, mA, tb , thdmtype)
test = Calc2HDM(mode  = mode,
                sqrts = sqrts,
                type  = thdmtype,
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
test.setpdf('NNPDF31_nnlo_as_0118_nf_4_mc_hessian')
test.computeBR()

cross_sections = test.getXsecFromSusHi()

xsec_ggh = cross_sections['full']['ggh'][0]
err_integration_ggh = cross_sections['full']['ggh'][1]

xsec_bbh_lo   = cross_sections['split']['bbh_lo'][0]
xsec_bbh_nlo  = cross_sections['split']['bbh_nlo'][0]
xsec_bbh_nnlo = cross_sections['split']['bbh_nnlo'][0]

err_integration_bbh_lo   =  cross_sections['split']['bbh_lo'][1]
err_integration_bbh_nlo  =  cross_sections['split']['bbh_nlo'][1]
err_integration_bbh_nnlo =  cross_sections['split']['bbh_nnlo'][1]


mH_list = [] 
xsectot_list = []
HtohhBR_list = []
htobbBR_list = [] 

if doplot_crosssection:
    while mH < 900:
        test.setmH(mH)
        test.setmA(mA)
        test.mhc = mhc
    
        xsec = test.getXsecFromSusHi()
        test.computeBR()
        
        mH_list.append(mH)
        xsectot_list.append(xsec_ggh*test.HtohhBR*test.htobbBR*test.htobbBR)
        mH+=100
    
    plt.plot(mH_list,xsectot_list, color='black')
    
    plt.ylabel(r'$\sigma$ ($pp \rightarrow H$) $\times$ BR($H \rightarrow hh \rightarrow bbbb$)')
    plt.xlabel(r'$m_H$ [GeV]')
    
    plt.savefig('test_hh.png')
