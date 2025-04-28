#!/usr/bin/python

import math
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np


from cp3_llbb.Calculators42HDM.Calc2HDM import *

# ggH (cosb-a, tanbeta) =  (-0.925, 3.103)
# bbH (cosb-a, tanbeta) = (-0.278, 13.003) 

use_default_scales    = False
split_to_subprocesses = True
process = 'H'
mode = 'H'
sqrts = 13000
type = 1
tb = 5.578
mh = 125.
mH = 335.4
mA = 82.14 #99.61
ZtollBR = 0.067264
mZ = 9.118760e+01
mhc = max(mH, mA)
m12 = math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
beta=math.atan(tb)
cba = -0.02
alpha=math.atan(tb)-math.acos(cba)
sba = math.sin(math.atan(tb)-alpha)
# or simply : sba = math.sqrt(1 - pow(cba, 2))
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
    muR4bbh = (mA + mZ + mb + mb__tilde__ )
    muF4bbh = muR4bbh

print( muR4ggh, muF4ggh, muR4bbh, muF4bbh )

outputFile = "forstep3_out_mH-{}_mA-{}_tb-{}_type-{}.dat".format(mH, mA, tb , type)

test = Calc2HDM(mode  = mode, 
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
                            
print( '(mH, mA)= ({}, {}) GeV , tb= {}, cosb-a= {}'.format( mH, mA, tb, cba))
print( 'cross-section ggH: {} +/- {} pb'.format(xsec_ggh, err_integration_ggh))
print( 'cross-section bbH: {}, {}, {} pb -> LO, NLO, NNLO'.format(xsec_bbh_lo, xsec_bbh_nlo, xsec_bbh_nnlo))
print( '       Uncer. bbH: {}, {}, {} pb -> LO, NLO, NNLO'.format(err_integration_bbh_lo, err_integration_bbh_nlo, err_integration_bbh_nnlo))

print( ' BR (H -> ZA): {}'.format(test.HtoZABR))
print( ' BR (A -> bb): {}'.format(test.AtobbBR))
print( ' BR (Z -> ll): {}'.format(ZtollBR))
print( ' BR (A -> bb)* BR (H -> ZA)* BR (Z -> ll): {}'.format(test.AtobbBR*test.HtoZABR*ZtollBR))
print( "*"*80)


mA_list = [] 
xsec_ggh_list = []
xsec_bbh_list = []
subprocess_gg_nlo_list = []
subprocess_qg_nlo_list = []
subprocess_qq_nlo_list = []

HtoZABR_list = []
AtobbBR_list = [] 

fig= plt.figure(figsize=(8,6))
ax = fig.add_subplot(111)

while mA < mH-mZ:
    test.setmA(mA)
    test.computeBR()
    mA_list.append(mA)
    
    xsec_ggh_X_BR = xsec_ggh*test.HtoZABR*test.AtobbBR
    xsec_bbh_X_BR = xsec_bbh*test.HtoZABR*test.AtobbBR
    
    xsec_bbh_list.append(xsec_bbh_X_BR)
    xsec_ggh_list.append(xsec_ggh_X_BR)
    
    subprocess_gg_nlo_X_BR = subprocess_gg_nlo*test.HtoZABR*test.AtobbBR
    subprocess_qg_nlo_X_BR = subprocess_qg_nlo*test.HtoZABR*test.AtobbBR
    subprocess_qq_nlo_X_BR = subprocess_qq_nlo*test.HtoZABR*test.AtobbBR
    
    subprocess_gg_nlo_list.append(subprocess_gg_nlo_X_BR)
    subprocess_qg_nlo_list.append(subprocess_qg_nlo_X_BR)
    subprocess_qq_nlo_list.append(subprocess_qq_nlo_X_BR)
    
    HtoZABR_list.append(test.HtoZABR)
    AtobbBR_list.append(test.AtobbBR)
    
    mA+=2

if split_to_subprocesses :
    plt.plot(mA_list,subprocess_gg_nlo_list, color='mediumturquoise', marker='o', label= 'gg [NLO]' )
    plt.plot(mA_list,subprocess_qg_nlo_list, color='indigo', marker='o', label= 'qg [NLO]' )
    plt.plot(mA_list,subprocess_qq_nlo_list, color='purple', marker='o', label= 'qq [NLO]' )
    plt.plot(mA_list,xsec_ggh_list, color='black', marker='o', label= 'gg-fusion from SUSHI [NLO]' )
    plt.plot(mA_list,xsec_bbh_list, color='gold', marker='o', label='bb-associated production [NNLO]')
else:
    plt.plot(mA_list,xsec_ggh_list, color='black', marker='o', label= 'gg-fusion from SUSHI' )
    plt.plot(mA_list,xsec_bbh_list, color='gold', marker='o', label='bb-associated production')

plt.plot(mA_list,HtoZABR_list, color='red', linestyle='dashed', linewidth=1.5, label=r'$BR(H\rightarrow ZA)$')
plt.plot(mA_list,AtobbBR_list, color='b', linestyle='dashed', linewidth=1.5, label=r'$BR(A\rightarrow bb)$')
ax.axvspan((172.76- 1)*2 , (172.76+ 1)*2 , alpha=0.5, color='lightgrey', label=r'$2 x Top_{mass} \pm 1\sigma $')
plt.ylabel(r'$\sigma* BR(H\rightarrow ZA)* BR(A\rightarrow bb)[pb]$')
plt.xlabel(r'$M_{A} [GeV]$')
plt.yscale('log')
plt.title(r'$2HDM-type%s: M_{H}=%s GeV, M_{H^+}=M_{H}, tan\beta= %s, cos(\beta-\alpha)= 0.01, mh= 125. GeV$'%(('I'if type==1 else('II')), mH, tb), fontsize=10.)
plt.xlim(min(mA_list), max(mA_list))
plt.legend()
plt.grid()
fig.savefig('test_2hdm-type%s_tb-%s_mH-%s_splitProd-%s.png'%(type, tb, mH, split_to_subprocesses))
plt.gcf().clear()
