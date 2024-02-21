#!/bin/env python
import math
from cp3_llbb.Calculators42HDM.Calc2HDM import *
import json

mode = 'H' # means A-> ZH  or 'H' means H ->ZA
sqrts = 13000
type = 2
mh = 125
tb = 1.
sba = 0.99 # math.sqrt(1 - pow(cba, 2))
cba = math.sqrt(1 + pow(sba, 2))
mZ = 91.1876
ZtollBR =  3.3658 * 2 / 100. # no taus

mb = 4.92 # mb(OS) pole mass
mb__tilde__ = 4.92 # mb~
MZ= 9.118760e+01

results = {
    'mH': [],
    'mA': [],
    'sigma': [],
    'sigma_errIntegration': [],
    'sigma_err_muRm': [],
    'sigma_err_muRp': [],
    'BR': []
    }

def float_to_str(x, digits=2):
    tmp = ':.{:d}f'.format(digits)
    tmp = ('{' + tmp + '}').format(x)
    return tmp.replace('.', 'p')

for mA, mH in [
#    (1000,600),
#    (2100,2000),
#    (900,350),
#    (500,400),
#    (2100,1000),
#    (1500,1400),
#    (1500,1000),
    (2100,1600),
#    (500,350),
#    (700,350),
#    (1200,850),
#    (1800,1600),
#    (500,370),
#    (700,370),
#    (900,370),
#    (700,400),
#    (900,400),
#    (1200,400),
#    (1500,400),
#    (1800,400),
#    (2100,400),
#    (800,600),
#    (800,650),
#    (1000,700),
#    (1000,800),
    (1200,1000) ]:
    
    mhc = max(mH, mA)
    m12 = math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
    
    muR4ggh = mH/2
    muF4ggh = muR4ggh
    muR4bbh = mH/2
    muF4bbh = (mH+2*mb)/4
    
    #outputFile = "2hdmc_results/2hdmc1.8.0_mA-{}_mH-{}.dat".format(mA, mH)
    outputFile = "2hdmc1.8.0_mA-{}_mH-{}.dat".format(mA, mH)
    x = Calc2HDM(mode = mode, sqrts = sqrts, type = type, tb = tb, m12 = m12, mh = mh, mH = mH, mA = mA, mhc = mhc, sba = sba, outputFile = outputFile, muR4ggh = muR4ggh, muF4ggh = muF4ggh,  muR4bbh=muR4bbh, muF4bbh=muF4bbh)
    x.setpdf('NNPDF31_nnlo_as_0118_nf_4_mc_hessian')
    x.computeBR()
    
    #xsec_ggH, err_integration_ggH, err_muRm_ggH, err_muRp_ggH, xsec_bbH, err_integration_bbH, mb_MSscheme_muR= x.getXsecFromSusHi()
    #print ( " xsec( gg-fusion), AtoZHBR, HtobbBR ", xsec_ggH, x.AtoZHBR, x.HtottBR )
    
    #results['mH'].append(mH)
    #results['mA'].append(mA)
    #results['sigma'].append(xsec_ggH)
    #results['sigma_errIntegration'].append(err_integration_ggH)
    #results['sigma_err_muRm'].append(err_muRm_ggH)
    #results['sigma_err_muRp'].append(err_muRp_ggH)
    #results['BR'].append(x.AtoZHBR * x.HtottBR * ZtollBR)

#with open('AToZH_xsc_br_results.json', 'w+') as f:
#    json.dump(results, f)
