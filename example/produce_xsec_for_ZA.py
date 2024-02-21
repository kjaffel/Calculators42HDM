#!/usr/bin/python
import os, os.path, sys
import argparse
import json
import math

from cp3_llbb.Calculators42HDM.Calc2HDM import *


def float_to_str(x, digits=2):
    tmp = ':.{:d}f'.format(digits)
    tmp = ('{' + tmp + '}').format(x)
    return tmp.replace('.', 'p')



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='finner grid', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-o", "--outdir", default='2hdmonslurm', required=True, help="output dir")
    parser.add_argument("--cmssw", required=False, default=os.path.dirname(os.path.abspath(__file__)), help="this script need to be run inside cmssw dir")

    options = parser.parse_args()

    type = 2
    sqrts = 13000
    mh = 125
    tb = 1.5
    cba = 0.01
    sba = math.sqrt(1 - pow(cba, 2))
    mZ = 91.1876
    ZtollBR =  3.3658 * 2 / 100. # no taus
    mb = 4.92 # mb(OS) pole mass
    mb__tilde__ = 4.92 # mb~
    cmssw= os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
    
    results = {
        'mH': [],
        'mA': [],
        'sigma': [],
        'sigma_ggh': [],
        'sigma_bbh_lo': [],
        'sigma_bbh_nlo': [],
        'sigma_bbh_nnlo': [],
        'sigma_errIntegration_ggh': [],
        'sigma_errIntegration_bbh_lo': [],
        'sigma_errIntegration_bbh_nlo': [],
        'sigma_errIntegration_bbh_nnlo': [],
        'BR':[],
        'AtoZHBR':[],
        'HtoZABR':[],
        'HtobbBR':[],
        'AtobbBR':[],
        }

    try:
        stageOutF = os.path.join(cmssw, options.outdir, 'outputs')
        if not os.path.exists(stageOutF):
            os.makedirs(stageOutF)
    except OSError as err:
        print(err)
    
    
    for mH in xrange(125, 1010, 10):
        for mA in xrange(30, 1010, 10):
            if mA >= mH:
                mode = 'A'
                light   = 'H'
                heavy   = 'A'
                m_heavy = mA
                m_light = mH
            else:
                mode = 'H'
                light   = 'A'
                heavy   = 'H'
                m_heavy = mH
                m_light = mA
            
            mhc = max(mH, mA)
            m12 = math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
            
            muR4ggh = mH/2
            muF4ggh = muR4ggh
            muR4bbh = (mA + mZ + mb + mb__tilde__ )
            muF4bbh = muR4bbh
            File = "out_m{}-{}_m{}-{}_tb-{}_cosba-{}_mode-{}.dat".format(heavy, m_heavy, light, m_light, tb, round(cba, 2), mode)
            outputFile = os.path.join(options.outdir, 'outputs', File)

            print ('mode::', mode)
            x = Calc2HDM(mode = mode, 
                         sqrts = sqrts, 
                         type = type, 
                         tb = tb, 
                         m12 = m12, 
                         mh = mh, 
                         mH = mH, 
                         mA = mA, 
                         mhc = mhc, 
                         sba = sba, 
                         outputFile = outputFile, 
                         muR4ggh = muR4ggh, 
                         muF4ggh = muF4ggh,
                         muR4bbh = muR4bbh, 
                         muF4bbh = muF4bbh 
                         )
            
            os.chdir(cmssw)
            x.setpdf('NNPDF31_nnlo_as_0118_nf_4_mc_hessian')
            x.computeBR()
            cross_sections = x.getXsecFromSusHi()
            
            if cross_sections:
                xsec_ggh = cross_sections['full']['ggh'][0]
                err_integration_ggh = cross_sections['full']['ggh'][1]
        
                xsec_bbh_lo   = cross_sections['split']['bbh_lo'][0]
                xsec_bbh_nlo  = cross_sections['split']['bbh_lo'][0]
                xsec_bbh_nnlo = cross_sections['split']['bbh_lo'][0]

                err_integration_bbh_lo   =  cross_sections['split']['bbh_lo'][1]
                err_integration_bbh_nlo  =  cross_sections['split']['bbh_nlo'][1]
                err_integration_bbh_nnlo =  cross_sections['split']['bbh_nnlo'][1]
    
                if mode == 'H': print ("xsec ggh & bbh (lo), HtoZABR, AtobbBR ", xsec_bbh_lo, xsec_bbh_lo, x.HtoZABR, x.AtobbBR)
                else: print ("xsec ggh & bbh (lo), AtoZHBR, HtobbBR ", xsec_bbh_lo, xsec_bbh_lo, x.AtoZHBR, x.HtobbBR)
                
                if mode == 'H':
                    if any ( x is None for x in [x.HtoZABR, x.AtobbBR]):
                        continue
                    results['BR'].append(x.HtoZABR * x.AtobbBR * ZtollBR)
                    results['HtoZABR'].append(x.HtoZABR)
                    results['AtobbBR'].append(x.AtobbBR)
                else:
                    if any ( x is None for x in [x.AtoZHBR, x.HtobbBR]):
                        continue
                    results['BR'].append(x.AtoZHBR * x.HtobbBR * ZtollBR)
                    results['AtoZHBR'].append(x.AtoZHBR)
                    results['HtobbBR'].append(x.HtobbBR)
                
                results['mH'].append(mH)
                results['mA'].append(mA)
                results['sigma_ggh'].append(xsec_ggh)
                results['sigma_bbh_lo'].append(xsec_bbh_lo)
                results['sigma_bbh_nlo'].append(xsec_bbh_nlo)
                results['sigma_bbh_nnlo'].append(xsec_bbh_nnlo)
                results['sigma_errIntegration_ggh'].append(err_integration_ggh)
                results['sigma_errIntegration_bbh_lo'].append(err_integration_bbh_lo)
                results['sigma_errIntegration_bbh_nlo'].append(err_integration_bbh_nlo)
                results['sigma_errIntegration_bbh_nnlo'].append(err_integration_bbh_nnlo)
    
    filename = 'sigmaBR_HZA_type-2_tb-%s_cba-%s_mirroring__khawla_ver1.json' % (float_to_str(tb, 1), float_to_str(cba, 2))
    
    with open(os.path.join(options.outdir, filename), 'w+') as f:
        json.dump(results, f)
