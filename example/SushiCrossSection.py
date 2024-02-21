#!/usr/bin/python
import sys
import math
import yaml 

inPath = '/home/ucl/cp3/kjaffel/bamboodev/ZA_FullAnalysis/ZAStatAnalysis'
sys.path.append(inPath)
import Constants as Constants

from cp3_llbb.Calculators42HDM.Calc2HDM import *

# some parameters aren't valid, be careful :  for example when m_H < m_h(125)
signal_grid  = Constants.get_SignalMassPoints('fullrun2', returnKeyMode=False, split_sig_reso_boo= False)
extra_masses = {'HToZA': [(200., 125.), (209.9, 37.34), (500., 250.), (442.63, 161.81), (510., 130.), (700., 200.), (750., 610.),
                          (780., 680.), (800., 140.), 
                          (240.0, 130.0), (300.0, 135.0), (220.0, 127.0), (670.0, 500.0), (550.0, 300.0)
                        ],
                'AToZH': [(240.0, 130.0), (300.0, 135.0), (700.0, 200.0), (250.0, 125.0), (750.0, 610.0), (500.0, 250.0), (800.0, 140.0), 
                          (200.0, 125.0), (510.0, 130.0), (780.0, 680.0), (220.0, 127.0), (670.0, 500.0), (550.0, 300.0)
                        ]
                } # you can always add more points here !

sqrts= 13000
type = 2
tb   = 1.5
cba  = 0.01
mh   = 125.
mb   = 4.92
mZ   = 9.118760e+01
mb__tilde__ = 4.92

dict_   = {}
dict_['tanbeta'] = tb
dict_['cos(b-a)']= cba
dict_['type']=type
dict_['model']='2HDM'

for process, masses_for_reg_mode in signal_grid.items():
    for decay, list_of_masses in masses_for_reg_mode['resolved'].items(): # that's okay resolved and boosted points are the same !
        dict_[decay] = {} 
        for m_heavy, m_light in list_of_masses+extra_masses[decay]:

            heavy   = decay[0]
            light   = decay[-1]
            mode    = decay[0]
            print(  m_heavy, m_light )
            
            dict_[decay]['M{}_{}_M{}_{}'.format(heavy, m_heavy, light, m_light)] = {}
            
            prod = 'gg%s'%mode if process == 'gg_fusion' else 'bb%s'%mode
            mhc  = max(m_heavy, m_light)
            m12  = math.sqrt(pow(mhc, 2) * tb / (1 + pow(tb, 2)))
            
            if prod == 'gg%s'%mode:
                muR = m_heavy/2
                muF = muR
            elif prod == 'bb%s'%mode:
                muR = (m_light + mZ + mb + mb__tilde__ )
                muF = muR

            if mode == 'H':
                mH = m_heavy 
                mA = m_light
            else:
                mH = m_light
                mA = m_heavy

            beta  =math.atan(tb)
            alpha =math.atan(tb)-math.acos(cba)
            sba   = math.sin(math.atan(tb)-alpha)
            # or simply : sba = math.sqrt(1 - pow(cba, 2))
            
            outputFile = "out_mH-{}_mA-{}_tb-{}_mode-{}.dat".format(mH, mA, tb, mode)
            
            calc = Calc2HDM(mode = mode, 
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
                           muR   = muR, 
                           muF   = muF)
            
            calc.setpdf('NNPDF31_nnlo_as_0118_nf_4_mc_hessian')
            calc.computeBR()
           
            if mode =='H':
                dict_[decay]['M{}_{}_M{}_{}'.format(heavy, m_heavy, light, m_light)]['branching-ratio']={ 
                        'HToZA': calc.HtoZABR, 
                        'ATobb': calc.AtobbBR}
            elif mode =='A':
                dict_[decay]['M{}_{}_M{}_{}'.format(heavy, m_heavy, light, m_light)]['branching-ratio']={ 
                        'AToZH': calc.AtoZHBR, 
                        'HTobb': calc.HtobbBR}


            cross_sections = calc.getXsecFromSusHi()

            if cross_sections:
                xsec_ggh = cross_sections['full']['ggh'][0]
                err_integration_ggh = cross_sections['full']['ggh'][1]

                xsec_bbh_lo   = cross_sections['split']['bbh_lo'][0]
                xsec_bbh_nlo  = cross_sections['split']['bbh_lo'][0]
                xsec_bbh_nnlo = cross_sections['split']['bbh_lo'][0]
                
                err_integration_bbh_lo   =  cross_sections['split']['bbh_lo'][1]
                err_integration_bbh_nlo  =  cross_sections['split']['bbh_nlo'][1]
                err_integration_bbh_nnlo =  cross_sections['split']['bbh_nnlo'][1]
                
                dict_[decay]['M{}_{}_M{}_{}_tb_{}'.format(heavy, m_heavy, light, m_light, tb)]['cross-section']= {
                    'gg%s'%mode: '{}  +/-  {}  pb'.format(xsec_ggh, err_integration_ggh),
                    'bb%s'%mode:
                        { 'LO'   : '{}  +/-  {}  pb'.format(xsec_bbh_lo, err_integration_bbh_lo),
                          'NLO'  : '{}  +/-  {}  pb'.format(xsec_bbh_nlo, err_integration_bbh_nlo),
                          'NNLO' : '{}  +/-  {}  pb'.format(xsec_bbh_nnlo, err_integration_bbh_nnlo)
                        }
                    }

print( dict_)

inFile = os.path.join(inPath, "data/sushi1.7.0-xsc_tanbeta-{}_2hdm-type{}.yml".format(tb, type))
with open(inFile, 'w') as _f:
    yaml.dump(dict_, _f, default_flow_style=False)

print( 'sushi file saved in :', inFile)   
