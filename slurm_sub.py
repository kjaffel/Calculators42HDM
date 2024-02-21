import os, os.path, sys
import math
import json
import yaml
import glob 
import argparse
import subprocess
import numpy as np

from CP3SlurmUtils.Configuration import Configuration
from CP3SlurmUtils.SubmitWorker import SubmitWorker



def generate_tb_ranges(start, end, step, num_ranges):
    """
    Generate a list of ranges for the parameter tb.
    Args:
    - start (float): The starting value of tb.
    - end (float): The ending value of tb.
    - step (float): The step size for tb.
    - num_ranges (int): The number of ranges to generate.
    Returns:
    - List of lists: A list containing the specified number of ranges for tb.
    """
    tb_ranges = []
    for i in range(num_ranges):
        tb_range = [start + i * step, start + (i + 1) * step]
        tb_ranges.append(tb_range)
    return tb_ranges


def postprocessing(mode, mH, mA, thdm, outdir, do):
    if do == 'tanb_vs_mass':
        cba    = run_opts[do]['cba'] 
        mTofix = run_opts[do]['m_tofix']
        _type  = run_opts[do]['type']

        full_dict = {'HToZA': {}, 'AToZH': {} }
        full_dict['cos(b-a)']= cba
        full_dict['type']=_type
        full_dict['model']='2HDM'

        for jsF in glob.glob( os.path.join(outdir, 'outputs', '*', "2hdmc1.8.0-br_cba-{}_mAorH-{}_2hdm-type{}.yml".format(cba, mTofix, _type))):
                
            if not os.path.isfile(jsF):
                continue

            with open(jsF) as f_:
                data = yaml.safe_load(f_)
            
            full_dict['HToZA'].update(data['HToZA'])
            full_dict['AToZH'].update(data['AToZH'])

        TotFile = os.path.join(outdir, "2hdmc1.8.0-br_cba-{}_mAorH-{}_2hdm-type{}.yml".format(cba, mTofix, _type))
        with open(TotFile, 'w') as _f:
            yaml.dump(full_dict, _f, default_flow_style=False)
        print( 'Done finalizing tasks, result file can be found in :: %s'%TotFile)
    
    elif do == 'tanb_vs_cosba':
        results = {
            'tb': [],
            'cba': [],
            'sigma_ggh': [],
            'sigma_errIntegration_ggh': [],
            'sigma_bbh_lo': [],
            'sigma_errIntegration_bbh_lo': [],
            'sigma_bbh_nlo': [],
            'sigma_errIntegration_bbh_nlo': [],
            'TotBR': [],
            'HtoZABR': [],
            'AtobbBR': [],
            'AtoZHBR': [],
            'HtobbBR': [],
        }
        
        for n, (ki, kf) in enumerate(generate_tb_ranges(run_opts[do]['tb'][0], run_opts[do]['tb'][-1], 5.0, 10)):
            for i, tb in enumerate(np.arange(ki, kf, 0.5)):
                for j, ba in enumerate(run_opts[do]['ba']):

                    cba = math.cos(ba) 
                    
                    filename  = 'sigmaBR_%sZ%s_type-2_M%s-%s_M%s-%s_tb-%s_cba-%s.json' %(
                            heavy, light, heavy, str(round(m_heavy,2)), light, str(round(m_light,2)), str(round(tb,2)), str(round(cba, 5)))
                    
                    jsF = os.path.join(outdir, 'results', filename)
                    if not os.path.exists(jsF):
                        continue
                    
                    print('processing ::', jsF)
                    with open(jsF) as f:
                        cfg = json.load(f)
                    
                    if not cfg:
                        continue
                    
                    for k in cfg.keys():
                        results[k.decode("utf-8")] += cfg[k]
        
        jsF_result = os.path.join(outdir, 'sigmaBR_%sZ%s_type-2_M%s-%s_M%s-%s_tb_cba.json'%(
                    heavy, light, heavy, str(round(m_heavy,2)), light, str(round(m_light,2)))
                    )    
        
        print( len(results['tb']), len(results['cba']))
        print( max(results['tb']), max(results['cba']) ,  min(results['tb']), min(results['cba']))
        
        with open(jsF_result, 'w+') as f:
           json.dump(results, f)
        print( 'result file is saved in ::' , jsF_result)
    return 


def SlurmRun2HDMCalculator(output, mode, mH, mA, finalize=False, test=False, do=''):
    config = Configuration()
    config.sbatch_partition = 'cp3'
    config.sbatch_qos = 'cp3'
    config.cmsswDir = os.path.dirname(os.path.abspath(__file__))
    config.sbatch_chdir = os.path.join(output, 'slurm', thdm)
    config.stageoutDir = config.sbatch_chdir
    config.sbatch_time = '28:59:00'
    config.sbatch_memPerCPU = '7000'
    config.apptainerImage  = "/cvmfs/cp3.uclouvain.be/singularity/centos77-lscsoft"
    #config.stageoutFiles = ['*.json']
    #config.environmentType = 'cms'
    #config.inputSandboxContent = []
    #config.numJobs=1
    config.inputParams = []
    cmssw  = config.cmsswDir
    outdir = config.sbatch_chdir
    
    
    if finalize:
        postprocessing(mode, mH, mA, thdm, outdir, do)
    
    else:
        test_executed = False
        if do == 'tanb_vs_cosba':
            config.inputParamsNames = ["cmssw", "outdir", "mode", "mH", "mA", "tb", "ba"]
            
            for n, (ki, kf) in enumerate(generate_tb_ranges(run_opts[do]['tb'][0], run_opts[do]['tb'][-1], 5.0, 10)):# to avoid running into max job submission allowed
                for i, tb in enumerate(np.arange(ki, kf, 0.5)):
                    for j, ba in enumerate(run_opts[do]['ba']):
                        
                        if test and test_executed and len(config.inputParams) ==1:
                            break  # Exit the innermost loop after running the test once
                        config.inputParams.append([cmssw, outdir, mode, mH, mA, tb, ba])
                        test_executed = True
                config.payload = \
                """
                python ${cmssw}/example/produce_xsec_for_ZA_tanb_cosba.py -o ${outdir} \\
                                                                          --cmssw ${cmssw} \\
                                                                          --mode ${mode} \\
                                                                          --mH ${mH} \\
                                                                          --mA ${mA} \\
                                                                          --tb ${tb} \\
                                                                          --ba ${ba}
                """
                submitWorker = SubmitWorker(config, submit=True, yes=True, debug=True, quiet=True)
                submitWorker()

        elif do =='tanb_vs_mass':
            config.inputParamsNames = ["cmssw", "outdir", "mpave", "mfix", "tanbeta", "slurm", "jobid"]
            
            mfix = run_opts[do]['m_tofix']
            jobid = 1
            for i, mpave in enumerate(run_opts[do]['m_topave']):
                for j, tanbeta in enumerate(run_opts[do]['tb']):
                    
                    if test and test_executed and len(config.inputParams) ==1:
                        break  # Exit the innermost loop after running the test once
                    config.inputParams.append([cmssw, outdir, mpave, mfix, tanbeta, 'true', jobid])
                    test_executed = True
                    jobid +=1
            config.payload = \
            """
            echo ${jobid} ${tanbeta} ${mfix} ${mpave} ${cmssw} ${outdir}  
            python ${cmssw}/example/runsushi_for_2hdm_onslurm.py -o ${outdir} \\
                                                                 --cmssw ${cmssw} \\
                                                                 --mpave ${mpave} \\
                                                                 --mfix ${mfix} \\
                                                                 --tanbeta ${tanbeta} \\
                                                                 --jobid ${jobid} \\
                                                                 --slurm
            """
            submitWorker = SubmitWorker(config, submit=True, yes=True, debug=True, quiet=True)
            submitWorker()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='slurm for 2hdmc and sushi', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument("-o", "--output"  , default=None, required=True, help="output dir")
    parser.add_argument("--mode", type=str, choices=['H', 'A'], default ='H', required=False,
                help="This will decide the mode (H-> ZA or A -> ZH) also the renormalisation and factorisation scale")
    parser.add_argument('--mH', type=float, default=500., help='mass of the Higgs bosons H')
    parser.add_argument('--mA', type=float, default=300., help='mass of the pseudo-scalar A')
    parser.add_argument("--finalize"      , default=False, action='store_true', 
                help='This is a post-processing step to combine all results in one json File')
    parser.add_argument("--test"          , default=False, action='store_true', help='')
    
    options = parser.parse_args()
    
    
    Extra  = [ 50., 100., 200., 250., 300., 400., 450., 1000.]
    Extra += [ 55.16, 566.51, 160.17, 87.1, 298.97, 186.51, 118.11, 137.54, 40.68, 47.37, 
               779.83, 350.77, 664.66, 74.8, 34.93, 254.82, 101.43, 217.19, 482.85, 411.54, 64.24, 30.0]
    
    run_opts = { 'tanb_vs_mass' : { 
                        'tb'      : np.arange(0.05, 50.1, 0.5),
                        'cba'     : 0.01,
                        'm_topave': sorted(Extra + np.arange(10., 1000., 50.).tolist()),
                        'm_tofix' : 379.,
                        'type'    : 2 },

                 'tanb_vs_cosba': { 
                        'tb': np.arange(0.01, 50.1, 0.5), 
                        'ba': np.arange(0.13, math.pi, 0.03) },
                 }

    if options.mode == 'H':
        light   = 'A'
        heavy   = 'H'
        m_heavy = options.mH
        m_light = options.mA
    elif options.mode == 'A':
        light   = 'H'
        heavy   = 'A'
        m_heavy = options.mA
        m_light = options.mH

    thdm = '%sToZ%s'%(heavy, light)
    for do_what in ['tanb_vs_mass']: #'tanb_vs_mass', 'tanb_vs_cosba']:

        if do_what == 'tanb_vs_mass': thdm ='' # will run over both

        SlurmRun2HDMCalculator(output  = options.output, 
                               mode    = options.mode, 
                               mH      = options.mH, 
                               mA      = options.mA, 
                               finalize= options.finalize, 
                               test    = options.test,
                               do      = do_what )
