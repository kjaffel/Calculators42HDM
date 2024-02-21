import os, os.path
import argparse

from CP3SlurmUtils.Configuration import Configuration
from CP3SlurmUtils.SubmitWorker import SubmitWorker


def submitSushi(outdir):
    config = Configuration()
    config.sbatch_partition = 'Def'
    config.sbatch_qos = 'normal'
    config.cmsswDir = os.path.dirname(os.path.abspath(__file__))
    config.sbatch_chdir = os.path.join(outdir, 'slurm')
    config.sbatch_time = '48:59:00'
    config.sbatch_memPerCPU = '7000'
    config.apptainerImage = "/cvmfs/cp3.uclouvain.be/singularity/centos77-lscsoft"
    config.stageoutDir = config.sbatch_chdir
    config.inputParamsNames = ['cmssw', 'outdir']
    config.sbatch_additionalOptions=['--exclude=mb-sky[002,005-014,016-018,020],mb-ivy220,mb-ivy213,mb-ivy212,mb-ivy211']
    #config.stageoutFiles = ['*.json']
    #config.environmentType = 'cms'
    #config.inputSandboxContent = []
    #config.numJobs=1
    config.inputParams = []

    cmssw  = config.cmsswDir
    output = config.sbatch_chdir

    config.inputParams.append([cmssw, output])
    config.payload = \
    """
    python ${cmssw}/example/produce_xsec_for_ZA.py -o ${outdir} --cmssw ${cmssw}  
    """
    submitWorker = SubmitWorker(config, submit=True, yes=True, debug=True, quiet=True)
    submitWorker()

    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='slurm for 2hdmc and sushi', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-o", "--outdir", default='2hdmonslurm', required=True, help="output dir")
    options = parser.parse_args()

    submitSushi(outdir=options.outdir)
