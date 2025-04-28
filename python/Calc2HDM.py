#!/usr/bin/python

import os
import sys
import subprocess
import math
import logging
import colorlog
    
from colorlog import ColoredFormatter

LOG_LEVEL = logging.DEBUG
stream    = logging.StreamHandler()
stream.setLevel(LOG_LEVEL)
logger    = logging.getLogger("Calculators42HDM")
logger.setLevel(LOG_LEVEL)
logger.addHandler(stream)
formatter = ColoredFormatter(
"%(log_color)s%(levelname)-8s%(reset)s %(log_color)s%(message)s",
datefmt=None,
reset=True,
log_colors={
        'DEBUG':    'cyan',
        'INFO':     'green',
        'WARNING':  'blue',
        'ERROR':    'red',
        'CRITICAL': 'red',
        },
secondary_log_colors={},
style='%')
stream.setFormatter(formatter)


class Calc2HDM:
    def __init__(self, mode = 'H', sqrts = 13000,  muR4ggh = 0.5, muF4ggh = 0.5, 
            muR4bbh = 1., muF4bbh = 0.25, type = 2, tb = 1, m12 = 0, mh = 125., 
            mH = 350., mA = 400., mhc = 420., sba = 0.99, outputFile = "out.dat"):

        # 11 = light Higgs (h), 12 = heavy Higgs (H), 21 = pseudoscalar (A)
        if mode == 'h' :
            self.mode = 11
        elif mode == 'H' :
            self.mode = 12
        elif mode == 'A' :
            self.mode = 21
        else :
            print ("Wrong mode selected! Please, use h, H or A." )

        self.tb     = tb
        self.m12    = m12
        self.m12_2  = m12*m12
        self.mh     = mh
        self.mH     = mH
        self.mA     = mA
        self.mhc    = mhc
        self.sba    = sba
        self.type   = type
        self.sqrts  = sqrts
        self.muR4ggh= muR4ggh
        self.muF4ggh= muF4ggh
        self.muR4bbh= muR4bbh
        self.muF4bbh= muF4bbh
        self.BRcomputed = 0
        self.outputFile = outputFile
        self.pdf    = "NNPDF31_nnlo_as_0118_nf_4_mc_hessian"
        self.HtohhBR =None
        self.HtoZABR =None
        self.AtoZHBR =None
        self.HtobbBR =None
        self.AtobbBR =None
        self.htobbBR =None
        self.HtottBR =None

    def __str__(self):
        return """
            tb= %.2f
            m12= %.2f
            mh= %.2f
            mH= %.2f
            mA= %.2f
            mhc= %.2f
            sba= %.2f
            type= %i
            outputFile= %s
            sqrts= %.2f
            muR4ggh= %.2f
            muF4ggh= %.2f
            muR4bbh= %.2f
            muF4bbh= %.2f
            BRcomputed= %.2f
            pdf= %s""" % (self.tb, 
                          self.m12, 
                          self.mh, 
                          self.mH, 
                          self.mA,
                          self.mhc, 
                          self.sba, 
                          self.type, 
                          self.outputFile, 
                          self.sqrts, 
                          self.muR4ggh, 
                          self.muF4ggh, 
                          self.muR4bbh, 
                          self.muF4bbh, 
                          self.BRcomputed, 
                          self.pdf)

    def setmA(self, M) :
        self.mA = M
    def settb(self, tb) :
        self.tb = tb
    def setsba(self, sba) :
        self.sba = sba
    def setmH(self, M) :
        self.mH = M
    def setmuF(self, muF4ggh, muF4bbh) :
        self.muF4ggh = muF4ggh
        self.muF4bbh = muF4bbh
    def setmuR(self, muR4ggh, muR4bbh) :
        self.muR4ggh = muR4ggh
        self.muR4bbh = muR4bbh
    def setmHc(self, M) :
        self.mhc = M
    def setm12(self, m12) :
        self.m12 = m12     
    def setm122(self, m122) :
        self.m122 = m122
    def setpdf(self,pdf) :
        self.pdf=pdf
    def setoutputFile(self,outputFile) :
        self.outputFile=outputFile

    def getXsecFromSusHi(self, sushiCardName = 'sushi.out') :
        
        def float_to_str(f):
            return str(f).replace('.', 'p')
        
        sushiCardName = 'mH_{}_mA_{}_tb_{}_cba_{}_mode_{}'.format(
                float_to_str(self.mH), float_to_str(self.mA), float_to_str(self.tb), float_to_str(self.sba), str(self.mode)) 
        
        sushiDefaultCardPath = "default_cards/default_sushi.in"
        sushiInputCardPath   = "Scan/" + str(self.pdf) +'/' + sushiCardName + ".in"
        sushiOutputCardPath  = "Scan/" + str(self.pdf) +'/' + sushiCardName + ".out"
        
        path_tosushi_output_cards="Scan/"+ str(self.pdf)
        if not os.path.exists(path_tosushi_output_cards):
            os.makedirs(path_tosushi_output_cards)
        
        # Replacements of variables into the input file
        replacements = {
            'MODE':str(self.mode),
            'TANBETA':str(self.tb),
            'M12':str(self.m12),
            'MSMH':str(self.mh),
            'MHEAVYH':str(self.mH),
            'MPSA':str(self.mA),
            'MCHARGEDH':str(self.mhc),
            'SINBA':str(self.sba),
            'MUR4ggh':str(self.muR4ggh),
            'MUF4ggh':str(self.muF4ggh),
            'MUR4bbh':str(self.muR4bbh),
            'MUF4bbh':str(self.muF4bbh),
            'TYPE':str(int(self.type)),
            'SQRTS':str(self.sqrts),
            'CUSTOMPDFNNLO':str(self.pdf)
            }

        #print replacements
        with open(sushiDefaultCardPath) as f1:
            with open(sushiInputCardPath,'w') as f2:
                for line in f1:
                    newline = line
                    for src, target in replacements.iteritems():
                        newline = newline.replace(src, target)
                    f2.write(newline)
     
        #running SusHi
        pwd = os.getcwd()
        
        #run_sushi = ["../../SusHi-1.7.0/bin/sushi", sushiCardName + ".in", sushiCardName + ".out"]
        #run_sushi = ["{}/SusHi-1.7.0/bin/sushi".format(pwd), sushiCardName + ".in", sushiCardName + ".out"]
        run_sushi  = "{}/SusHi-1.7.0/bin/sushi {}.in {}.out".format(pwd, sushiCardName, sushiCardName)
        print( "Running command:", run_sushi )

        sushi_out = 'Scan/sushilogs'
        if not os.path.isdir(sushi_out):
            os.makedirs(sushi_out)

        # overwrite these files 
        #outputExists = False 
        outputExists  = os.path.isfile(os.path.join(pwd, sushiOutputCardPath) )
        if not outputExists:
            logFile = os.path.join(pwd, sushi_out, 'sushi1.7.0_'+sushiCardName+'.log')
            print("Current working directory:", os.getcwd())
            with open(logFile, 'w+') as f:
                p = subprocess.Popen(run_sushi, stdout=f, stderr=f, shell=True, cwd=os.path.join(pwd, path_tosushi_output_cards))
                p.communicate()
        else:
            print('# SusHi was run already, looking for results in %s' % sushiOutputCardPath)
        
        os.chdir(pwd)
        integerror_gg = 0 # +/- integ. error: ggh XS in pb
        muRm_gg = 0 # - from muR variation: ggh XS in pb
        muRp_gg = 0 # + from muR variation: ggh XS in pb
        Xsec_gg = 0 # ggh XS in pb
        
        integerror_bb = 0 # +/- integ. error: ggh XS in pb
        muRm_bb = 0 # - from muR variation: ggh XS in pb
        muRp_bb = 0 # + from muR variation: ggh XS in pb
        Xsec_bb = 0 # bbh XS in pb
        
        ggh_lo= 0              # LO w/ NLO PDFs
        ggh_nlo= 0             # NLO
        subprocess_gg_nlo = 0  # NLO gg
        subprocess_qg_nlo = 0  # NLO qg
        subprocess_qq_nlo = 0  # NLO qq
                
        integerror_ggh_lo = 0  # +/-: LO w/ NLO PDFs
        integerror_ggh_nlo= 0  # +/-: NLO
        integerror_gg_nlo = 0  # +/-: NLO gg
        integerror_qg_nlo = 0  # +/-: NLO qg
        integerror_qq_nlo = 0  # +/-: NLO qq
        
        bbh_lo =  0              # LO
        bbh_nlo = 0              # NLO
        bbh_nnlo = 0             # NNLO
        integerror_bbh_lo = 0    # +/-: LO
        integerror_bbh_nlo = 0   # +/-: NLO
        integerror_bbh_nnlo = 0  # +/-: NNLO            
        
        mb_MSscheme_muR = 0.
        
        # extracting xsec from the output file
        cross_section = None
        if os.path.exists(os.path.join(pwd, sushiOutputCardPath)):
            with open(os.path.join(pwd, sushiOutputCardPath),'r') as f:
                Block = None
                for line in f:

                    if line.startswith('Block'):
                            Block = line.split()[1]

                    if '# m_b for bottom Yukawa' in line:
                        mb_MSscheme_muR = line.split()[1]
                    
                    if Block == 'SUSHIggh':
                        if 'ggh XS in pb' in line:
                            if ' 1 ' in line:
                                Xsec_gg = line.split()[1]
                            elif ' 101 ' in line:
                                integerror_gg = line.split()[1]
                            elif ' 102 ' in line:
                                muRm_gg = line.split()[1]
                            elif ' 103 ' in line:
                                muRp_gg = line.split()[1]

                    if Block == 'SUSHIbbh':
                        if 'bbh XS in pb' in line:
                            if ' 1 ' in line:
                                Xsec_bb = line.split()[1]
                            elif ' 101 ' in line:
                                integerror_bb = line.split()[1]

                    if Block == 'XSGGH':
                        if ' 1 ' in line:
                            ggh_lo = line.split()[1] ## LO w/ NLO PDFs
                        elif ' 101 ' in line:
                            integerror_ggh_lo = line.split()[1]
                        
                        elif ' 2 ' in line:
                            ggh_nlo = line.split()[1]
                        elif ' 102 ' in line:
                            integerror_ggh_nlo = line.split()[1]
                        
                        if ' NLO gg' in line:
                            if ' 21 ' in line:
                                subprocess_gg_nlo = line.split()[1]
                            elif ' 121 ' in line:
                                integerror_gg_nlo = line.split()[1]
                        
                        if ' NLO qg' in line:
                            if ' 22 ' in line:
                                subprocess_qg_nlo = line.split()[1]
                            elif ' 122 ' in line:
                                integerror_qg_nlo = line.split()[1]

                        if ' NLO qq' in line:
                            if ' 23 ' in line:
                                subprocess_qq_nlo = line.split()[1]
                            elif ' 123 ' in line:
                                integerror_qq_nlo = line.split()[1]
                        
                    if Block == 'XSBBH':
                        if ' 1 ' in line:
                            bbh_lo = line.split()[1]
                        elif ' 2 ' in line:
                            bbh_nlo = line.split()[1]
                        elif ' 3 ' in line:
                            bbh_nnlo = line.split()[1]
                        elif ' 101 ' in line:
                            integerror_bbh_lo = line.split()[1]
                        elif ' 102 ' in line:
                            integerror_bbh_nlo = line.split()[1]
                        elif ' 103 ' in line:
                            integerror_bbh_nnlo = line.split()[1]
                    
                    if Block =='HGGSUSY':
                        # we don't want to go there 
                        # this will overwrite the values we already get
                        break
            
            cross_section = {"split": 
                                    { 'gg_nlo'  : (float(subprocess_gg_nlo), float(integerror_gg_nlo) ), 
                                      'qg_nlo'  : (float(subprocess_qg_nlo), float(integerror_qg_nlo) ),
                                      'qq_nlo'  : (float(subprocess_qq_nlo), float(integerror_qq_nlo) ),
                                      'bbh_lo'  : (float(bbh_lo),   float(integerror_bbh_lo) ),
                                      'bbh_nlo' : (float(bbh_nlo),  float(integerror_bbh_nlo) ),
                                      'bbh_nnlo': (float(bbh_nnlo), float(integerror_bbh_nnlo) )
                                      },
                            "full": { 'ggh': (float(Xsec_gg), float(integerror_gg)), 
                                      'bbh': (float(Xsec_bb), float(integerror_bb)) }
                            }
        return cross_section 


    def computeBR(self):
        pwd = os.getcwd()
        print(pwd)
        command = ["./2HDMC-1.8.0/CalcPhys", 
                   str(self.mh), 
                   str(self.mH), 
                   str(self.mA), 
                   str(self.mhc), 
                   str(self.sba), 
                   "0", 
                   "0", 
                   str(self.m12_2), 
                   str(self.tb), 
                   str(self.type), 
                   self.outputFile]
        
        # overwrite these files 
        #outputExists = False
        outputExists = os.path.isfile(self.outputFile)
        logFile = self.outputFile.replace('.dat', '.log')
        if not outputExists:
            print (' '.join(command))
            with open(logFile, 'w+') as f:
                p = subprocess.Popen(command, stdout=f, stderr=f, cwd=pwd)
                p.communicate() # wait until process finishes
        else:
            print ('# 2HDMC was run already, looking for results in %s' % self.outputFile)
        
        if os.path.isfile(logFile):
            with open(os.path.join(logFile)) as f:
                for line in f:
                    if 'A  -> b  b' in line:
                        myLine=line.replace('A  -> b  b     ', '')
                        myLine2=myLine.split()[0]
                        self.wh3tobb = float(myLine2)
                    if 'H  -> b  b' in line:
                        myLine=line.replace('H  -> b  b     ', '')
                        myLine2=myLine.split()[0]
                        self.wh2tobb = float(myLine2)

        UnitLine = None
        PertLine = None
        StabLine = None
        with open(os.path.join(self.outputFile), 'r') as f:
            for line in f:
                if 'unitarity' in line:
                    UnitLine = line
                    self.unitarity = bool(int(line.split()[1]))
                if 'Perturbativity' in line:
                    PertLine = line
                    self.perturbativity = bool(int(line.split()[1]))
                if 'Stability' in line:
                    StabLine = line
                    self.stability = bool(int(line.split()[1]))
                if UnitLine and PertLine and StabLine:
                # stop parsing the file if we have all the informations...
                    break 
            
        modeh =0
        modeH =0
        modeA =0
        modeHc=0  
        with open(os.path.join(self.outputFile), 'r') as f:
            for line in f:
                if "DECAY  25" in line :
                    self.hwidth = float(line.split()[2])
                    modeh = 1
                    modeH = 0
                    modeA = 0
                elif "DECAY  35" in line :
                    self.Hwidth = float(line.split()[2])
                    modeh = 0
                    modeH = 1
                    modeA = 0
                elif "DECAY  36" in line :
                    self.Awidth = float(line.split()[2])
                    modeh = 0
                    modeA = 1
                    modeH = 0
                elif "DECAY  37" in line :
                    self.chargedHwidth = float(line.split()[2])
                    modeh  = 0
                    modeA  = 0
                    modeH  = 0
                    modeHc = 1

                elif "23    36" in line :
                    ZABRLine2 = line.replace("       ","")
                    ZABRLine3 = ZABRLine2.replace("     2      23    36","")
                    if modeH == 1 :
                        self.HtoZABR = float(ZABRLine3)
                    elif modeh == 1 :
                        self.htoZABR = float(ZABRLine3)
    
                elif "23    36" in line :
                    ZABRLine2 = line.replace("       ","")
                    ZABRLine3 = ZABRLine2.replace("     2      23    36","")
                    if modeH == 1 :
                        self.HtoZABR = float(ZABRLine3)
                    elif modeh == 1 :
                        self.htoZABR = float(ZABRLine3)
    
                elif "23    35" in line and modeA == 1 :
                    ZHBRLine2 = line.replace("       ","")
                    ZHBRLine3 = ZHBRLine2.replace("     2      23    35","")
                    self.AtoZHBR = float(ZHBRLine3)
            
                elif "36    36" in line :
                    AABRLine2 = line.replace("       ","")
                    AABRLine3 = AABRLine2.replace("     2      36    36","")
                    if modeH == 1 :
                        self.HtoAABR = float(AABRLine3)
                    elif modeh == 1 :
                        self.htoAABR = float(AABRLine3)
    
                elif "35    35" in line and modeH == 1 :
                  HHBRLine2 = line.replace("       ","")
                  HHBRLine3 = HHBRLine2.replace("     2      35    35","")
                  self.AtoHHBR = float(HHBRLine3)

                elif "3    -3" in line :
                  ssBRLine2 = line.replace("     2       3    -3","")
                  ssBRLine3 = ssBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtossBR = float(ssBRLine3)
                  elif modeA == 1 :
                    self.AtossBR = float(ssBRLine3)
                  elif modeh == 1 :
                    self.htossBR = float(ssBRLine3)

                elif "4    -4" in line :
                  ccBRLine2 = line.replace("     2       4    -4","")
                  ccBRLine3 = ccBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoccBR = float(ccBRLine3)
                  elif modeA == 1 :
                    self.AtoccBR = float(ccBRLine3)
                  elif modeh == 1 :
                    self.htoccBR = float(ccBRLine3)

                elif "11   -11" in line :
                  eeBRLine2 = line.replace("     2      11   -11","")
                  eeBRLine3 = eeBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoeeBR = float(eeBRLine3)
                  elif modeA == 1 :
                    self.AtoeeBR = float(eeBRLine3)
                  elif modeh == 1 :
                    self.htoeeBR = float(eeBRLine3)

                elif "13   -13" in line :
                  mumuBRLine2 = line.replace("     2      13   -13","")
                  mumuBRLine3 = mumuBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtomumuBR = float(mumuBRLine3)
                  elif modeA == 1 :
                    self.AtomumuBR = float(mumuBRLine3)
                  elif modeh == 1 :
                    self.htomumuBR = float(mumuBRLine3)

                elif "23    22" in line :
                  ZgaBRLine2 = line.replace("     2      23    22","")
                  ZgaBRLine3 = ZgaBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoZgaBR = float(ZgaBRLine3)
                  elif modeA == 1 :
                    self.AtoZgaBR = float(ZgaBRLine3)
                  elif modeh == 1 :
                    self.htoZgaBR = float(ZgaBRLine3)
                
                elif "5    -5" in line :
                  bbBRLine2 = line.replace("     2       5    -5","")
                  bbBRLine3 = bbBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtobbBR = float(bbBRLine3)
                  elif modeA == 1 :
                    self.AtobbBR = float(bbBRLine3)
                  elif modeh == 1 :
                    self.htobbBR = float(bbBRLine3)
            
                elif "23    25" in line :
                  ZhBRLine2 = line.replace("     2      23    25","")
                  ZhBRLine3 = ZhBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoZhBR = float(ZhBRLine3)
                  elif modeA == 1 :
                    self.AtoZhBR = float(ZhBRLine3)
            
                elif "6    -6" in line :
                  ttBRLine2 = line.replace("     2       6    -6","")
                  ttBRLine3 = ttBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtottBR = float(ttBRLine3)
                  elif modeA == 1 :
                    self.AtottBR = float(ttBRLine3)
                  elif modeh == 1 :
                    self.htottBR = float(ttBRLine3)
            
                elif "15   -15" in line :
                  tautauBRLine2 = line.replace("     2      15   -15","")
                  tautauBRLine3 = tautauBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtotautauBR = float(tautauBRLine3)
                  elif modeA == 1 :
                    self.AtotautauBR = float(tautauBRLine3)
                  elif modeh == 1 :
                    self.htotautauBR = float(tautauBRLine3)
            
                elif "25    25" in line :
                  hhBRLine2 = line.replace("     2      25    25","")
                  hhBRLine3 = hhBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtohhBR = float(hhBRLine3)
                  elif modeA == 1 :
                    self.AtohhBR = float(hhBRLine3)
            
                elif "23    23" in line :
                  ZZBRLine2 = line.replace("     2      23    23","")
                  ZZBRLine3 = ZZBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoZZBR = float(ZZBRLine3)
                  elif modeA == 1 :
                    self.AtoZZBR = float(ZZBRLine3)
                  elif modeh == 1 :
                    self.htoZZBR = float(ZZBRLine3)
            
                elif "24   -24" in line :
                  WWBRLine2 = line.replace("     2      24   -24","")
                  WWBRLine3 = WWBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoWWBR = float(WWBRLine3)
                  elif modeA == 1 :
                    self.AtoWWBR = float(WWBRLine3)
                  elif modeh == 1 :
                    self.htoWWBR = float(WWBRLine3)
            
                elif "22    22" in line :
                  ggBRLine2 = line.replace("     2      22    22","")
                  ggBRLine3 = ggBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoggBR = float(ggBRLine3)
                  elif modeA == 1 :
                    self.AtoggBR = float(ggBRLine3)
                  elif modeh == 1 :
                    self.htoggBR = float(ggBRLine3)
           
                elif "21    21" in line :
                  glugluBRLine2 = line.replace("     2      21    21","")
                  glugluBRLine3 = glugluBRLine2.replace("       ","")
                  if modeH == 1 :
                    self.HtoglugluBR = float(glugluBRLine3)
                  elif modeA == 1 :
                    self.AtoglugluBR = float(glugluBRLine3)
                  elif modeh == 1 :
                    self.htoglugluBR = float(glugluBRLine3)
                
                elif "lambda_1" in line:
                    self.lambda_1 = line.split()[1]
                elif "lambda_2" in line:
                    self.lambda_2 = line.split()[1]
                elif "lambda_3" in line:
                    self.lambda_3 = line.split()[1]
                elif "lambda_4" in line:
                    self.lambda_4 = line.split()[1]
                elif "lambda_5" in line:
                    self.lambda_5 = line.split()[1]
                elif "lambda_6" in line:
                    self.lambda_6 = line.split()[1]
                elif "lambda_7" in line:
                    self.lambda_7 = line.split()[1]
