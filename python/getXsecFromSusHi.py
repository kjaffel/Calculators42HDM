#!/usr/bin/python

import os
import sys


def getXsecFromSusHi( mode = 'H', sqrts = 13000, tb = 1, m12 = 0, mh = 125, mH = 350, mA = 350, 
        mhc = 350, sba = 0.99, muR4ggh = 0.5, muF4ggh = 0.5, muR4bbh = 1., muF4bbh = 0.25):

    sushiDefaultCardPath = "default_cards/default_sushi.in"
    sushiInputCardPath = "Scan/"+str(mH)+"_"+str(mA)+".in"
    sushiOutputCardPath = "Scan/"+str(mH)+"_"+str(mA)+".out"

    # Replacements of variables into the input file
    replacements = {'TANBETA':str(tb),
                    'M12':str(m12),
                    'MSMH':str(mh),
                    'MHEAVYH':str(mH), 
                    'MPSA':str(mA), 
                    'MCHARGEDH':str(mhc), 
                    'SINBA':str(sba), 
                    'MUR4ggh':str(muR4ggh), 
                    'MUF4ggh':str(muF4ggh),
                    'MUR4bbh':str(muR4bbh), 
                    'MUF4bbh':str(muF4bbh),
                    }
    sushiDefaultCard = open(sushiDefaultCardPath)
    sushiInputCard = open(sushiInputCardPath, 'w')
    for line in sushiDefaultCard:
        for src, target in replacements.iteritems():
            line = line.replace(src, target)
        sushiInputCard.write(line)
    sushiDefaultCard.close()
    sushiInputCard.close()

    #running SusHi
    run_sushi="./SusHi-1.7.0/bin/sushi "+sushiInputCardPath+" "+sushiOutputCardPath+" "
    os.system(run_sushi)

    Xsec = None
    # extracting xsec from the output file
    with open(sushiOutputCardPath,'r') as f:
        for line in f:
            if '# ggh XS in pb' not in f:
                continue
            print (line.split())
            Xsec = line.split()[1]
            break
    return Xsec 
    

