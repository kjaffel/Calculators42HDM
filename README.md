
# 2HDMC - Two-Higgs-Doublet Model Caclulators -
2HDMC is a general-purpose calculator for the two-Higgs doublet model. It allows parametrization of the Higgs potential in many different ways, convenient specification of generic Yukawa sectors, the evaluation of decay widths (including higher-order QCD corrections), theoretical constraints and much more.

Version used: ``2HDMC1.8.0`` and ``Sushi1.7.0``.
## Install in a CMSSW release:
```bash
    # setup your env : for cp3 ingrid-ui1 cluster 
    module load gcc/gcc-7.3.0-sl7_amd64 lhapdf/6.1.6-sl7_gcc73
    
    # Install a CMSSW release .eg. CMSSW_10_2_22
    cmsrel CMSSW_10_2_22
    cd CMSSW_10_2_22/src
    cmsenv
    git cms-init

    # Get and execute the install script
    wget https://raw.githubusercontent.com/kjaffel/Calculators42HDM/master/sushi_2hdmc_cmssw.sh 
    source sushi_2hdmc_cmssw.sh
```
## Install With Conda (from [@pdavid](https://github.com/pieterdavid)):
- [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html)
- [CONDA-FORGE](https://conda-forge.org/#page-top)

```bash
    # setup your env
    conda create -n SusHi170 gsl gfortran_linux-64 gxx_linux-64 gcc_linux-64
    conda activate SusHi170

    wget https://github.com/kjaffel/Calculators42HDM/master/sushi_2hdmc_conda.sh
    source sushi_2hdmc_conda.sh
```

## Run the test:
   
    python example/test.py
 
## Trouble-Shooting: 

```bash
lib/libsushi2HDMC.a(runthdmc.o): In function `runthdmc_':
runthdmc.f:(.text+0x13d): undefined reference to `thdmc_set_param_'
collect2: error: ld returned 1 exit status
make[1]: *** [Makefile:197: bin/sushi.2HDMC] Error 1
```
In ``2HDMC-1.8.0/Makefile`` replace these lines:
```bash
  -CFLAGS= -std=c++11 -Wall $(DEBUG) $(OPT)
  +CFLAGS= -std=c++11 -Wall -fPIE $(DEBUG) $(OPT)
  -SOURCES=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp
  +SOURCES=THDM.cpp SM.cpp DecayTable.cpp Constraints.cpp Util.cpp runTHDM.cpp
```    

## Keep an eye on the latest versions: 
- [2HDMC](https://2hdmc.hepforge.org)
- [SUSHI](https://sushi.hepforge.org/)
- [LHAPDF](https://lhapdf.hepforge.org/index.html)
- It can also be useful to download HiggsBounds/HiggsSignals (Note: 2HDMC work without HB/HS), since 2HDMC can interface the LEP, Tevatron and LHC constraints implemented in these codes. HiggsBounds/HiggsSignals are available: 
- [HiggsBounds](http://higgsbounds.hepforge.org)
- [HiggsSignals](http://higgsbounds.hepforge.org)
## Useful Links:  
- [2HDME: Two-Higgs-Doublet Model Evolver](https://arxiv.org/pdf/1811.08215.pdf) 
- [Flavour Les Houches Accord: Interfacing Flavour related Codes](https://arxiv.org/pdf/1008.0762.pdf)
