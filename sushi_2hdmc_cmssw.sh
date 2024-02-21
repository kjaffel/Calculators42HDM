#!/bin/bash

DIR='cp3_llbb/Calculators42HDM'
set -x
# inside a CMSSW release
if [[ ! -d "$DIR" ]]; then
    echo "clone ..."
    git clone -o upstream git@github.com:kjaffel/Calculators42HDM.git $DIR
fi
scram b
eval `scramv1 runtime -sh` # This is equivalent to cmsenv

pwd
pushd $DIR

#echo "--> Building LHAPDF"
#wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.3.0.tar.gz -O LHAPDF-6.3.0.tar.gz
### ^ or use a web browser to download, which will get the filename correct
#tar xf LHAPDF-6.3.0.tar.gz
#pushd LHAPDF-6.3.0
#./configure #--prefix=/path/for/installation
#make
#make install
#popd
#
#echo "--> Building Higgs Bounds"
#wget https://higgsbounds.hepforge.org/downloads?f=HiggsBounds-5.3.2beta.tar.gz -O HiggsBounds-5.3.2beta.tar.gz
#tar -zxvf HiggsBounds-5.3.2beta.tar.gz
#pushd HiggsBounds-5.3.2beta
#./configure
#make
#popd
#
#echo "--> Building Higgs Signals"
#wget https://higgsbounds.hepforge.org/downloads?f=HiggsSignals-2.2.3beta.tar.gz -O HiggsSignals-2.2.3beta.tar.gz
#tar -zxvf HiggsSignals-2.2.3beta.tar.gz
#pushd HiggsSignals-2.2.3beta
#./configure
#make
#popd
#
echo "--> Building 2HDMC"
wget http://www.hepforge.org/archive/2hdmc/2HDMC-1.8.0.tar.gz
tar -zxvf 2HDMC-1.8.0.tar.gz
pushd 2HDMC-1.8.0
make
popd

#echo "--> Building FeynHiggs"
#echo " please download the tar.gz file from http://www.feynhiggs.de/ and then continue !"
#tar -zxvf FeynHiggs-2.18.0.tar.gz
#pushd FeynHiggs-2.18.0
#./configure
#make
#make install
#popd
#
echo "--> Building SUSHI"
echo " NOTE: Sushi should be the last to configure: in order to link the previous packages !"
wget http://www.hepforge.org/archive/sushi/SusHi-1.7.0.tar.gz
tar -zxvf SusHi-1.7.0.tar.gz
pushd SusHi-1.7.0
./configure
## Edit the makefile, to link the proper 2HDMC version 
#sed -i -e 's|^HiggsBounds_DIR =.*$|HiggsBounds_DIR = ../HiggsBounds-5.3.2beta|' Makefile
#sed -i -e 's|^HiggsSignals_DIR =.*$|HiggsSignals_DIR = ../HiggsSignals-2.2.3beta|' Makefile
sed -i -e 's|^2HDMCPATH =.*$|2HDMCPATH = ../2HDMC-1.8.0|' Makefile
sed -i -e 's|^2HDMCVERSION =.*$|2HDMCVERSION = 1.8.0|' Makefile
#sed -i -e 's|^FHPATH =.*$|FHPATH = ../FeynHiggs-2.18.0|' Makefile
#sed -i -e 's|^FHPATH =.*$|FHPATH = 2.18.0|' Makefile
make predef=2HDMC
#make predef=FH
popd

set +x
