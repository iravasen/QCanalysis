#!/bin/bash

todo(){
  read answer
  case "$answer" in
    1) echo -e "\n=> Starting fake-hit rate analysis run by run"
       root -l -q AnalyzeStaveHitmaps.C+ ;;
    2) echo -e "\n=> Starting comparison between runs"
       root -l -q CompareNoisyPixelsInRuns.C+ ;;
    *) echo -e "Invalid option \n"
       echo -e "Retype an option \c"
       todo ;;
  esac
}

analysismenu(){
  echo -e "\n=> Choose the analysis you want to perform \n"
  echo "[Analyses on noisy pixels]"
  echo -e "\t 1. Fake-hit rate run by run"
  echo -e "\t 2. Compare noisy pixels between runs"
  echo -e "\n"
  echo -e "Enter option \c"
  cd analysismacros
  todo
  cd ..
}

doanalysisinflp(){
  echo -e "\n=> Which ITS Layer do you want to analyse [0,1,2,3,4,5,6]? \c"
  read layernum
  case "$layernum" in
    0) cd /home/its/QCNew/QCanalysis
       echo -e "\n=> Starting analysis on flp"
       echo -e "Load QualityControl environment on flp"
       eval $(alienv load QualityControl/latest)
       echo -e "\nPreparation of the sample (may take several minutes)"
       cd analysismacros
       find /data/L0_shifts/ -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
       python readthrdata.py
       rm datatoanalyse.txt
       ;;
    *) echo -e "\nLayer not yet available"
       doanalysisinflp ;;
  esac
}

todooption(){
  read answerfirst
  case "$answerfirst" in
    1) echo -e "\n=> Compiling the software for the database"
       make
       echo -e "\n=> Downloading files to analyse"
       ./getObject
       echo -e "\n"
       analysismenu ;;
    2) analysismenu ;;
    3) username=$(whoami)
       if [ "$username"!="its" ]; then
	echo -e "\nInsert you CERN username: \c"
	read usercern
	echo -e "... Connecting to lxplus and to flpits1"
	ssh -o "ProxyCommand ssh $usercern@lxplus.cern.ch -q nc %h 22" its@flpits1.cern.ch "$(typeset -f doanalysisinflp); doanalysisinflp"
       else
	doanalysisinflp
       fi; ;;
    *) echo -e "Invalid option \n"
       echo -e "Retype an option \c"
       todooption ;;
  esac
}

# MAIN TASK

echo "==== Loading environment modules ===="
export ALIBUILD_WORK_DIR="/home/alidock/.sw"
eval $(alienv load QualityControl/latest 2> /dev/null)

echo -e "\n==== What to do ====\n"
echo -e "=> Chose an option: \n\n"
echo -e "\t 1. Download and analyse data"
echo -e "\t 2. Analyse data only"
echo -e "\t 3. Analyse data on flp (for threshold scan)"
echo -e "\n"
echo -e "Enter option \c"
todooption
