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
echo -e "\n"
echo -e "Enter option \c"
todooption
