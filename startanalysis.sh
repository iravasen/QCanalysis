#!/bin/bash

remove(){
  rm *.d
  rm *.pcm
  rm *.so
}

todo(){
  read answer
  case "$answer" in
    1) echo -e "\n=> Starting fake-hit rate analysis run by run"
       root -l -b -q AnalyzeStaveHitmaps.C++ 
       remove ;;
    2) echo -e "\n=> Starting comparison between runs"
       root -l -b -q CompareNoisyPixelsInRuns.C++
       remove ;;
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

todoinflp(){
  read answerflp
  case "$answerflp" in
    1) echo -e "\n=> Starting thresholds analysis run by run"
       root -l -b -q AnalyzeThrScanAvgThr.C++
       remove ;;
    2) echo -e "\n=> Starting dead pixel comparison between runs"
       root -l -b -q CompareDeadPixelsInRuns.C++
       remove ;;
    *) echo -e "Invalid option \n"
       echo -e "Retype an option \c"
       todoinflp ;;
  esac
}

analysismenuonflp(){
  echo -e "\n=> Choose the analysis you want to perform \n"
  echo "[Threshold scan analyses]"
  echo -e "\t 1. Average stave thresholds run by run"
  echo -e "\t 2. Compare dead pixels between runs"
  echo -e "\n"
  echo -e "Enter option \c"
  todoinflp
}

updategitrepo(){
  echo -e "\n=> Updating the git repository (your modification will be kept!)"
  git stash #in case there is something modified by the user
  git pull --rebase
  git stash pop #apply last modifications from the user
}

doanalysisinflp(){
  echo -e "\n=> Which ITS Layer do you want to analyse [0,1,2,3,4,5,6]? \c"
  read layernum
  case "$layernum" in
    0) cd $1
       updategitrepo
       echo -e "\n=> Starting analysis on flp"
       echo -e "Load QualityControl environment on flp"
       eval $(alienv load QualityControl/latest)
       echo -e "\n=> Preparation of the sample (may take several minutes)\n"
       cd analysismacros
       find /data/L0_shifts/ -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
       python readthrdata.py "IB" 0
       rm datatoanalyse.txt
       analysismenuonflp
       ;;
    1) cd $1
       updategitrepo
       echo -e "\n=> Starting analysis on flp"
       echo -e "Load QualityControl environment on flp"
       eval $(alienv load QualityControl/latest)
       echo -e "\n=> Preparation of the sample (may take several minutes)\n"
       cd analysismacros
       echo -e "In which folder the data are saved? \c"
       read foldername
       find $foldername -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
       python readthrdata.py "IB" 1
       rm datatoanalyse.txt
       analysismenuonflp
       ;;
     
     2) cd $1
       updategitrepo
       echo -e "\n=> Starting analysis on flp"
       echo -e "Load QualityControl environment on flp"
       eval $(alienv load QualityControl/latest)
       echo -e "\n=> Preparation of the sample (may take several minutes)\n"
       cd analysismacros
       echo -e "In which folder the data are saved? \c"
       read foldername
       find $foldername -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
       python readthrdata.py "IB" 2
       rm datatoanalyse.txt
       analysismenuonflp
       ;;


    *) echo -e "\nLayer not yet available"
       doanalysisinflp ;;
  esac
}

chooseflp(){
  echo -e "\nIn which flp your data are [type the number only]? \c"
  read flpnum
  case "$flpnum" in
    1) path="/home/its/QCNew/QCanalysis" ;;
    7) path="/home/its/QC/QCanalysis" ;;
    *) echo "Invalid number (for the moment)"
       chooseflp ;;
   esac
}


connect(){
  echo -e "\nInsert you CERN username: \c"
  read usercern
  echo -e "... Connecting to lxplus and to $1"
  ssh -o "ProxyCommand ssh $usercern@lxplus.cern.ch -q nc %h 22" its@$1.cern.ch "$(typeset -f); doanalysisinflp $path"
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
       hname=$(hostname)
       chooseflp
       if [ $hname == "flpits$flpnum" ]; then
	doanalysisinflp $path
       else
	connect "flpits$flpnum"
       fi; ;;

       #if [ $username != "its" ] && [ $hname != "flpits1" ] && [ $hname != "flpits7" ]; then
	#echo -e "\nTo which flp you want to connect [type the number only]? \c"
	#read flpnum
	#connect "flpits$flpnum"
       #else
	#doanalysisinflp
       #fi; ;;
    *) echo -e "Invalid option \n"
       echo -e "Retype an option \c"
       todooption ;;
  esac
}

# MAIN TASK
echo "==== Updating git repository (your modification are kept!) ===="
updategitrepo

echo "\n==== Loading environment modules ===="
#export ALIBUILD_WORK_DIR="/home/alidock/.sw"
eval $(alienv load QualityControl/latest 2> /dev/null)

echo -e "\n==== What to do ====\n"
echo -e "=> Chose an option: \n\n"
echo -e "\t 1. Download and analyse data"
echo -e "\t 2. Analyse data only"
echo -e "\t 3. Analyse data on flp (for threshold scan)"
echo -e "\n"
echo -e "Enter option \c"
todooption
