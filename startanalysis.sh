#!/bin/bash

remove(){
  rm *.d
  rm *.pcm
  rm *.so
}

todo(){
  read answer
  case "$answer" in
    1) echo -e "\n\e[32m=> Starting fake-hit rate analysis run by run\e[39m"
       root -l -b -q AnalyzeLayerOccupancy.C++
       remove ;;
    2) echo -e "\n\e[32m=> Starting comparison between runs\e[39m"
       root -l -b -q CompareNoisyPixelsInRuns.C++
       remove ;;
    3) echo -e "\n\e[32m=> Starting fake-hit rate analysis masking hot pixels\e[39m"
       root -l -b -q MaskNoisyPixelsInRuns.C++
       remove ;;
    4) echo -e "\n\e[32m=> Starting fake-hit rate correlation analysis\e[39m"
       root -l -b -q CompareLayerOccupancy.C++
       remove ;;
    5) echo -e "\n\e[32m=> Starting analysis of the errors for all runs\e[39m"
       root -l -b -q AnalyzeErrorsFHR.C++
       remove ;;
    6) echo -e "\n\e[32m=> Starting threshold analysis run by run\e[39m"
       root -l -b -q AnalyzeLayerThresholds.C++
       remove ;;
    7) echo -e "\n\e[32m=> Starting dead pixel analysis run by run\e[39m"
       root -l -b -q AnalyzeLayerDeadPixels.C++
       remove ;;
    8) echo -e "\n\e[32m=> Starting threshold correlation analysis\e[39m"
       root -l -b -q CompareLayerThresholds.C++
       remove ;;
    9) echo -e "\n\e[32m=> Starting dead-pixels correlation analysis\e[39m"
       root -l -b -q CompareLayerDeadPixels.C++
       remove ;;
    10) echo -e "\n\e[32m=> Starting comparison of dead-pixels between runs\e[39m"
       root -l -b -q CompareDeadPixelsInRuns.C++
       remove ;;
    11) echo -e "\n\e[32m=> Starting preparation of dead-pixel map\e[39m"
       root -l -b -q MakeDeadPixelMap.C++
       remove ;;
    12) echo -e "\n\e[32m=> Starting analysis of the lane status flags for all runs\e[39m"
       root -l -b -q AnalyzeLaneStatusFlag.C++
       remove ;;
    13) echo -e "\n\e[32m=> Starting dump of lanes into error, fault, warning\e[39m"
       root -l -b -q DumpLaneStatusFlag.C++
       remove ;;
    14) echo -e "\n\e[32m=> Starting analysis of the trigger&flags for all runs\e[39m"
       root -l -b -q AnalyzeTrgFlg.C++
       remove ;;
    15) echo -e "\n\e[32m=> Starting average cluster size analysis run by run\e[39m"
       root -l -b -q AnalyzeLayerAverageClusterSize.C++
       remove ;;
    16) echo -e "\n\e[32m=> Starting cluster occupation analysis run by run\e[39m"
       root -l -b -q AnalyzeLayerClusterOccupancy.C++
       remove ;;
    17) echo -e "\n\e[32m=> Starting full analysis of track and vertex parameters\e39m"
       root -l -b -q Track_trending_postprocessing.C++
       remove ;;
    *) echo -e "Invalid option \n"
       echo -e "Retype an option \c"
       todo ;;
  esac
}

analysismenu(){
  echo -e "\n\e[32m=> Choose the analysis you want to perform \n"
  echo "[Analyses on Fake-Hit Rate]"
  echo -e "\t 1.  Fake-hit rate run by run"
  echo -e "\t 2.  Compare number of noisy pixels between runs"
  echo -e "\t 3.  Fake-hit rate study with hot pixels masking"
  echo -e "\t 4.  Fake-hit rate correlation analysis (reference run to be chosen)"
  echo -e "\t 5.  Error analysis for all runs"
  echo -e "\n"
  echo "[Analyses on Threshold runs]"
  echo -e "\t 6.  Average threshold run by run"
  echo -e "\t 7.  Total dead pixels run by run"
  echo -e "\t 8.  Threshold correlation analysis (reference run to be chosen)"
  echo -e "\t 9. Dead-pixels correlation analysis (reference run to be chosen)"
  echo -e "\t 10. Compare number of dead-pixels between runs"
  echo -e "\t 11. Make dead-pixels map for layer(s)"
  echo -e "\n"
  echo "[Analyses on FEE]"
  echo -e "\t 12. FEE Post Processing Offline: Lane Status Flags (ERROR,FAULT,WARNING)"
  echo -e "\t 13. Dump in txt file of lanes into ERROR, FAULT, WARNING"
  echo -e "\t 14. Trigger Flags analysis for all runs"
  echo -e "\n"
  echo "[Analyses on Clusters]"
  echo -e "\t 15.  Average cluster size run by run"
  echo -e "\t 16.  Cluster occupation run by run"
  echo -e "\n"
  echo "[Analyses on Tracks]"
  echo -e "\t 17.  Full trending of track and vertex parameters\e[39m"
  echo -e "\n"
  echo -e "Enter option \c"
  cd analysismacros
  todo
  cd ..
}

directanalysisoption(){
  echo -e "\n\e[32m=> Choose what to do for the analysis\n"
  echo -e "\t 1. Preparation of data sample and analysis"
  echo -e "\t 2. Direct analysis (I have already a data sample)\e[39m"
  echo -e "\n"
  echo -e "Enter option \c"
  read directoptan
}

todoinflp(){
  read answerflp
  case "$answerflp" in
    1) directanalysisoption
       case "$directoptan" in
         1) echo -e "\n\e[32m=> Preparation of the sample (may take several minutes)\n\e[39m"
	    find $foldername -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
            python readthrdata.py $1 $2
            rm datatoanalyse.txt
            echo -e "\n\e[32m=> Starting thresholds analysis run by run\e[39m"
            root -l -b -q AnalyzeThrScanAvgThr.C++
            remove
            ;;
         2) echo -e "\n\e[32m=> Starting thresholds analysis run by run\e[39m"
            root -l -b -q AnalyzeThrScanAvgThr.C++
            remove
            ;;
       esac
       ;;
    2) directanalysisoption
       case "$directoptan" in
         1) echo -e "\n\e[32m=> Preparation of the sample (may take several minutes)\n\e[39m"
	    find $foldername -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
            python readthrdata.py $1 $2
            rm datatoanalyse.txt
            echo -e "\n\e[32m=> Starting dead pixel comparison between runs\e[39m"
            root -l -b -q CompareDeadPixelsInRuns.C++
            remove
            ;;
         2) echo -e "\n\e[32m=> Starting dead pixel comparison between runs\e[39m"
            root -l -b -q CompareDeadPixelsInRuns.C++
            remove
            ;;
       esac
       ;;
    3) directanalysisoption
       case "$directoptan" in
         1) echo -e "\n\e[32m=> Preparation of the sample (may take several minutes)\n\e[39m"
	    find $foldername -name "thresholds.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
            python readthrdata.py $1 $2
            echo -e "\n\e[32m=> Starting thresholds analysis run by run\e[39m"
            root -l -b -q AnalyzeThrScanAvgThr.C++
            remove
            echo -e "\n\e[32m=> Starting dead pixel comparison between runs\e[39m"
            root -l -b -q CompareDeadPixelsInRuns.C++
            remove
            rm datatoanalyse.txt
            ;;
         2) echo -e "\n\e[32m=> Starting thresholds analysis run by run\e[39m"
            root -l -b -q AnalyzeThrScanAvgThr.C++
            remove
            echo -e "\n\e[32m=> Starting dead pixel comparison between runs\e[39m"
            root -l -b -q CompareDeadPixelsInRuns.C++
            remove
            ;;
        esac
       ;;
    4)  directanalysisoption
        case "$directoptan" in
          1)  echo -e "\n\e[32m=> Preparation of the sample (may take several minutes)\n\e[39m"
	      find $foldername -name "hitmap.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
              python readfhitdata.py $1 $2
              rm datatoanalyse.txt
              echo -e "\n\e[32m=> Starting fake-hit rate analysis run by run\e[39m"
              root -l -b -q AnalyzeStaveHitmaps_flp.C++
              remove
              ;;
          2)  echo -e "\n\e[32m=> Starting fake-hit rate analysis run by run\e[39m"
              root -l -b -q AnalyzeStaveHitmaps_flp.C++
              remove
              ;;
        esac
       ;;
    5)  directanalysisoption
        case "$directoptan" in
          1)  echo -e "\n\e[32m=> Preparation of the sample (may take several minutes)\n\e[39m"
	      find $foldername -name "hitmap.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
              python readfhitdata.py $1 $2
              rm datatoanalyse.txt
              echo -e "\n\e[32m=> Starting noisy pixels comparison between runs\e[39m"
              root -l -b -q CompareNoisyPixelsInRuns_flp.C++
              remove
              ;;
          2)  echo -e "\n\e[32m=> Starting noisy pixels comparison between runs\e[39m"
              root -l -b -q CompareNoisyPixelsInRuns_flp.C++
              remove
              ;;
        esac
       ;;
    6) directanalysisoption
       case "$directoptan" in
         1) echo -e "\n\e[32m=> Preparation of the sample (may take several minutes)\n\e[39m"
	    find $foldername -name "hitmap.npy.gz" -print0 | sort -z | xargs -r0 | tr " " "\n" > datatoanalyse.txt
            python readfhitdata.py $1 $2
            echo -e "\n\e[32m=> Starting fake-hit rate analysis run by run\e[39m"
            root -l -b -q AnalyzeStaveHitmaps_flp.C++
            remove
            echo -e "\n\e[32m=> Starting noisy pixels comparison between runs\e[39m"
            root -l -b -q CompareNoisyPixelsInRuns_flp.C++
            remove
            rm datatoanalyse.txt
            ;;
         2) echo -e "\n\e[32m=> Starting fake-hit rate analysis run by run\e[39m"
            root -l -b -q AnalyzeStaveHitmaps_flp.C++
            remove
            echo -e "\n\e[32m=> Starting noisy pixels comparison between runs\e[39m"
            root -l -b -q CompareNoisyPixelsInRuns_flp.C++
            remove
            ;;
        esac
       ;;

    *) echo -e "Invalid option \n"
       echo -e "Retype an option \c"
       todoinflp ;;
  esac
}

analysismenuonflp(){
  echo -e "\n\e[32m=> Choose the analysis you want to perform \n"
  echo "[Threshold scan analyses]"
  echo -e "\t 1. Average stave thresholds run by run"
  echo -e "\t 2. Compare dead pixels between runs"
  echo -e "\t 3. Option 1 and 2 together"
  echo "[Analyses on noisy pixels]"
  echo -e "\t 4. Fake-hit rate run by run"
  echo -e "\t 5. Compare noisy pixels between runs"
  echo -e "\t 6. Option 4 and 5 together\e[39m"
  echo -e "\n"
  echo -e "Enter option \c"
  todoinflp $1 $2 # $1 = "IB or OB" and $2 = layer number
}

updategitrepo(){
  echo -e "\n\e[32m=> Updating the git repository (your modification will be kept!)\e[39m"
  git stash #in case there is something modified by the user
  git pull --rebase
  git stash pop #apply last modifications from the user
}

startflp(){
  echo -e "\n\e[32m=> Starting analysis on flp\e[39m"
  echo -e "\e[32mLoading QualityControl environment on flp\e[39m"
  eval $(alienv load QualityControl/latest)
  cd analysismacros
  echo -e "\e[32mIn which folder the data are saved?\e[39m \c"
  read foldername
  analysismenuonflp $1 $2 # $1 = "IB or OB" and $2 = layer number
}


doanalysisinflp(){
  echo -e "\n\e[32m=> Which ITS Layer do you want to analyse [0,1,2,3,4,5,6]?\e[39m \c"
  read layernum
  case "$layernum" in
    0) cd $1
       updategitrepo
       startflp "IB" 0
       ;;
    1) cd $1
       updategitrepo
       startflp "IB" 1
       ;;

     2) cd $1
       updategitrepo
       startflp "IB" 2
       ;;

    *) echo -e "\nLayer not yet available"
       doanalysisinflp ;;
  esac
}

chooseflp(){
  echo -e "\n\e[32mIn which flp your data are [type the number only]?\e[39m \c"
  read flpnum
  case "$flpnum" in
    1) path="/home/its/QCNew/QCanalysis" ;;
    7) path="/home/its/QC/QCanalysis" ;;
    *) echo "Invalid number (for the moment)"
       chooseflp ;;
   esac
}


connect(){
  echo -e "\n\e[32mInsert you CERN username:\e[39m \c"
  read usercern
  echo -e "\n\e[32mConnecting to lxplus and to $1\e[39m"
  ssh -o "ProxyCommand ssh $usercern@lxplus.cern.ch -q nc %h 22" its@$1.cern.ch "$(typeset -f); doanalysisinflp $path"
}

todooption(){
  read answerfirst
  case "$answerfirst" in
    1) echo -e "\n\e[32m=> Choose what to do with the analysis:\n"
       echo -e "\t 1. Download data from CCDB and perform the analysis"
       echo -e "\t 2. Analyse data only (I have already a data sample)\e[39m"
       echo -e "\n"
       echo -e "Enter option: \c"
       read optccdb
       case "$optccdb" in
         1) echo -e "\n\e[32m=> Compiling the software for the database\e[39m"
            make
            echo -e "\n\e[32m=> Choose what to download\n"
            echo -e "1. Fake-hit rate data"
            echo -e "2. Threshold data"
            echo -e "3. Track data"
	    echo -e "4. FEE data"
	    echo -e "5. Cluster data\e[39m"
            echo -e "Enter option \c"
            read downloadoption
            case "$downloadoption" in
              1) echo -e "\n\e[32m=> Downloading data for fake-hit rate to analyse\e[39m"
                 ./getObject expert 1
                 echo -e "\n"
                 analysismenu ;;
              2) echo -e "\n\e[32m=> Downloading data for thresholds to analyse\e[39m"
                 ./getObject expert 2
                 echo -e "\n"
                 analysismenu ;;
              3) echo -e "\n\e[32m=> Downloading data for tracks to analyse\e[39m"
                 ./getObject expert 3
                 echo -e "\n"
                 analysismenu ;;
	      4) echo -e "\n\e[32m=> Downloading data for FEE to analyse\e[39m"
		 ./getObject expert 4
		 echo -e "\n"
                 analysismenu ;;
              5) echo -e "\n\e[32m=> Downloading data for clusters to analyse\e[39m"
      		 ./getObject expert 5
      		 echo -e "\n"
                 analysismenu;;
              *) echo -e "\n\e[32m=> Option not valid. Rerun the script. \e[39m" ;;
            esac
            ;;

         2) analysismenu ;;
       esac
       ;;
    2) username=$(whoami)
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

doexpert(){
  echo -e "\n\e[32m===== What to do ====\n"
  echo -e "=> Chose an option: \n\n"
  echo -e "\t 1. Download and analyse data from CCDB"
  echo -e "\t 2. Analyse data on FLP\e[39m"
  echo -e "\n"
  echo -e "\e[39mEnter option \c"
  todooption
}

doshifters(){
     echo -e "\n\e[32m=> Compiling the software for the database\e[39m"
     make
     echo -e "\n\e[32m=> Choose what to download\n"
     echo -e "1. Fake-hit rate runs"
     echo -e "2. Threshold runs\e[39m"
     echo -e "Enter option \c"
     read downloadoption
     case "$downloadoption" in
      1) echo -e "\n\e[32m=> Downloading files to analyse\e[39m"
         ./getObject shifter 1
         echo -e "\n"
         echo -e "\n\e[32m=> Starting analysis run by run and preparation of the 24h Report for fake-hit rate runs \e[39m"
         cd analysismacros
         root -l -b -q ReportFHR_shifters.C++
         cd ..
         remove ;;
      2) echo -e "\n\e[32m=> Downloading files to analyse\e[39m"
         ./getObject shifter 2
         echo -e "\n"
         echo -e "\n\e[32m=> Starting analysis run by run and preparation of the 24h Report for Threshold runs \e[39m"
         cd analysismacros
         root -l -b -q ReportTHR_shifters.C++
         cd ..
         remove ;;
      *) echo -e "\n\e[32m=> Option not valid. Rerun the script. \e[39m" ;;
    esac
}

# MAIN TASK
echo -e "\n\e[32m============================================\n"
echo -e "========== Welcome to QCanalysis ===========\e[39m"
echo -e "\e[31m== Any issue? Write to Ivan Ravasenga: ivan.ravasenga@cern.ch ==\e[39m"
echo -e "\n\e[32m============================================\e[39m"

updategitrepo

echo -e "\n\e[32m==== Loading environment modules ===="
eval $(alienv load QualityControl/latest 2> /dev/null)

#Choose mode: shifter / expert mode
echo -e "\n==== Choose mode ====\n"
echo -e "1. Shifter mode"
echo -e "2. Expert mode"
echo -e "\n"
echo -e "\e[39mEnter option \c"
read answeroption
case "$answeroption" in
  1) doshifters ;;
  2) doexpert ;;
esac
