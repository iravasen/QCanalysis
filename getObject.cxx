#include <string>
#include <iostream>
#include <time.h>
#include <algorithm>
#include "QualityControl/DatabaseFactory.h"
#include "QualityControl/CcdbDatabase.h"
#include "QualityControl/MonitorObject.h"
#include "CCDB/CcdbApi.h"
#include <TBufferJSON.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TFile.h>
#include <TTree.h>

using namespace std;
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;

//functions to download data
string GetCorrectTS(string selrun, vector<string> runs, vector<string> timestamps);
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname);
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string timestamp);
vector<string> GetGoodRunList(o2::ccdb::CcdbApi ccdbApi, string run1, string run2, string runtype);
bool RunShifter(auto *ccdb, string myname, int opt);
bool RunExpert(auto *ccdb, string myname, int opt);
void DownloadTimestamps(auto *ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, long int ts_start, long int ts_end, int lnum);
void DownloadRuns(auto *ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, int lnum);
bool GetListOfHisto(auto* ccdb, string myname, string taskname, string objname, vector<long int> timestamps, vector<long int> timestamps2, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers, vector<int>runnumbers2);
bool Download(int choice, auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, long int ts_start, long int ts_end, int lnum);
string GetOptName(int opt);
string GetListName(int opt, int ilist);

const int nStavesInLay[7] = {12, 16, 20, 24, 30, 42, 48};
TFile *outputfile;

//to which CCDB we have to connect
// For P2 operations put: alio2-cr1-flp187.cern.ch:8083
string ccdbport = "ccdb-test.cern.ch:8080";


int main(int argc, char **argv)
{

  string myname = argv[0];

  std::unique_ptr<DatabaseInterface> mydb = DatabaseFactory::create("CCDB");
  if(mydb == nullptr) {
    cerr << myname << ": failed to create DatabaseInterface" << endl;
    return -1;
  }

  auto* ccdb = dynamic_cast<CcdbDatabase*>(mydb.get());
  if(ccdb == nullptr) {
    cerr << myname << ": ccdb pointer is null" << endl;
    return -1;
  }

  ccdb->connect(ccdbport.c_str(), "", "", "");

  if(strcmp(argv[1],"expert")==0)
    RunExpert(ccdb, myname, atoi(argv[2]));
  else if (strcmp(argv[1],"shifter")==0)
    RunShifter(ccdb, myname, atoi(argv[2]));


  ccdb->disconnect();

  return 1;
}//end main

//
//Get Last (most recent) Run
//
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname){

  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  stringstream ss(objectlist);
  string word;
  array<string,2> runts;

  while(ss>>word){
    if(word=="Created:"){// take the one related to file creation
      ss>>word;
      runts[0] = word;
    }
    if(word=="Run"){
      ss>>word;
      ss>>word;
      runts[1] = word;
      break;
    }
  }

  return runts;
}

//
// Get Run 24h depending on the actual time stamp
//
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string timestamp){

  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  stringstream ss(objectlist);
  string word;
  array<string,2> runts;
  bool islast = false;

  long int stamp_int_actual = stol(timestamp);
  long int stamp_int_24hago = stamp_int_actual - 86400000; //remove number of milliseconds in 1 day

  while(ss>>word){
    if(word=="Created:"){// take the one related to file creation
      ss>>word;
      if(stol(word)>=stamp_int_24hago){
        runts[0] = word;
      }
      else islast=true;
    }
    if(word=="Run"){
      ss>>word;
      ss>>word;
      if(islast) break;
      if(stol(runts[0])>=stamp_int_24hago) runts[1] = word;
    }
  }

  return runts;
}

//
// Expert mode
//
bool RunShifter(auto *ccdb, string myname, int opt){

  //Choose the layer number
  int layernum;
  cout<<endl;
  cout<<endl;
  cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
  cin>>layernum;

  //add error and trigger-flags plot to data
  bool adderrordata = true;
  if(opt==2) adderrordata = false;

  //taskname
  string taskname[4] = {"qc/ITS/MO/ITSFHR", "qc/ITS/MO/ITSFHR", "qc/ITS/MO/ITSFHR", "qc/ITS/MO/ITSFHR"};

  //set the task name
  switch(opt){
    case 1: {// fake-hit
      taskname[0] = "qc/ITS/MO/ITSFHR";//L0T, L0B
      taskname[1] = "qc/ITS/MO/ITSFHR";//L1T, L1B
      taskname[2] = "qc/ITS/MO/ITSFHR"; //L2T
      taskname[3] = "qc/ITS/MO/ITSFHR"; //L2B
      break;
    }

    case 2: { //thr scan
      taskname[0] = "qc/ITS/MO/ITSTHRTask0";
      taskname[1] = "qc/ITS/MO/ITSTHRTask1";
      taskname[2] = "qc/ITS/MO/ITSTHRTask2T";
      taskname[3] = "qc/ITS/MO/ITSTHRTask2B";
      break;
    }

    default: break;
  }//end of switch



  //CCDB api initialization
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init(ccdbport.c_str());

  //set variables (run interval of timestamp interval)
  array<string,2> runts1;
  array<string,2> runts2;

  //Decide how many different elements are needed
  int nListElements = 2;
  switch(opt){
    case 1: nListElements = 2; break;
    case 2: nListElements = 3; break;
    default: nListElements = 2;
  }
  if(adderrordata && opt==1)
    nListElements+=2;//2 because we add trigger and error plots

  //Output file
  string layername;
  if(layernum==-1)
    layername = "all-IB-layers";
  else
    layername = Form("Layer%d",layernum);

  string suffix = "run";
  string optname = GetOptName(opt);

  //Download depending on the option (opt)
  switch(opt){
    case 1: {

      //run interval definition
      cout<<"Finding runs in: "<<taskname[0]<<"/Occupancy/Layer0/Layer0ChipStave"<<endl;
      runts2 = GetLastRunWithTS(ccdbApi, taskname[0], "Occupancy/Layer0/Layer0ChipStave"); //take a random object name since run-list is the same.
      runts1 = GetRunWithTS24hAgo(ccdbApi, taskname[0], "Occupancy/Layer0/Layer0ChipStave", runts2[0]);
      cout<<"Run interval selected:       "<<runts1[1]<<"-"<<runts2[1]<<endl;
      cout<<"Timestamp interval selected: "<<runts1[0]<<"-"<<runts2[0]<<endl;

      //output file
      outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),runts1[1].c_str(), suffix.c_str(), runts2[1].c_str(), adderrordata? "_w_error_and_trig_data":""), "RECREATE");
      outputfile->cd();

      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum],  objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum);
              break;
            }

            case 1: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9){
                    Download(1, ccdb, ccdbApi, myname, taskname[layernum+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum);
                  }
                  else{
                    Download(1, ccdb, ccdbApi, myname, taskname[layernum], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum);
                  }
                }
                else {
                  Download(1, ccdb, ccdbApi, myname, taskname[layernum], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum);
                }

              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum);
              break;
            }

            case 3: {//error files
              string objname = "General/TriggerVsFeeid";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum);
              break;
            }
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){

        vector<string> goodrunlist = GetGoodRunList(ccdbApi, runts1[1], runts2[1], "Fhr");

        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, myname, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay);
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9){
                      Download(1, ccdb, ccdbApi, myname, taskname[ilay+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay);
                    }
                    else{
                      Download(1, ccdb, ccdbApi, myname, taskname[ilay], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay);
                    }
                  }
                  else{
                    Download(1, ccdb, ccdbApi, myname, taskname[ilay], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]),ilay);
                  }
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[0], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), 0);
              break;
            }

            case 3: {//trigger files
              string objname = "General/TriggerVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[0], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), 0);
              break;
            }
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }//end case 1

    case 2: {//thresholds
      //run interval definition
      cout<<"Finding runs in: "<<taskname[0]<<"/Threshold/Layer0/Threshold_Vs_Chip_and_Stave"<<endl;
      runts2 = GetLastRunWithTS(ccdbApi, taskname[0], "Threshold/Layer0/Threshold_Vs_Chip_and_Stave"); //take a random object name since run-list is the same.
      runts1 = GetRunWithTS24hAgo(ccdbApi, taskname[0], "Threshold/Layer0/Threshold_Vs_Chip_and_Stave", runts2[0]);
      vector<string> goodrunlist = GetGoodRunList(ccdbApi, runts1[1], runts2[1], "Thr");
      cout<<"Run interval selected:       "<<runts1[1]<<"-"<<runts2[1]<<endl;
      cout<<"Timestamp interval selected: "<<runts1[0]<<"-"<<runts2[0]<<endl;

      outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),runts1[1].c_str(), suffix.c_str(), runts2[1].c_str(), adderrordata? "_w_error_and_trig_data":""), "RECREATE");
      outputfile->cd();

      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]),layernum);
              break;
            }

            case 1: {
              string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum], objname, runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum);
              break;
            }

            case 2: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9)
                    Download(1, ccdb, ccdbApi, myname, taskname[layernum+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum);
                  else
                    Download(1, ccdb, ccdbApi, myname, taskname[layernum], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum);
                }
                else{
                  Download(1, ccdb, ccdbApi, myname, taskname[layernum],  Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), runts1[1], runts2[1], vector<string>(), stol(runts1[0]), stol(runts2[0]), layernum);
                }
              }
              break;
            }

            default: break;
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, myname, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay);
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, myname, taskname[ilay], objname, runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay);
              }
              break;
            }

            case 2: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9)
                      Download(1, ccdb, ccdbApi, myname, taskname[ilay+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay);
                    else
                      Download(1, ccdb, ccdbApi, myname, taskname[ilay], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay);
                  }
                  else{
                    Download(1, ccdb, ccdbApi, myname, taskname[ilay], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), runts1[1], runts2[1], goodrunlist, stol(runts1[0]), stol(runts2[0]), ilay);
                  }
                }
              }
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }// end of case 2

  }//end switch

  outputfile->Close();
  delete outputfile;

  return 1;
}

//
// Expert mode
//
bool RunExpert(auto *ccdb, string myname, int opt){

  //choose IB, OB or both
  int IBorOB;
  cout << endl;
  cout << endl;
  cout << "Choose beteen IB (0), OB (1) or both (2, = all layers of IB and OB)" << endl;
  cin>> IBorOB;

  //Choose the layer number
  int layernum;
  cout<<endl;

  switch(IBorOB){
  case 0: {
    cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
    cin >>layernum;
    break;
  }
  case 1: {
    cout<<"Enter the layer number [put -1 for all OB layers]"<<endl;
    cin >>layernum;
    break;
  }
  case 2: {
    layernum=-1;
    break;
  }
  default: break;
  }
  //  cin>>layernum;

  //taskname
  string taskname[8]    = {"qc/ITS/MO/ITSFHR", "qc/ITS/MO/ITSFHR", "qc/ITS/MO/ITSFHR", "qc/ITS/MO/ITSFHR",  "qc/ITS/MO/FHRTask",  "qc/ITS/MO/FHRTask",  "qc/ITS/MO/FHRTask",  "qc/ITS/MO/FHRTask"};

  //set the task name
  switch(IBorOB){
  case 0: {//Inner Barrel
    switch(opt){
    case 1: {// fake-hit
      taskname[0] = "qc/ITS/MO/ITSFHR";//L0T, L0B
      taskname[1] = "qc/ITS/MO/ITSFHR";//L1T, L1B
      taskname[2] = "qc/ITS/MO/ITSFHR"; //L2T
      taskname[3] = "qc/ITS/MO/ITSFHR"; //L2B
      break;
    }

    case 2: { //thr scan
      taskname[0] = "qc/ITS/MO/ITSTHRTask0";
      taskname[1] = "qc/ITS/MO/ITSTHRTask1";
      taskname[2] = "qc/ITS/MO/ITSTHRTask2T";
      taskname[3] = "qc/ITS/MO/ITSTHRTask2B";

      break;
    }
    default: break;
    }//end of switch
    break;
  }
  case 1: { //Outer Barrel
    switch(opt){
    case 1: {// fake-hit
      taskname[4] = "qc/ITS/MO/FHRTask";//L3
      taskname[5] = "qc/ITS/MO/FHRTask";//L4
      taskname[6] = "qc/ITS/MO/FHRTask";//L5
      taskname[7] = "qc/ITS/MO/FHRTask";//L6

      break;
    }

    case 2: { //thr scan
      //path name to be changed for OB
      taskname[4] = "qc/ITS/MO/ITSTHRTask0";
      taskname[5] = "qc/ITS/MO/ITSTHRTask1";
      taskname[6] = "qc/ITS/MO/ITSTHRTask2T";
      taskname[7] = "qc/ITS/MO/ITSTHRTask2B";

      break;
    }
    default: break;
    }//end of switch
  }
  case 2: { //Inner and Outer Barrel
    switch(opt){
    case 1: {// fake-hit
      taskname[0] = "qc/ITS/MO/FHRTask";//L0
      taskname[1] = "qc/ITS/MO/FHRTask";//L1
      taskname[2] = "qc/ITS/MO/FHRTask";//L2
      taskname[3] = "qc/ITS/MO/FHRTask";//not used for FHR
      taskname[4] = "qc/ITS/MO/FHRTask";//L3
      taskname[5] = "qc/ITS/MO/FHRTask";//L4
      taskname[6] = "qc/ITS/MO/FHRTask";//L5
      taskname[7] = "qc/ITS/MO/FHRTask";//L6

      break;
    }

    case 2: { //thr scan
      // path name to be changed
      taskname[0] = "qc/ITS/MO/ITSTHRTask0";
      taskname[1] = "qc/ITS/MO/ITSTHRTask1";
      taskname[2] = "qc/ITS/MO/ITSTHRTask2T";
      taskname[3] = "qc/ITS/MO/ITSTHRTask2B";
      taskname[4] = "qc/ITS/MO/ITSTHRTask0";
      taskname[5] = "qc/ITS/MO/ITSTHRTask1";
      taskname[6] = "qc/ITS/MO/ITSTHRTask2T";
      taskname[7] = "qc/ITS/MO/ITSTHRTask2B";

      break;
    }
    default: break;
    }//end of switch
  }
  default: break;
  }//end of switch

  //Ask whether attach a error report to the results
  bool adderrordata = false;
  string erranswer = "n";
  if(opt==1){//only is user selects to download FHR
    cout<<endl;
    cout<<"Error and trigger plots needed? [y/n]"<<endl;
    cin>>erranswer;
    if(erranswer=="y")
      adderrordata = true;
  }

  //Chose about run number or time interval
  int choice;
  cout<<endl;
  cout<<"Download data - Choose an option:\n"<<endl;
  cout<<"1. Enter run numbers"<<endl;
  cout<<"2. Enter timestamps"<<endl;
  cout<<"[Type 1 or 2]: ";
  cin>>choice;

  //set variables (run interval of timestamp interval)
  string run1="0", run2="0";
  time_t ts_start=0, ts_end=0;
  vector<string> nums; //can contain selected run interval or timestamp interval
  switch(choice){
    case 1: {
      cout<<"\n"<<"Run interval [just run number, NO ZEROs in front]"<<endl;
      cout<<"Enter first  run: ";
      cin>> run1;
      cout<<"Enter second run: ";
      cin>> run2;
      nums.push_back(run1);
      nums.push_back(run2);
      break;
    }

    case 2:{
      //Enter time period of the data to be downloaded
      struct tm datetime_start;
      struct tm datetime_end;
      string day, month, hour, min, sec;
      int year;
      cout<<"\n"<<"Time period [format (hour 24H): DD MM YYYY HH MM SS]"<<endl;
      cout<<"Enter start date and time: ";
      cin>> day >> month >> year >> hour >> min >> sec;
      string datetime1 = Form("%d%s%s_%s%s%s",year,month.c_str(),day.c_str(),hour.c_str(),min.c_str(),sec.c_str());
      datetime_start.tm_mday = stoi(day);
      datetime_start.tm_mon = stoi(month)-1;//-1 because january = 0;
      datetime_start.tm_year = year-1900;
      datetime_start.tm_hour = stoi(hour);
      datetime_start.tm_min = stoi(min);
      datetime_start.tm_sec = stoi(sec);
      datetime_start.tm_isdst = -1; //-1 means DST unknown
      ts_start = mktime(&datetime_start);//get timestamp start
      ts_start*=1000; //convert in milliseconds

      cout<<"Enter end date and time  : ";
      cin>> day >> month >> year >> hour >> min >> sec;
      string datetime2 = Form("%d%s%s_%s%s%s",year,month.c_str(),day.c_str(),hour.c_str(),min.c_str(),sec.c_str());
      datetime_end.tm_mday = stoi(day);
      datetime_end.tm_mon = stoi(month)-1;//-1 because january = 0;
      datetime_end.tm_year = year-1900;
      datetime_end.tm_hour = stoi(hour);
      datetime_end.tm_min = stoi(min);
      datetime_end.tm_sec = stoi(sec);
      datetime_end.tm_isdst = -1; //-1 means DST unknown
      ts_end = mktime(&datetime_end);//get timestamp start
      ts_end*=1000; //convert in milliseconds
      printf("Timestamp start (ms): %ld\n", (long) ts_start);
      printf("Timestamp end   (ms): %ld\n", (long) ts_end);
      nums.push_back(datetime1);
      nums.push_back(datetime2);
      break;
    }

    default: break;
  }

  //Decide how many different elements are needed
  int nListElements = 2;
  switch(opt){
    case 1: nListElements = 2; break;
    case 2: nListElements = 3; break;
    default: nListElements = 2;
  }
  if(adderrordata)
    nListElements+=2;

  //CCDB api initialization
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init(ccdbport.c_str());

  //Output file
  string layername;
  if(layernum==-1)
    if(IBorOB==0)    layername = "all-IB-layers";
    else if (IBorOB==1) layername = "all-OB-layers";
    else layername = "all-layers";
  else
    layername = Form("Layer%d",layernum);

  int layernumEff=layernum;
  if(IBorOB==1) layernumEff=layernum+1;

  int ilayMin = 0;
  int ilayMax = 2;
  if (IBorOB==1) {
    ilayMin = 3;
    ilayMax =6;
  }
  else if (IBorOB==2) {
    ilayMin = 0;
    ilayMax =6;
  }

  string suffix;
  switch(choice){
    case 1: suffix = "run"; break;
    case 2: suffix = ""; break;
    default: suffix="";
  }
  string optname = GetOptName(opt);
  outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),nums[0].c_str(), suffix.c_str(), nums[1].c_str(), adderrordata? "_w_error_and_trig_data":""), "RECREATE");
  outputfile->cd();

  //Download depending on the option (opt)
  switch(opt){
    case 1: {
      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end,layernum);
              break;
            }

            case 1: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9){//L2B
                    Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
                  }
                  else{
                    Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
                  }
                }
                else{
                  Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
              break;
            }

            case 3: {//error files
              string objname = "General/TriggerVsFeeid";
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
              break;
            }

            default: break;
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
                cout<<"\nAll data in "<<taskname[ilayEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay);
              }
              break;
            }

            case 1: {
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9){//L2B
                      Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff+1], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay);
                    }
                    else{
                      Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay);
                    }
                  }
                  else{
                    Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff],  Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, ilay);
                  }
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[0], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0);
              break;
            }

            case 3: {//trigger files
              string objname = "General/TriggerVsFeeid";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[0], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0);
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }//end case 1

    case 2: {// thresholds
      vector<string> goodrunlist = GetGoodRunList(ccdbApi, run1, run2, "Thr"); //good run list
      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff],  objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end,layernum);
              break;
            }

            case 1: {
              string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernumEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], objname, run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
              break;
            }

            case 2: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                if(layernum==2){
                  if(istave>9)
                    Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
                  else
                    Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
                }
                else{
                  Download(choice, ccdb, ccdbApi, myname, taskname[layernumEff], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",layernum,istave), run1, run2, vector<string>(), (long)ts_start, (long)ts_end, layernum);
                }

              }
              break;
            }

            default: break;
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
	      //for(int ilay=0; ilay<=2; ilay++){
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilayEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff], objname, run1, run2, goodrunlist, (long)ts_start, (long)ts_end,ilay);
              }
              break;
            }

            case 1: {
	      // for(int ilay=0; ilay<=2; ilay++){
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilayEff]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff],  objname, run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay);
              }
              break;
            }

            case 2: {
	      //              for(int ilay=0; ilay<=2; ilay++){
              for(int ilay=ilayMin; ilay<=ilayMax; ilay++){
		int ilayEff=ilay;
		if (IBorOB==1) ilayEff=ilay+1;
		else if (IBorOB==2 && ilay>=3) ilayEff=ilay+1;
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  if(ilay==2){
                    if(istave>9)
                      Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff+1], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay);
                    else
                      Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff], Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay);
                  }
                  else{
                    Download(choice, ccdb, ccdbApi, myname, taskname[ilayEff],  Form("DeadPixel/Layer%d/Stave%d/HIC0/DeadPixelHITMAP",ilay,istave), run1, run2, goodrunlist, (long)ts_start, (long)ts_end, ilay);
                  }
                }
              }
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }// end of case 2

    case 3: {//download a tree (only an example now!)
      Download(choice, ccdb, ccdbApi, myname, "qc/ITS/QCFHRTest", "Occupancy/PixelTree", run1, run2, vector<string>(), (long)ts_start, (long)ts_end, 0);
      break;
    }//end of case 3
  }//end switch

  outputfile->Close();
  delete outputfile;

  return 1;
}

//
// Return run list of GOOD runs = runs in all CCDB paths for TH2
//
vector<string> GetGoodRunList(o2::ccdb::CcdbApi ccdbApi, string run1, string run2, string runtype){//run type can be "Thr" or "Fhr"

  string objnames[4] = {"Threshold/Layer0/Threshold_Vs_Chip_and_Stave", "Threshold/Layer1/Threshold_Vs_Chip_and_Stave", "Threshold/Layer2/Threshold_Vs_Chip_and_Stave", "Threshold/Layer2/Threshold_Vs_Chip_and_Stave"};
  string tasknames[4] = {"qc/ITS/MO/ITSTHRTask0", "qc/ITS/MO/ITSTHRTask1", "qc/ITS/MO/ITSTHRTask2T", "qc/ITS/MO/ITSTHRTask2B"};

  //in case of FHR
  if(runtype=="Fhr"){
    objnames[0] = "Occupancy/Layer0/Layer0ChipStave";
    objnames[1] = "Occupancy/Layer1/Layer1ChipStave";
    objnames[2] = "Occupancy/Layer2/Layer2ChipStave";
    objnames[3] = "Occupancy/Layer2/Layer2ChipStave";

    tasknames[0] = "qc/ITS/MO/ITSFHR";
    tasknames[1] = "qc/ITS/MO/ITSFHR";
    tasknames[2] = "qc/ITS/MO/ITSFHR";
    tasknames[3] = "qc/ITS/MO/ITSFHR";
  }

  vector<vector<string>> runlists;
  for(int i=0; i<4; i++){//get run list between run1 and run2 for all the paths defined above
    string objectlist = ccdbApi.list(tasknames[i] + "/" + objnames[i],false,"text/plain");
    stringstream ss(objectlist);
    string word;
    vector<string>runlist;
    while(ss>>word){
      if(word=="Created:"){// take the one related to file creation
        ss>>word;
        //alltimestampsALT.push_back(word);
      }
      if(word=="Run"){
        ss>>word;
        ss>>word;
        if(stoi(word)>=stoi(run1) && stoi(word)<=stoi(run2)){
	        if(std::find(runlist.begin(), runlist.end(), word) == runlist.end())//if element doesn't exist already
		        runlist.push_back(word);
	      }
        if(stoi(word)==stoi(run1)) break;
      }
    }
    runlists.push_back(runlist);
  }

  //Extract runs present in all runlists
  vector<string> runsfiltered;
  for(int il=0; il<(int)runlists[0].size(); il++){
    string run = runlists[0][il];
    int cntrun = 1;
    if((int)runlists[2].size()==0) cntrun++;// in case half L2 is missing
    if((int)runlists[3].size()==0) cntrun++;// in case half L2 is missing
    //look for this run in the other lists
    for(int cntlist=1; cntlist<4; cntlist++){
      for(int il2=0; il2<(int)runlists[cntlist].size(); il2++){
        if(run==runlists[cntlist][il2]) {
          cntrun++;
        }
      }
    }
    if(cntrun==4){
      runsfiltered.push_back(run);
      cout<<"Run: "<<run<<endl;
    }
  }

  return runsfiltered;
}


//
// Download data based on timestamps
//
void DownloadTimestamps(auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, long int ts_start, long int ts_end, int lnum){

  //Extract all the time stamps of the object
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  stringstream ss(objectlist);
  string word;
  vector <string> timestamps;
  while(ss>>word){
    if(word=="Created:"){// take the one related to file creation
      ss>>word;
      timestamps.push_back(word);
    }
  }

  //filter all the timestamps to get the ones within start and end dates (get most recent!!)
  vector <long int> timestamps_selperiod;
  int counter=0;
  for(int its=0; its<(int)timestamps.size()-1; its++){
    if(stol(timestamps[its])>=ts_start && stol(timestamps[its])<=ts_end){
      counter++;
      if(counter==1) {
        if(its==0)
          timestamps_selperiod.push_back(stol(timestamps[its]));//the first (most recent) is taken by definition if it is really the first in the list
        else if(stol(timestamps[its-1])-stol(timestamps[its])>65000) timestamps_selperiod.push_back(stol(timestamps[its]));//if not the real first, check if it is the most recent (in this way the user can insert whatever date interval)
      }
      if(stol(timestamps[its])-stol(timestamps[its+1])>65000) timestamps_selperiod.push_back(stol(timestamps[its+1]));
    }

    else if(stol(timestamps[its])>ts_end) continue;
    else break;
  }
  cout<<"\n"<<"Selected timestamps: "<<endl;
  for(int i=0; i<(int)timestamps_selperiod.size();i++){
    cout<<timestamps_selperiod[i]<<endl;
  }

  bool isperstave = 0;
  if(objname.find("HITMAP")!=string::npos) isperstave = 1;
  GetListOfHisto(ccdb, myname, taskname, objname, timestamps_selperiod, vector<long int>(), lnum, 0, isperstave, vector<int>(), vector<int>());

  timestamps_selperiod.clear();
}


//
// Download data based on run numbers - available in metadata from 03/07/2019 21.49 (run 582 --> fake hit scan)
//
void DownloadRuns(auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, int lnum){

  //Extract all the time stamps and run numbers of the object
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  string objectlist2 = " ";
  /* alt*
  if(tasknamealternative!=taskname){
    objectlist2 = ccdbApi.list(tasknamealternative + "/" + objname,false,"text/plain");
  }
  */
  string objectlistL2B = " ";
  if(/*objname.find("Layer2ChipStave")!=string::npos || objname.find("ErrorFile")!=string::npos || objname.find("TriggerFile")!=string::npos
     ||*/ objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos){
    objectlistL2B = ccdbApi.list(Form("qc/ITS/MO/ITS%sTask2B/%s",objname.find("Chip_and_Stave")!=string::npos ? "THR":"Raw",objname.c_str()), false, "text/plain");
  }
  cout<<endl;
  cout<<endl;
  cout<<"Ready to get files from "<<taskname<<"/"<<objname<<endl;
  /* alt*
  if(tasknamealternative!=taskname){
    cout<<"... And eventually from (alternative path for backward compatibility): "<< tasknamealternative<<"/"<<objname<<endl;
  }
  */
  stringstream ss(objectlist);
  //alt*  stringstream ss2(objectlist2);
  stringstream ss3(objectlistL2B);
  string word;
  vector<string> alltimestamps, timestamps, runs;
  //alt*  vector<string> alltimestampsALT, timestampsALT, runsALT;//for alternative path
  vector<string> alltimestampsL2B, timestampsL2B, runsL2B;//for L2B

  //filter string from alternative path (at the moment: only L1 from before run 300134)
  /*
  while(ss2>>word){
    if(word=="Created:"){// take the one related to file creation
      ss2>>word;
      alltimestampsALT.push_back(word);
    }
    if(word=="Run"){
      ss2>>word;
      ss2>>word;
      runsALT.push_back(word);
      timestampsALT.push_back(alltimestampsALT[alltimestampsALT.size()-1]);//this keep only the timestamps connected to a run number
      if(stoi(word)==stoi(run1)) break;
    }
  }
  */

  //filter normal path but correct timestamps with the one from alternative path (if needed)
  while(ss>>word){
    if(word=="Created:"){// take the one related to file creation
      ss>>word;
      alltimestamps.push_back(word);
    }
    if(word=="Run"){
      ss>>word;
      ss>>word;
      /* alt*
      if(tasknamealternative!=taskname && stol(word)<300134){//for L1 backward compatibility
        for(int it=0; it<(int)runsALT.size(); it++){
          if(stol(runsALT[it])<300134){
            runs.push_back(runsALT[it]);
            timestamps.push_back(timestampsALT[it]);
          }
        }
        break;
      }
      else if(tasknamealternative!=taskname && stol(word)>=300134){
        runs.push_back(word);
        timestamps.push_back(alltimestamps[alltimestamps.size()-1]);//this keep only the timestamps connected to a run number
      }
      else if(tasknamealternative==taskname){
        runs.push_back(word);
        timestamps.push_back(alltimestamps[alltimestamps.size()-1]);//this keep only the timestamps connected to a run number
      }
      */
      runs.push_back(word);
      timestamps.push_back(alltimestamps[alltimestamps.size()-1]);//this keep only the timestamps connected to a run number
      if(stoi(word)==stoi(run1)) break;
    }
  }


  //save all runs and timestamps for L2B
  while(ss3>>word){
    if(word=="Created:"){// take the one related to file creation
      ss3>>word;
      alltimestampsL2B.push_back(word);
    }
    if(word=="Run"){
      ss3>>word;
      ss3>>word;
      runsL2B.push_back(word);
      timestampsL2B.push_back(alltimestampsL2B[alltimestampsL2B.size()-1]);//this keep only the timestamps connected to a run number
      if(stoi(word)==stoi(run1)) break;
    }
  }



  //filter all runs to get the ones within run1 and run2 (get most recent for each run!!) --> take from good run list if defined
  vector <long int> timestamps_selperiod;
  vector <int> runs_selperiod;
  int counter=0;
  for(int irun=0; irun<(int)runs.size(); irun++){
    if((stoi(runs[irun])>=stoi(run1) && stoi(runs[irun])<=stoi(run2)) && (goodrunlist.size()==0 || std::find(goodrunlist.begin(), goodrunlist.end(), runs[irun]) != goodrunlist.end())){
      counter++;
      if(counter==1){
        runs_selperiod.push_back(std::stoi(runs[irun]));
        timestamps_selperiod.push_back(std::stol(timestamps[irun]));
      }
      else{
        if(runs[irun]==runs[irun-1]) continue;
        else{
          runs_selperiod.push_back(std::stoi(runs[irun]));
          timestamps_selperiod.push_back(std::stol(timestamps[irun]));
        }
      }
    }
  }

  //do the same for L2B (when it's its turn)
  vector <long int> timestamps_selperiodL2B;
  vector <int> runs_selperiodL2B;
  counter=0;
  for(int irun=0; irun<(int)runsL2B.size(); irun++){
    if(stoi(runsL2B[irun])>=stoi(run1) && stoi(runsL2B[irun])<=stoi(run2) && (goodrunlist.size()==0 || std::find(goodrunlist.begin(), goodrunlist.end(), runsL2B[irun]) != goodrunlist.end())){
      counter++;
      if(counter==1){
        runs_selperiodL2B.push_back(std::stoi(runsL2B[irun]));
        timestamps_selperiodL2B.push_back(std::stol(timestampsL2B[irun]));
      }
      else{
        if(runsL2B[irun]==runsL2B[irun-1]) continue;
        else{
          runs_selperiodL2B.push_back(std::stoi(runsL2B[irun]));
          timestamps_selperiodL2B.push_back(std::stol(timestampsL2B[irun]));
        }
      }
    }
  }

  cout<<"\n"<<"Selected runs and corresponding timestamps: "<<endl;
  for(int i=0; i<(int)runs_selperiod.size();i++){
    cout<<"run"<<runs_selperiod[i]<<" - "<<timestamps_selperiod[i]<<endl;
  }

  bool isperstave = 0;
  if(objname.find("HITMAP")!=string::npos) isperstave = 1;
  GetListOfHisto(ccdb, myname, taskname, objname, timestamps_selperiod, timestamps_selperiodL2B, lnum, 1, isperstave, runs_selperiod, runs_selperiodL2B);

  timestamps_selperiod.clear();
  runs_selperiod.clear();
  timestamps.clear();
  runs.clear();
  alltimestamps.clear();

}

//
// Get correct timestamp
//
string GetCorrectTS(string selrun, vector<string> runs, vector<string> timestamps){

  for(int it=0; it<(int)timestamps.size(); it++){
    if(selrun==runs[it]){
      return timestamps[it];
    }
  }
  return "nots";
}

//
//Get list of histogram inside an object
//
bool GetListOfHisto(auto* ccdb, string myname, string taskname, string objname, vector<long int> timestamps, vector<long int> timestamps2, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers, vector<int>runnumbers2){

  //Getting root files from the database and write them to file
  cout<<"\n"<<"... Getting files from the database"<<endl;

  string stvnum = "0";
  if(isperstave){
    stvnum = objname.substr(objname.find("Stave")+5,2);
    if(stvnum.find("/")!=string::npos)
      stvnum = objname.substr(objname.find("Stave")+5,1);
  }

 // cout<<"timestamps  size: "<<timestamps.size()<<endl;
 // cout<<"timestamps2 size: "<<timestamps2.size()<<endl;

  //in case L2T is missing (excluded)
  if((int)timestamps.size()==0 && (objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos)){
    for(int i=0; i<(int)timestamps2.size();i++){

      //MonitorObject *monitor = ccdb->retrieve(taskname, objname, timestamps[i]);

      auto monitor = ccdb->retrieveMO(Form("qc/ITS/MO/ITS%sTask2B",objname.find("Chip_and_Stave")!=string::npos ? "THR":"Raw"), objname, timestamps2.size()>0 ? timestamps2[i]:timestamps[i]);

      if (monitor == nullptr) {
        cerr << myname << ": failed to get MonitorObject for timestamp: " << timestamps[i]<< endl;
        return 0;
      }

      TObject *obj = nullptr;
      obj = monitor->getObject();
      monitor->setIsOwner(false);
      string c = obj->ClassName();
      TH2 *h2s = 0x0;

      //if(strstr(c,"TH2")!=nullptr){
      if(c.find("TH2")!=string::npos){
        string histname = "";
        histname = Form("h2_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers2[i]) : "", timestamps2[i]);

        h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
        outputfile->cd();
        h2s->Write();
      }
    }
  }

  //General case
  for(int i=0; i<(int)timestamps.size();i++){

    //MonitorObject *monitor = ccdb->retrieve(taskname, objname, timestamps[i]);

    auto monitor = ccdb->retrieveMO(taskname, objname, timestamps[i]);

    if (monitor == nullptr) {
      cerr << myname << ": failed to get MonitorObject for timestamp: " << timestamps[i]<< endl;
      return 0;
    }

    TObject *obj = nullptr;
    obj = monitor->getObject();
    monitor->setIsOwner(false);
    //for L2B only
    TObject *obj2 = nullptr;
    if(/*objname.find("Layer2ChipStave")!=string::npos || objname.find("ErrorFile")!=string::npos || objname.find("TriggerFile")!=string::npos
       || */objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos){
      auto monitor2 = ccdb->retrieveMO(Form("qc/ITS/MO/ITS%sTask2B",objname.find("Chip_and_Stave")!=string::npos ? "THR":"Raw"), objname, timestamps2.size()>0 ? timestamps2[i]:timestamps[i]);
      obj2 = monitor2->getObject();
      monitor2->setIsOwner(false);
    }
    string c = obj->ClassName();
    TH2 *h2s = 0x0;
    TH2 *h2sbis = 0x0; //for L2B only
    TH1 *h1s = 0x0;
    THnSparse *hSp = 0x0;
    TTree *tree;

    //if(strstr(c,"TH1")!=nullptr){
    if(c.find("TH1")!=string::npos){
      string histname = "";
      histname = Form("h1_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      h1s = dynamic_cast<TH1*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      h1s->Write();
    }
    //if(strstr(c,"TH2")!=nullptr){
    if(c.find("TH2")!=string::npos){
      string histname = "";
      if(objname.find("Error")!=string::npos)
        histname = Form("h2_err%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      else if(objname.find("Trigger")!=string::npos)
        histname = Form("h2_trg%s_%ld", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      else
        histname = Form("h2_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
      if(/*objname.find("Layer2ChipStave")!=string::npos || objname.find("ErrorFile")!=string::npos || objname.find("TriggerFile")!=string::npos
        || */objname.find("Layer2/Threshold_Vs_Chip_and_Stave")!=string::npos || objname.find("Layer2/DeadPixel_Vs_Chip_and_Stave")!=string::npos){
        h2sbis = dynamic_cast<TH2*>(obj2->Clone(Form("%s_L2B",histname.c_str())));
        h2s->Add(h2sbis);
      }
      outputfile->cd();
      h2s->Write();
    }

    if(c.find("THnSparse")!=string::npos){
      string histname = "";
      histname = Form("hsparse_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      hSp = dynamic_cast<THnSparse*>(obj->Clone(histname.c_str()));
      outputfile->cd();
      hSp->Write();
    }

    if(c.find("TTree")!=string::npos){
      string treename = Form("ttree_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);;
      //obj->Print();
      tree = dynamic_cast<TTree*>(obj->Clone(treename.c_str()));
      //tree->Print();
      outputfile->cd();
      tree->Write();
    }

    delete h1s;
    delete h2s;
    delete hSp;
    delete obj;
  }

  return 1;
}

//
// Download depending on the choice
//
bool Download(int choice, auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2, vector<string> goodrunlist, long int ts_start, long int ts_end, int lnum){

  switch(choice){
    case 1: {
      DownloadRuns(ccdb, ccdbApi, myname, taskname, objname, run1, run2, goodrunlist, lnum);//download with runs
      break;//download runs data
    }

    case 2: {
      DownloadTimestamps(ccdb, ccdbApi, myname, taskname, objname, ts_start, ts_end, lnum); //download with timestamps
      break;
    }

    default: return 0;
  }

  return 1;

}

//
// Get name depending on the option
//
string GetOptName(int opt){
  switch(opt){
    case 1: return "FHRMAPS_HITMAPS";
    case 2: return "THRMAPS_DEADPIXMAPS";
    case 3: return "NOISYPIX_TREE";
    default: return "0";
  }
}

//
// Get name to assign to the output list
//
string GetListName(int opt, int ilist){
  switch(opt){
    case 1:{
      if(!ilist)
        return "fhrmaps";
      else
        return "hitmaps";
    }
    default: return "0";
  }
}
