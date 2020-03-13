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

using namespace std;
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;

//functions to download data
string GetCorrectTS(string selrun, vector<string> runs, vector<string> timestamps);
array<string,2> GetLastRunWithTS(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname);
array<string,2> GetRunWithTS24hAgo(o2::ccdb::CcdbApi ccdbApi, string taskname, string objname, string timestamp);
bool RunShifter(auto *ccdb, string myname);
bool RunExpert(auto *ccdb, string myname);
void DownloadTimestamps(auto *ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, long int ts_start, long int ts_end, int lnum);
void DownloadRuns(auto *ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string tasknamealternative, string objname, string run1, string run2, int lnum);
bool GetListOfHisto(auto* ccdb, string myname, string taskname, string tasknamealternative, string objname, vector<long int> timestamps, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers);
bool Download(int choice, auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string tasknamealternative, string objname, string run1, string run2, long int ts_start, long int ts_end, int lnum);
string GetOptName(int opt);
string GetListName(int opt, int ilist);

const int nStavesInLay[3] = {12, 16, 20};
TFile *outputfile;


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

  ccdb->connect("ccdb-test.cern.ch:8080", "", "", "");

  if(strcmp(argv[1],"expert")==0)
    RunExpert(ccdb, myname);
  else if (strcmp(argv[1],"shifter")==0)
    RunShifter(ccdb, myname);


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

  cout<<"now: "<<stamp_int_actual<<"  24hago: "<<stamp_int_24hago<<endl;

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
bool RunShifter(auto *ccdb, string myname){
  //Choose what to download
  int opt;
  cout<<endl;
  cout<<endl;
  cout<<"Choose what to download:"<<endl;
  cout<<"1. FakeHitRate runs"<<endl;
  cout<<endl;
  cout<<"Enter the option: ";
  cin>>opt;
  if(opt<1 || opt>1){
    cout<<"Invalid option"<<endl;
    return 0;
  }

  //Choose the layer number
  int layernum;
  cout<<endl;
  cout<<endl;
  cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
  cin>>layernum;

  //taskname
  string taskname[3] = {"qc/ITS/ITSRawTask","qc/ITS/ITSRawTask", "qc/ITS/ITSRawTask"};

  //Choose the side: top or bottom
  string side;
  bool adderrordata = false;
  cout<<endl;
  cout<<"Top or Bottom? [T/B]"<<endl;
  cin>>side;
  if(side=="T" || side=="t"){
    adderrordata = true;
  }
  else if(side=="B" || side=="b"){
    if(layernum<0) {
      taskname[0] = "qc/ITS/ITSRawTaskIBB2";
      taskname[1] = "qc/ITS/ITSRawTaskIBB3";
      taskname[2] = "qc/ITS/ITSRawTaskIBB1";
    }
    else if(layernum==1) taskname[1] = "qc/ITS/ITSRawTaskIBB3";
    else if(layernum==2) taskname[2] = "qc/ITS/ITSRawTaskIBB1";
    else taskname[0] = "qc/ITS/ITSRawTaskIBB2";
    adderrordata = false;
  }

  //CCDB api initialization
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init("ccdb-test.cern.ch:8080");

  //set variables (run interval of timestamp interval)
  array<string,2> runts1;
  array<string,2> runts2;

  //Decide how many different elements are needed
  int nListElements = 2;
  switch(opt){
    case 1: nListElements = 2; break;
    default: nListElements = 2;
  }
  if(adderrordata)
    nListElements+=1;//2 because we add trigger and error plots

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
      runts2 = GetLastRunWithTS(ccdbApi, taskname[0], "Occupancy/Layer0/Layer0ChipStave"); //take a random object name since run-list is the same.
      runts1 = GetRunWithTS24hAgo(ccdbApi, taskname[0], "Occupancy/Layer0/Layer0ChipStave", runts2[0]);

      //output file
      outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),runts1[1].c_str(), suffix.c_str(), runts2[1].c_str(), adderrordata? "_w_error_data":""), "RECREATE");
      outputfile->cd();

      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum], taskname[layernum], objname, runts1[1], runts2[1], stol(runts1[0]), stol(runts2[0]),layernum);
              break;
            }

            case 1: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                Download(1, ccdb, ccdbApi, myname, taskname[layernum], taskname[layernum], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), runts1[1], runts2[1], stol(runts1[0]), stol(runts2[0]),layernum);
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorFile";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[layernum], taskname[layernum], objname, runts1[1], runts2[1], stol(runts1[0]), stol(runts2[0]),layernum);
              break;
            }
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
                Download(1, ccdb, ccdbApi, myname, taskname[ilay], taskname[ilay], objname, runts1[1], runts2[1], stol(runts1[0]), stol(runts2[0]),ilay);
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  Download(1, ccdb, ccdbApi, myname, taskname[ilay], taskname[ilay], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), runts1[1], runts2[1], stol(runts1[0]), stol(runts2[0]),ilay);
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorFile";
              cout<<"\nAll data in "<<taskname[0]+"/"+objname<<" between run"<<runts1[1]<<" and run"<<runts2[1]<<" are going to be downloaded."<<endl;
              Download(1, ccdb, ccdbApi, myname, taskname[0], taskname[0], objname, runts1[1], runts2[1], stol(runts1[0]), stol(runts2[0]),layernum);
              break;
            }
          }
        }//end loop on lists
      }//end if layernum==-1
    }//end case 1

    case 2: {

    }

  }//end switch

  outputfile->Close();
  delete outputfile;

  return 1;
}

//
// Expert mode
//
bool RunExpert(auto *ccdb, string myname){
  //Choose what to download
  int opt;
  cout<<endl;
  cout<<endl;
  cout<<"Choose what to download:"<<endl;
  cout<<"1. Fake-hit scan data"<<endl;
  cout<<"2. Threshold scan data"<<endl;
  cout<<endl;
  cout<<"Enter the option: ";
  cin>>opt;
  if(opt<1 || opt>2){
    cout<<"Invalid option. Doing nothing."<<endl;
    return 0;
  }

  //Choose the layer number
  int layernum;
  cout<<endl;
  cout<<endl;
  cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
  cin>>layernum;

  //taskname
  string taskname[3]    = {"qc/ITS/ITSRawTask", "qc/ITS/ITSRawTask", "qc/ITS/ITSRawTask"};
  string tasknamealt[3] = {"qc/ITS/ITSRawTask", "qc/ITS/ITSRawTask", "qc/ITS/ITSRawTask"};//alternative tasks (backward compatibility)

  //Choose the side: top or bottom
  string side;
  cout<<endl;
  cout<<"Top or Bottom? [T/B]"<<endl;
  cin>>side;

  //set the task name
  cout<<"OPT: "<<opt<<endl;
  switch(opt){
    case 1: {// fake-hit
      if(side=="B" || side=="b"){
        if(layernum<0) {
          taskname[0] = "qc/ITS/ITSRawTaskIBB2";
          taskname[1] = "qc/ITS/ITSRawTaskIBB3";
          taskname[2] = "qc/ITS/ITSRawTaskIBB1";

          tasknamealt[0] = "qc/ITS/ITSRawTaskIBB2";
          tasknamealt[1] = "qc/ITS/ITSRawTaskIBB2";
          tasknamealt[2] = "qc/ITS/ITSRawTaskIBB1";
        }
        else if(layernum==1) {taskname[1] = "qc/ITS/ITSRawTaskIBB3"; tasknamealt[1] = "qc/ITS/ITSRawTaskIBB2";}
        else if(layernum==2) {taskname[2] = "qc/ITS/ITSRawTaskIBB1"; tasknamealt[2] = "qc/ITS/ITSRawTaskIBB1";}
        else {taskname[0] = "qc/ITS/ITSRawTaskIBB2"; tasknamealt[0] = "qc/ITS/ITSRawTaskIBB2";}
      }
      break;
    }

    case 2: { //thr scan
      if(side=="B" || side=="b"){
        cout<<"HERE HERE HERE"<<endl;
        taskname[0] = "qc/ITS/THTest2";
        taskname[1] = "qc/ITS/THTest3";
        taskname[2] = "qc/ITS/THTest";

        tasknamealt[0] = "qc/ITS/THTest2";
        tasknamealt[1] = "qc/ITS/THTest3";
        tasknamealt[2] = "qc/ITS/THTest";
      }
      break;
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
  ccdbApi.init("ccdb-test.cern.ch:8080");

  //Output file
  string layername;
  if(layernum==-1)
    layername = "all-IB-layers";
  else
    layername = Form("Layer%d",layernum);
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
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], objname, run1, run2, (long)ts_start, (long)ts_end,layernum);
              break;
            }

            case 1: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, (long)ts_start, (long)ts_end, layernum);
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorFile";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], objname, run1, run2, (long)ts_start, (long)ts_end, layernum);
              break;
            }

            case 3: {//error files
              string objname = "General/TriggerFile";
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], objname, run1, run2, (long)ts_start, (long)ts_end, layernum);
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
                string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], objname, run1, run2, (long)ts_start, (long)ts_end, ilay);
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, (long)ts_start, (long)ts_end, ilay);
                }
              }
              break;
            }

            case 2: {//error files
              string objname = "General/ErrorFile";
              for(int ilay=0; ilay<=2; ilay++){
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], objname, run1, run2, (long)ts_start, (long)ts_end, ilay);
              }
              break;
            }

            case 3: {//error files
              string objname = "General/TriggerFile";
              for(int ilay=0; ilay<=2; ilay++){
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], objname, run1, run2, (long)ts_start, (long)ts_end, ilay);
              }
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
      break;
    }//end case 1

    case 2: {// thresholds
      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          switch(il){
            case 0: {
              string objname = Form("Threshold/Layer%d/Threshold_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], objname, run1, run2, (long)ts_start, (long)ts_end,layernum);
              break;
            }

            case 1: {
              string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",layernum);
              cout<<"\nAll data in "<<taskname[layernum]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], objname, run1, run2, (long)ts_start, (long)ts_end, layernum);
              break;
            }

            case 2: {
              for(int istave=0; istave<nStavesInLay[layernum]; istave++){
                Download(choice, ccdb, ccdbApi, myname, taskname[layernum], tasknamealt[layernum], Form("DeadPixel/Layer%d/Stave%d/DeadPixelHITMAP",layernum,istave), run1, run2, (long)ts_start, (long)ts_end, layernum);
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
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], objname, run1, run2, (long)ts_start, (long)ts_end,ilay);
              }
              break;
            }

            case 1: {
              for(int ilay=0; ilay<=2; ilay++){
                string objname = Form("DeadPixel/Layer%d/DeadPixel_Vs_Chip_and_Stave",ilay);
                cout<<"\nAll data in "<<taskname[ilay]+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
                Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], objname, run1, run2, (long)ts_start, (long)ts_end, ilay);
              }
              break;
            }

            case 2: {
              for(int ilay=0; ilay<=2; ilay++){
                for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                  Download(choice, ccdb, ccdbApi, myname, taskname[ilay], tasknamealt[ilay], Form("DeadPixel/Layer%d/Stave%d/DeadPixelHITMAP",ilay,istave), run1, run2, (long)ts_start, (long)ts_end, ilay);
                }
              }
              break;
            }

            default: break;
          }
        }//end loop on lists
      }//end if layernum==-1
    }// end of case 2
    break;
  }//end switch

  outputfile->Close();
  delete outputfile;

  return 1;
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
  GetListOfHisto(ccdb, myname, taskname, " ", objname, timestamps_selperiod, lnum, 0, isperstave, vector<int>());

  timestamps_selperiod.clear();
}


//
// Download data based on run numbers - available in metadata from 03/07/2019 21.49 (run 582 --> fake hit scan)
//
void DownloadRuns(auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string tasknamealternative, string objname, string run1, string run2, int lnum){

  //Extract all the time stamps and run numbers of the object
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  string objectlist2 = " ";
  if(tasknamealternative!=taskname){
    objectlist2 = ccdbApi.list(tasknamealternative + "/" + objname,false,"text/plain");
  }
  cout<<endl;
  cout<<endl;
  cout<<"Ready to get files from "<<taskname<<"/"<<objname<<endl;
  if(tasknamealternative!=taskname){
    cout<<"... And from (alternative path for backward compatibility): "<< tasknamealternative<<"/"<<objname<<endl;
  }
  stringstream ss(objectlist);
  stringstream ss2(objectlist2);
  string word;
  vector<string> alltimestamps, timestamps, runs;
  vector<string> alltimestampsALT, timestampsALT, runsALT;//for alternative path

  //filter string from alternative path (at the moment: only L1 from before run 300134)
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

  //filter normal path but correct timestamps with the one from alternative path (if needed)
  while(ss>>word){
    if(word=="Created:"){// take the one related to file creation
      ss>>word;
      alltimestamps.push_back(word);
    }
    if(word=="Run"){
      ss>>word;
      ss>>word;
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
      if(stoi(word)==stoi(run1)) break;
    }
  }

  //filter all runs to get the ones within run1 and run2 (get most recent for each run!!)
  vector <long int> timestamps_selperiod;
  vector <int> runs_selperiod;
  int counter=0;
  for(int irun=0; irun<(int)runs.size(); irun++){
    if(stoi(runs[irun])>=stoi(run1) && stoi(runs[irun])<=stoi(run2)){
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

  cout<<"\n"<<"Selected runs and corresponding timestamps: "<<endl;
  for(int i=0; i<(int)runs_selperiod.size();i++){
    cout<<"run"<<runs_selperiod[i]<<" - "<<timestamps_selperiod[i]<<endl;
  }

  bool isperstave = 0;
  if(objname.find("HITMAP")!=string::npos) isperstave = 1;
  GetListOfHisto(ccdb, myname, taskname, tasknamealternative, objname, timestamps_selperiod, lnum, 1, isperstave, runs_selperiod);

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
bool GetListOfHisto(auto* ccdb, string myname, string taskname, string tasknamealternative, string objname, vector<long int> timestamps, int lnum, bool isrunknown, bool isperstave, vector<int>runnumbers){

  //Getting root files from the database and write them to file
  cout<<"\n"<<"... Getting files from the database"<<endl;

  string stvnum = "0";
  if(isperstave){
    stvnum = objname.substr(objname.find("Stave")+5,2);
    if(stvnum.find("/")!=string::npos)
      stvnum = objname.substr(objname.find("Stave")+5,1);
  }

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
    string c = obj->ClassName();
    TH2 *h2s = 0x0;
    TH1 *h1s = 0x0;
    THnSparse *hSp = 0x0;

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
        histname = Form("h2_L%d_err%s_%ld", lnum, isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      else if(objname.find("Trigger")!=string::npos)
        histname = Form("h2_L%d_trg%s_%ld", lnum, isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);
      else
        histname = Form("h2_L%d%s%s_%ld", lnum, isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]);

      h2s = dynamic_cast<TH2*>(obj->Clone(histname.c_str()));
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

    delete h1s;
    delete h2s;
    delete hSp;
    //delete monitor;
    delete obj;
  }

  return 1;
}

//
// Download depending on the choice
//
bool Download(int choice, auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string tasknamealternative, string objname, string run1, string run2, long int ts_start, long int ts_end, int lnum){

  switch(choice){
    case 1: {
      DownloadRuns(ccdb, ccdbApi, myname, taskname, tasknamealternative, objname, run1, run2, lnum);//download with runs
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
