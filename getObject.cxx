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
#include <TFile.h>

using namespace std;
using namespace o2::quality_control::repository;
using namespace o2::quality_control::core;

//functions to download data
void DownloadTimestamps(auto *ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, time_t ts_start, time_t ts_end);
void DownloadRuns(auto *ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2);
bool GetListOfHisto(auto* ccdb, string myname, string taskname, string objname, vector<long int> timestamps, bool isrunknown, bool isperstave, vector<int> runnumber);
bool Download(int choice, auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2, time_t ts_start, time_t ts_end);
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

  //taskname
  string taskname = "qc/ITS/ITSRawTask";

  //Choose what to download
  int opt;
  cout<<endl;
  cout<<endl;
  cout<<"Choose what to download:"<<endl;
  cout<<"1. Fake-hit scan data"<<endl;
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
  }

  //Decide how many different elements are needed
  int nListElements = 2;
  switch(opt){
    case 1: nListElements = 2; break;
    default: nListElements = 2;
  }

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
  outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),nums[0].c_str(), suffix.c_str(), nums[1].c_str()), "RECREATE");
  outputfile->cd();

  //Download depending on the option (opt)
  switch(opt){
    case 1: {
      if(layernum>=0){
        for(int il=0; il<nListElements; il++){//loop on lists
          if(!il){
            string objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
            cout<<"\nAll data in "<<taskname+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
            Download(choice, ccdb, ccdbApi, myname, taskname, objname, run1, run2, ts_start, ts_end);
          }
          else{
            for(int istave=0; istave<nStavesInLay[layernum]; istave++){
              Download(choice, ccdb, ccdbApi, myname, taskname, Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), run1, run2, ts_start, ts_end);
            }
          }

        }//end loop on lists
      }//end if layernum>=0

      else if(layernum==-1){
        for(int il=0; il<nListElements; il++){//loop on lists
          if(!il){
            for(int ilay=0; ilay<=2; ilay++){
              string objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
              cout<<"\nAll data in "<<taskname+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
              Download(choice, ccdb, ccdbApi, myname, taskname, objname, run1, run2, ts_start, ts_end);
            }
          }
          else{
            for(int ilay=0; ilay<=2; ilay++){
              for(int istave=0; istave<nStavesInLay[ilay]; istave++){
                Download(choice, ccdb, ccdbApi, myname, taskname, Form("Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), run1, run2, ts_start, ts_end);
              }
            }
          }
        }//end loop on lists
      }//end if layernum==-1
    }//end case 1

  }//end switch

  outputfile->Close();
  delete outputfile;

  ccdb->disconnect();

  return 1;
}//end main

//
// Download data based on timestamps
//
void DownloadTimestamps(auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, time_t ts_start, time_t ts_end){

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
    if(stol(timestamps[its])>=(long)ts_start && stol(timestamps[its])<=(long)ts_end){
      counter++;
      if(counter==1) {
        if(its==0)
          timestamps_selperiod.push_back(stol(timestamps[its]));//the first (most recent) is taken by definition if it is really the first in the list
        else if(stol(timestamps[its-1])-stol(timestamps[its])>65000) timestamps_selperiod.push_back(stol(timestamps[its]));//if not the real first, check if it is the most recent (in this way the user can insert whatever date interval)
      }
      if(stol(timestamps[its])-stol(timestamps[its+1])>65000) timestamps_selperiod.push_back(stol(timestamps[its+1]));
    }

    else if(stol(timestamps[its])>(long)ts_end) continue;
    else break;
  }
  cout<<"\n"<<"Selected timestamps: "<<endl;
  for(int i=0; i<(int)timestamps_selperiod.size();i++){
    cout<<timestamps_selperiod[i]<<endl;
  }

  bool isperstave = 0;
  if(objname.find("HITMAP")!=string::npos) isperstave = 1;
  GetListOfHisto(ccdb, myname, taskname, objname, timestamps_selperiod, 0, isperstave, vector<int>());

  timestamps_selperiod.clear();
}


//
// Download data based on run numbers - available in metadata from 03/07/2019 21.49 (run 582 --> fake hit scan)
//
void DownloadRuns(auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2 ){

  //Extract all the time stamps and run numbers of the object
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  cout<<endl;
  cout<<endl;
  cout<<"Ready to get files from "<<taskname<<"/"<<objname<<endl;

  stringstream ss(objectlist);
  string word;
  vector<string> alltimestamps, timestamps, runs;

  while(ss>>word){
    if(word=="Created:"){// take the one related to file creation
      ss>>word;
      alltimestamps.push_back(word);
    }
    if(word=="Run"){
      ss>>word;
      ss>>word;
      runs.push_back(word);
      timestamps.push_back(alltimestamps[alltimestamps.size()-1]);//this keep only the timestamps connected to a run number
      if(stoi(word)==stoi(run1)) break;
    }
  }

  //filter all runs to get the ones within run1 and run2 (get most recent for each run!!)
  vector <long int> timestamps_selperiod;
  vector <int> runs_selperiod;
  int counter=0;
  for(int irun=0; irun<(int)runs.size()-1; irun++){
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
  GetListOfHisto(ccdb, myname, taskname, objname, timestamps_selperiod, 1, isperstave, runs_selperiod);

  timestamps_selperiod.clear();
  runs_selperiod.clear();
  timestamps.clear();
  runs.clear();
  alltimestamps.clear();

}




//
//Get list of histogram inside an object
//
bool GetListOfHisto(auto* ccdb, string myname, string taskname, string objname, vector<long int> timestamps, bool isrunknown, bool isperstave, vector<int>runnumbers){

  //Getting root files from the database and write them to file
  cout<<"\n"<<"... Getting files from the database"<<endl;

  //Get layer number from object name
  string lnum = objname.substr(objname.find("Layer")+5,1);
  string stvnum = "0";
  if(isperstave){
    stvnum = objname.substr(objname.find("Stave")+5,2);
    if(stvnum.find("/")!=string::npos)
      stvnum = objname.substr(objname.find("Stave")+5,1);
  }

  for(int i=0; i<(int)timestamps.size();i++){

    MonitorObject *monitor = ccdb->retrieve(taskname, objname, timestamps[i]);

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

    //if(strstr(c,"TH1")!=nullptr){
    if(c.find("TH1")!=string::npos){
      h1s = dynamic_cast<TH1*>(obj->Clone(Form("h1_L%s%s%s_%ld", lnum.c_str(), isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i])));
      outputfile->cd();
      h1s->Write();
    }
    //if(strstr(c,"TH2")!=nullptr){
    if(c.find("TH2")!=string::npos){
      h2s = dynamic_cast<TH2*>(obj->Clone(Form("h2_L%s%s%s_%ld", lnum.c_str(), isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i])));
      outputfile->cd();
      h2s->Write();
    }

    delete h1s;
    delete h2s;
    delete monitor;
    delete obj;
  }

  return 1;
}

//
// Download depending on the choice
//
bool Download(int choice, auto* ccdb, o2::ccdb::CcdbApi ccdbApi, string myname, string taskname, string objname, string run1, string run2, time_t ts_start, time_t ts_end){

  switch(choice){
    case 1: {
      DownloadRuns(ccdb, ccdbApi, myname, taskname, objname, run1, run2);
      break;//download runs data
    }

    case 2: {
      DownloadTimestamps(ccdb, ccdbApi, myname, taskname, objname, ts_start, ts_end); //download
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
