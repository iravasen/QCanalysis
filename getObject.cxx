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
#include <TList.h>
//#include "QualityControl/QcInfoLogger.h"

using namespace std;
using namespace o2::quality_control::repository;

//functions to download data
void DownloadTimestamps(auto *ccdb, string myname, string taskname, string objname, TList *list, time_t ts_start, time_t ts_end);
void DownloadRuns(auto *ccdb, string myname, string taskname, string objname, TList *list, string run1, string run2);
bool GetListOfHisto(auto* ccdb, string myname, string taskname, TList *list, string objname, int layernum, vector<long int> timestamps, bool isrunknown, bool isperstave, vector<int> runnumber);
bool Download(int choice, auto* ccdb, string myname, string taskname, string objname, TList *list, string run1, string run2, time_t ts_start, time_t ts_end);
string GetOptName(int opt);

const int nStavesInLay[3] = {12, 16, 20};


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

  //std::vector<string> tasks = ccdb->getListOfTasksWithPublications();
  std::vector<string> tasks = ccdb->getListing();

  //Get list of tasks
  std::vector<string>::iterator iTask;
  /*cout << "\n\nTasks with publications:" << endl;
  for(iTask = tasks.begin(); iTask != tasks.end(); iTask++) {
    string task = *iTask;
    cout << task << endl;
  }*/
  string taskname = "qc/ITS/ITSRawTask";
  string objname;
  //Choose what to download
  int opt;
  cout<<endl;
  cout<<endl;
  cout<<"Choose what to download:"<<endl;
  cout<<"1. Layer fake-hit rate maps (fake-hit rate runs)"<<endl;
  cout<<"2. Layer hitmaps (fake-hit rate runs)"<<endl;
  cout<<endl;
  cout<<"Enter the option: ";
  cin>>opt;
  if(opt<1 || opt>2){
    cout<<"Invalid option"<<endl;
    return 0;
  }

  //Choose the layer number
  int layernum;
  cout<<endl;
  cout<<endl;
  cout<<"Enter the layer number [put -1 for all IB layers]"<<endl;
  cin>>layernum;

  //Get list of objects inside the task
  /*std::vector<string> objects = ccdb->getPublishedObjectNames(taskname);

  std::vector<string>::iterator iObj;
  cout << "\n\nObjects for task " << taskname << endl;
  for(iObj = objects.begin(); iObj != objects.end(); iObj++) {
    string obj = *iObj;
    obj.erase(std::remove(obj.begin(), obj.end(), '\\'), obj.end());//remove backslash from object name
    obj.erase(0,1);//remove the first slash
    cout << obj << endl; //cout list of objects
  }

  string objname;
  cout << "Enter a object name to get its content [quit to exit]: ";
  cin >> objname;
  if(objname == "quit") return 0;
  */
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
      cout<<"\nAll data in "<<taskname+"/"+objname<<" between run"<<run1<<" and run"<<run2<<" are going to be downloaded."<<endl;
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

  TList *list = new TList();
  list->SetName("mylist");
  list->SetOwner();

  //Download depending on the option (opt)
  switch(opt){
    case 1: {
      if(layernum>=0){
        objname = Form("Occupancy/Layer%d/Layer%dChipStave",layernum,layernum);
        Download(choice, ccdb, myname, taskname, objname, list, run1, run2, ts_start, ts_end);
      }
      else if(layernum==-1){
        for(int ilay=0; ilay<=2; ilay++){
          objname = Form("Occupancy/Layer%d/Layer%dChipStave",ilay,ilay);
          Download(choice, ccdb, myname, taskname, objname, list, run1, run2, ts_start, ts_end);
        }
      }
    }//end case 1

    case 2: {

      if(layernum>=0){
        //std::vector<string> objects = ccdb->getPublishedObjectNames(taskname+Form("/Occupancy/Layer%d",layernum));
        //std::vector<string>::iterator iObj;
        //cout << "\n\nObjects for " << taskname+Form("/Occupancy/Layer%d",layernum) << endl;
        for(int istave=0; istave<nStavesInLay[layernum]; istave++){
          Download(choice, ccdb, myname, taskname, Form("/Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",layernum,istave,layernum,istave), list, run1, run2, ts_start, ts_end);
        }
      }
      else if(layernum==-1){
        for(int ilay=0; ilay<=2; ilay++){
          //std::vector<string> objects = ccdb->getPublishedObjectNames(taskname+Form("/Occupancy/Layer%d",ilay));
          //std::vector<string>::iterator iObj;
          //cout << "\n\nObjects for " << taskname+Form("/Occupancy/Layer%d",layernum) << endl;
          for(int istave=0; istave<nStavesInLay[ilay]; istave++){
            Download(choice, ccdb, myname, taskname, Form("/Occupancy/Layer%d/Stave%d/Layer%dStave%dHITMAP",ilay,istave,ilay,istave), list, run1, run2, ts_start, ts_end);
          }
        }//end loop on layers
      }
    }//end case 2
  }

  //other option may be included here (FUTURE)


  //Save list of objects into a file
  //Save files
  string objname_mod = objname;
  string taskname_mod = taskname;
  std::replace(objname_mod.begin(),objname_mod.end(),'/','-');
  std::replace(taskname_mod.begin(),taskname_mod.end(),'/','-');
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
  TFile *outputfile = new TFile(Form("Data/Output_%s_%s_from_%s%s_to_%s%s.root",layername.c_str(), optname.c_str(), suffix.c_str(),nums[0].c_str(), suffix.c_str(), nums[1].c_str()), "RECREATE");
  list->Write();
  outputfile->Close();

  ccdb->disconnect();

  return 1;
}//end main

//
// Download data based on timestamps
//
void DownloadTimestamps(auto* ccdb, string myname, string taskname, string objname, TList* list, time_t ts_start, time_t ts_end){

  //Extract all the time stamps of the object
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init("ccdb-test.cern.ch:8080");
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  //cout<<objectlist<<endl;
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
  GetListOfHisto(ccdb, myname, taskname, list, objname, timestamps_selperiod, 0, isperstave, vector<int>());
}


//
// Download data based on run numbers - available in metadata from 03/07/2019 21.49 (run 582 --> fake hit scan)
//
void DownloadRuns(auto* ccdb, string myname, string taskname, string objname, TList *list, string run1, string run2 ){

  //Extract all the time stamps and run numbers of the object
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init("ccdb-test.cern.ch:8080");
  string objectlist = ccdbApi.list(taskname + "/" + objname,false,"text/plain");
  cout<<"Ready to get files from "<<taskname<<"/"<<objname<<endl;
  //cout<<objectlist<<endl;
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
      //cout<<runs[runs.size()-1]<<", "<<timestamps[timestamps.size()-1]<<endl;
    }

    /*if(timestamps.size() - runs.size() > 1){//
      timestamps.erase(timestamps.end() - 2, timestamps.end() - 1);//remove second-to-last and last timestamp
      break;
    }*/
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
  GetListOfHisto(ccdb, myname, taskname, list, objname, timestamps_selperiod, 1, isperstave, runs_selperiod);

}




//
//Get list of histogram inside an object
//
bool GetListOfHisto(auto* ccdb, string myname, string taskname, TList *list, string objname, vector<long int> timestamps, bool isrunknown, bool isperstave, vector<int>runnumbers){
  //Getting root files from the database and save them in 1 file
  //TList *list = new TList();
  //list->SetOwner();
  //list->SetName("mylist");
  TH2 *h2s;// = new TH2();
  TH1 *h1s;// = new TH1();
  cout<<"\n"<<"... Getting files from the database"<<endl;
  o2::quality_control::core::MonitorObject *monitor;

  //Get layer number from object name
  string lnum = objname.substr(objname.find("Layer")+5,1);
  string stvnum = "0";
  if(isperstave){
    stvnum = objname.substr(objname.find("Stave")+5,2);
    if(stvnum.find("/")!=string::npos)
      stvnum = objname.substr(objname.find("Stave")+5,1);
  }


  for(int i=0; i<(int)timestamps.size();i++){
    monitor = ccdb->retrieve(taskname, objname, timestamps[i]);
    if (monitor == nullptr) {
      cerr << myname << ": failed to get MonitorObject for timestamp: " << timestamps[i]<< endl;
      return 0;
    }
    //std::unique_ptr<TObject> obj(monitor->getObject());
    TObject *obj = monitor->getObject();
    monitor->setIsOwner(false);
    //obj->SaveAs(Form("%ld_file.root", timestamps_selperiod[i]));
    //TH2S *h = new TH2S();
    const char* c = obj->ClassName();
    //const char type[] = "TH1S";
    //h1s->Reset();
    //h2s->Reset();
    if(strstr(c,"TH1")!=nullptr){
      h1s = (TH1*)obj->Clone(Form("h1_L%s%s%s_%ld", lnum.c_str(), isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]));
      list->Add(h1s);
    }
    //cout<<"BEFORE IF"<<endl;
    if(strstr(c,"TH2")!=nullptr){
      h2s = (TH2*)obj->Clone(Form("h2_L%s%s%s_%ld", lnum.c_str(), isperstave ? Form("_Stv%s",stvnum.c_str()) : "", isrunknown ? Form("_run%d",runnumbers[i]) : "", timestamps[i]));
      list->Add(h2s);
      //cout<<"INSIDE IF"<<endl;
    }
    //cout<<"AFTER IF"<<endl;
    //h->Delete();
  }
  //cout<<"END OF LOOP"<<endl;
  return 1;
  //return list;
}

//
// Download depending on the choice
//
bool Download(int choice, auto* ccdb, string myname, string taskname, string objname, TList *list, string run1, string run2, time_t ts_start, time_t ts_end){

  switch(choice){
    case 1: {
      DownloadRuns(ccdb, myname, taskname, objname, list, run1, run2);
      break;//download runs data
    }

    case 2: {
      DownloadTimestamps(ccdb, myname, taskname, objname, list, ts_start, ts_end); //download
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
    case 1: return "FHRMAPS";
    case 2: return "HITMAPS";
    default: return "0";
  }
}
