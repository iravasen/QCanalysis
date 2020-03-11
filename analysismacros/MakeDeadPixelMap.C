#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TText.h>
#include <TSystem.h>
#include <TKey.h>
#include <THnSparse.h>
#include <TLatex.h>

#include "inc/constants.h"

using namespace std;

//Functions
std::array<long int,5> CompareTwoRuns(THnSparse *href, THnSparse *h2);
void SetStyle(TGraphErrors *ge, Color_t col);
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns);

void MakeDeadPixelMap(){
  string fpath;
  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*THRMAPS_DEADPIXMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;

  bool isIB;
  if(fpath.find("IB")!=string::npos){
    isIB = kTRUE;
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    if(stoi(layernum)>=0 && stoi(layernum)<=2) nchips = 9;
    else if (stoi(layernum)==3 && stoi(layernum)==4) nchips = 54*2;
    else nchips = 98*2;
    if(nchips==9) isIB=kTRUE;
    else isIB=kFALSE;
  }

  //Choose whether to skip runs
  string skipans, skipruns;
  cout<<endl;
  cout<<"Would you like to skip some run(s)? [y/n] ";
  cin>>skipans;
  if(skipans=="y" || skipans=="Y"){
    cout<<endl;
    cout<<"Specify run number(s) separated by comma (no white spaces!):";
    cin>>skipruns;
    cout<<endl;
  }
  else
    skipruns=" ";


  DoAnalysis(fpath, nchips, isIB, skipruns);
}

//
// Analysis
//
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns){

  gStyle->SetOptStat(0000);

  std::vector<THnSparse*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  vector<int> posrefrun;
  int nLayers=1, nRuns=1;

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TObject *obj;
  TKey *key;
  THnSparse *hsparse;
  while((key=((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D, sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    hsparse = (THnSparse*)obj;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    hmaps.push_back(hsparse);

    if((int)laynums.size()>1 && laynum!=laynums[laynums.size()-1]){
      nLayers++;
    }

    if((int)stavenums.size()>1 && stvnum==stavenums[stavenums.size()-1])
      nRuns++;
    else nRuns=1;

    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);
    stavenums.push_back(stvnum);
  }

  //Compare all the runs (non-empty ones) with the reference run chosen by the user
  vector<string> runlabel;
  int istave = posrefrun.size()-1;
  int ilayer = nLayers-1;
  TH2F *hHotMap[nLayers][20];
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int istave=0; istave<20; istave++)
      hHotMap[ilay][istave] = new TH2F(Form("hHotMap_L%s_Stv%d",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(), istave), "; ; ", 9216,0.5,9216.5,512,0.5,512.5);

  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    hHotMap[ilayer][stoi(stavenums[ihist])] = (TH2F*)hmaps[ihist]->Projection(1,0);
    hHotMap[ilayer][stoi(stavenums[ihist])]->SetName(Form("hHotMap_L%s_Stv%s",laynums[ihist].c_str(), stavenums[ihist].c_str()));
    if(ihist>0){
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
    }
  }//end loop on histograms


  //Draw hot pixel maps for each layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay),800,1200);
    cnv.SetTopMargin(0.4);
    cnv.Divide(1,nStavesInLay[ilay],0,0);
    for(int istave=0; istave<nStavesInLay[ilay]; istave++){
      hHotMap[ilay][istave]->SetTitle(" ");
      hHotMap[ilay][istave]->SetMarkerStyle(20);
      hHotMap[ilay][istave]->SetMarkerSize(0.6);
      hHotMap[ilay][istave]->SetMarkerColor(kRed);
      hHotMap[ilay][istave]->SetLineColor(kRed);

      cnv.cd(istave+1);
      cnv.GetPad(istave+1)->SetTickx();
      cnv.GetPad(istave+1)->SetTicky();
      cnv.GetPad(istave+1)->SetRightMargin(0.01);
      if(!istave) cnv.GetPad(istave+1)->SetTopMargin(0.1);

      hHotMap[ilay][istave]->Draw("P X+");
      hHotMap[ilay][istave]->GetXaxis()->SetRangeUser(0.,9216.);
      hHotMap[ilay][istave]->GetYaxis()->SetRangeUser(0.,512.);

      hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.005);
      hHotMap[ilay][istave]->GetYaxis()->SetTickLength(0.005);
      hHotMap[ilay][istave]->GetYaxis()->SetLabelSize(0.13);
      hHotMap[ilay][istave]->GetXaxis()->SetLabelSize(0.13);
      if(istave>0){
        hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(999);
        hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
        hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
      }
      else{
        hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(0.003);
        hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
        hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
      }

      TLatex lat;
      lat.SetTextAngle(90);
      lat.SetNDC();
      lat.SetTextSize(0.15);
      lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
    }
    cnv.cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.01,0.98,Form("L%s",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str()));

    cnv.SaveAs(Form("../Plots/Layer%s_Deadpixmap_%s.pdf", nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    cnv.SaveAs(Form("../Plots/Layer%s_Deadpixmap_%s.root", nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

}
