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

  //Draw hot pixel maps for each layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay),800,1200);
    cnv.SetTopMargin(0.4);
    cnv.Divide(1,nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])],0,0);
    for(int istave=0; istave<nStavesInLay[nLayers>1 ? ilay:stoi(laynums[0])]; istave++){
      TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%d",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(), istave), "; ; ", 9216,0.5,9216.5,512,0.5,512.5);
      int cnt = 0;

      for(int ihist=0; ihist<(int)hmaps.size(); ihist++){ //start from the bottom in order to start with the oldest run
        if(nLayers>1){
          if(stoi(laynums[ihist])==ilay && stoi(stavenums[ihist]) != istave) continue;
          if(stoi(laynums[ihist])!=ilay) continue;
        }
        else if(stoi(stavenums[ihist]) != istave) {
          continue;
        }
        if(!cnt) {hHotMap = (TH2F*)hmaps[ihist]->Projection(1,0);}
        else {
          TH2F *htemp = (TH2F*)hmaps[ihist]->Projection(1,0);
          hHotMap->Add(htemp);
          delete htemp;
        }

        if(ihist>0 && stavenums[ihist+1]!=stavenums[ihist]) break;
        cnt++;
      }
      cnv.cd(istave+1);
      cnv.GetPad(istave+1)->SetTickx();
      cnv.GetPad(istave+1)->SetTicky();
      cnv.GetPad(istave+1)->SetRightMargin(0.01);
      if(!istave) cnv.GetPad(istave+1)->SetTopMargin(0.1);

      hHotMap->SetTitle(" ");
      hHotMap->SetMarkerStyle(20);
      hHotMap->SetMarkerSize(0.6);
      hHotMap->SetMarkerColor(kRed);
      hHotMap->SetLineColor(kRed);

      hHotMap->GetXaxis()->SetRangeUser(0.,9216.);
      hHotMap->GetYaxis()->SetRangeUser(0.,512.);

      hHotMap->GetXaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetLabelSize(0.13);
      hHotMap->GetXaxis()->SetLabelSize(0.13);
      if(istave>0){
        hHotMap->GetXaxis()->SetLabelOffset(999);
        hHotMap->GetXaxis()->SetTickLength(0.05);
        hHotMap->GetXaxis()->SetNdivisions(530);
      }
      else{
        hHotMap->GetXaxis()->SetLabelOffset(0.003);
        hHotMap->GetXaxis()->SetNdivisions(530);
        hHotMap->GetXaxis()->SetTickLength(0.05);
      }

      hHotMap->DrawCopy("P X+");

      TLatex lat;
      lat.SetTextAngle(90);
      lat.SetNDC();
      lat.SetTextSize(0.15);
      lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));

      delete hHotMap;
    }//end loop on staves
    cnv.cd();
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.03);
    lat.DrawLatex(0.01,0.98,Form("L%s",nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str()));

    cnv.SaveAs(Form("../Plots/Layer%s_Deadpixmap_%s.pdf", nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    //cnv.SaveAs(Form("../Plots/Layer%s_Deadpixmap_%s.root", nLayers>1 ? to_string(ilay).c_str():laynums[0].c_str(),filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

}
