#include <string>
#include <iostream>
#include <vector>
#include <array>
#include <TH2.h>
#include <TH3.h>
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
#include <TLatex.h>
#include <THnSparse.h>

#include "inc/constants.h"

using namespace std;

//Functions
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap);
int GetNchipsActive(THnSparse *hmap, int maxchip);
int GetNrunsWOhits(TH2 *hFhrStv);
void SetStyle(TH1 *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath_hit, const int nChips, string skipruns, bool isIB);


void MaskNoisyPixelsInRuns(){
  string fpath;

  int nchips=9;
  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
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

  DoAnalysis(fpath, nchips, skipruns, isIB);
}

//
// Analysis
//
void DoAnalysis(string filepath_hit, const int nChips, string skipruns, bool isIB){

  gStyle->SetOptStat(0000);

  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  std::vector<THnSparse*> hmaps;
  std::vector<string> timestamps, runnumbers, stavenums, laynums;
  vector<int> posrefrun;
  int nLayers=1, nRuns=1;
  vector<string> runlabel;

  //Read the file and the list of plots with entries (hitmaps!)
  TFile *infile=new TFile(filepath_hit.c_str());
  TList *list = (TList*)infile->GetListOfKeys();
  TIter next(list);
  TKey *key;
  TObject *obj;
  while((key=(TKey*)next())){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1")
         && (!obj->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("Stv")==string::npos) continue;
    THnSparse *h2 = (THnSparse*)obj;
    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    string stvnum = objname.substr(objname.find("Stv")+3,2);
    if(stvnum.find("_")!=string::npos){
      stvnum = objname.substr(objname.find("Stv")+3,1);
    }
    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    hmaps.push_back(h2);

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

    if((int)stavenums.size()>1 && stvnum!=stavenums[stavenums.size()-2])
      continue;
    runlabel.push_back(runnum);

  }

  //Read file with fhr maps for each layer --> ONLY TO EXTRACT NUMBER OF TRIGGERS!!
  TList *list2 = (TList*)infile->GetListOfKeys();
  TIter next2(list2);
  TH2 *h2_2[nRuns];
  TKey *key2;
  TObject *obj2;
  int cntrun=0;
  while((key2=((TKey*)next2()))){
    obj2 = key2->ReadObj();
    if ((strcmp(obj2->IsA()->GetName(),"TProfile")!=0)
         && (!obj2->InheritsFrom("TH2"))
	       && (!obj2->InheritsFrom("TH1")
         && (!obj2->InheritsFrom("THnSparse")))
       ) {
            cout<<"<W> Object "<<obj2->GetName()<<" is not 1D, 2D or sparse histogram : will not be converted"<<endl;
       }
    string objname2 = (string)obj2->GetName();
    if(objname2.find("Stv")!=string::npos) break;
    //TH2 *h2temp = (TH2*)obj2;
    //if(!h2temp->GetEntries()) continue;
    string runnum =  objname2.find("run")==string::npos ? "norun":objname2.substr(objname2.find("run")+3, 6);
    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    h2_2[cntrun] = (TH2*)obj2;
    cntrun++;
    if(cntrun==nRuns+1) break;
    //if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj2->GetName()<<endl;
  }

  //Calculate number of triggers for each run
  cout<<endl;
  cout<<"... Extract the number of triggers (events) from for each run"<<endl;
  vector<double> ntrig;
  for(int ir=0; ir<nRuns; ir++){
    double fhr_run = h2_2[ir]->GetBinContent(1,1);
    int stavefound = 0;
    int chipfound = 0;
    if(fhr_run<1e-15){
      for(int ibinx=1; ibinx<=h2_2[ir]->GetNbinsX(); ibinx++){
        for(int ibiny=h2_2[ir]->GetNbinsY(); ibiny>=1; ibiny--){
          fhr_run = h2_2[ir]->GetBinContent(ibinx,ibiny);
          if(fhr_run>1e-15) {
            stavefound = ibiny-1;
            chipfound = ibinx-1;
            break;
          }
        }
        if(fhr_run>1e-15) break;
      }
    }
    if(fhr_run<1e-15){
      fhr_run = -1.;
      ntrig.push_back(-1.);
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers (ignored in the calculation of the average fhr)"<<endl;
      continue;
      //cout<<"INVALID FHR... setting it to -1"<<endl;
    }


    hmaps[stavefound*nRuns+ir] -> GetAxis(0) ->SetRange(1+1024*chipfound, 1024+1024*chipfound);
    TH2F *hprojsparse = (TH2F*)hmaps[stavefound*nRuns+ir]->Projection(1,0);
    hmaps[stavefound*nRuns+ir] -> GetAxis(0) ->SetRange(1, 9216);//reset the range
    double hits_chip = hprojsparse->Integral(1,1024,1,512);

    if(hits_chip/(512.*1024.*fhr_run) < 1e-15){//to avoid bugs due to bad runs
       ntrig.push_back(-1.);
       cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers (ignored in the calculation of the average fhr)"<<endl;
    }
    else{
      ntrig.push_back(hits_chip/(512.*1024.*fhr_run));
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers"<<endl;
    }
    delete hprojsparse;
  }

  //Start masking hottest pixels for each stave in each run, Fill also the histo with the hot pixel maps for each stave
  cout<<endl;
  cout<<"... Analysing FHR with (hot) pixel masking (Making also hot pixel map)"<<endl;
  vector<array<float,nMasked+1>> fhrall;
  TH2F *hHotMap[nLayers][20];
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int istave=0; istave<20; istave++)
      hHotMap[ilay][istave] = new TH2F(Form("hHotMap_L%s_Stv%d",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(), istave), "; ; ", 9216,0,9216, 512,0,512);
  int irun = nRuns-1;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    int nchipsactive = GetNchipsActive(hmaps[ihist],nChips);
    if(nchipsactive<nChips)
      cout<<"Layer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<<" Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nchipsactive<<endl;
    fhrall.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : stoi(laynums[ihist])][stoi(stavenums[ihist])]));
    irun--;
    if(ihist>0){
      if(stavenums[ihist-1]!=stavenums[ihist]){
        irun=nRuns-1;
      }
    }
  }

  //special binning
  double binstart = 0.4;
  double binsmasked[nMasked+2];
  binsmasked[0] = 0.1;
  binsmasked[1] = 0.5;
  for(int i=2; i<=nMasked+1; i++){
    binsmasked[i] = binsmasked[i-1]+1.;
  }

  TH2F *hFhrStv[nLayers][100];
  for(int ilay=0; ilay<nLayers; ilay++)
    for(int is=0; is<nStavesInLay[ilay]; is++)
      hFhrStv[ilay][is] = new TH2F(Form("h2FhrStv_%s_%d", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),is), Form("Layer-%s - Stave-%d; # Hot Pixel Masked;Run", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),is),nMasked+1, binsmasked, nRuns, 0.5, nRuns+0.5);

  //Fill histogram
  int ilayer = nLayers-1;
  int istave = nStavesInLay[ilayer]-1;
  irun=0;
  for(int i=0; i<(int)fhrall.size(); i++){
    for(int ifhr=0; ifhr<(int)fhrall[i].size(); ifhr++){
      hFhrStv[ilayer][istave]->SetBinContent(ifhr+1, irun+1, fhrall[i][ifhr]);
    }
    irun++;
    if(i<(int)fhrall.size()-1){// in case ref run is the first into the list of runs
      if(stavenums[fhrall.size()-i-2]!=stavenums[fhrall.size()-i-1]){
        istave--;
        irun=0;
      }
      if(laynums[fhrall.size()-i-2]!=laynums[fhrall.size()-i-1]){
        ilayer--;
        istave = nStavesInLay[ilayer]-1;
      }
    }
  }

  //Make FHR (averaged on all runs) vs #masked pix for all staves in a layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay));
    cnv.cd();
    cnv.SetLogy();
    cnv.SetLogx();
    cnv.SetTickx();
    cnv.SetTicky();
    cnv.SetMargin(0.0988,0.1,0.1,0.0993);
    TH1F *hframe = cnv.DrawFrame(0.1,7e-15,3*(nMasked),1e-3, Form("Layer %s - Average FHR %s; # Hot Pixels masked ; FHR (/event/pixel)",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
    hframe->SetBins(nMasked+1,binsmasked);
    //legend
    TLegend leg(0.904, 0.127,0.997,0.898);
    for(int is=0; is<nStavesInLay[ilay];is++){
      TH1F *proj = (TH1F*)hFhrStv[ilay][is]->ProjectionX(Form("proj_%d%d",ilay,is));
      //proj->SetTitle(Form("Layer %s - Average FHR %s; # Hot Pixels masked ; FHR (/event/pixel)",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
      int runswohits = GetNrunsWOhits(hFhrStv[ilay][is]);
      //cout<<runswohits<<endl;
      proj->Scale(1./(nRuns-runswohits)); //Divide by the number of runs minus the ones without hits
      SetStyle(proj, col[is<nStavesInLay[ilay]/2 ? is : is-nStavesInLay[ilay]/2],is<nStavesInLay[ilay]/2 ? 24:26);

      proj->Draw("PL same");

      leg.AddEntry(proj, Form("Stv%d",is),"pl");
    }
    leg.Draw("same");
    cnv.SaveAs(Form("../Plots/Layer%s_FHRpixmask_%s.pdf", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
    cnv.SaveAs(Form("../Plots/Layer%s_FHRpixmask_%s.root", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
  }

  //Draw hot pixel maps for each layer
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay),800,1200);
    cnv.SetTopMargin(0.4);
    cnv.Divide(1,nStavesInLay[ilay],0,0);
    for(int istave=0; istave<nStavesInLay[ilay]; istave++){
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
      //hHotMap[ilay][istave]->Draw("COLZ");
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
    lat.DrawLatex(0.01,0.98,Form("L%s",laynums[ilay*nRuns*nStavesInLay[ilay]].c_str()));

    cnv.SaveAs(Form("../Plots/Layer%s_Hotpixmap_%s.pdf", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
    cnv.SaveAs(Form("../Plots/Layer%s_Hotpixmap_%s.root", laynums[ilay*nRuns*nStavesInLay[ilay]].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
  }

}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap){

  array<float,nMasked+1> fhrstave;

  //clone
  THnSparse *hmapclone = (THnSparse*)hmap->Clone(Form("%s_clone",hmap->GetName()));
  cout<<hmap->GetName()<<" --> "<<hmapclone->GetNbins()<<" noisy pixels"<<endl;

  for(int iter=0; iter<nMasked+1; iter++){

    TH1F *hproj = (TH1F*)hmapclone->Projection(1);
    long int totalhits = hproj->Integral();
    float fhr = nchips==0 ? 0. : (float)totalhits / (512.*1024.*nchips*ntrig);


    //cout<<fhr<<endl;
    if(ntrig<0) fhr=0.;
    fhrstave[iter] = fhr;

    int coord[2];
    double max = -1.;
    int x=0,y=0;
    long int binwithmax = 0;
    for(int ibin=0; ibin<hmapclone->GetNbins(); ibin++){
      double bincontent = hmapclone->GetBinContent(ibin, coord);
      if(bincontent>max){
        max=bincontent;
        binwithmax = ibin;
        x=coord[0];
        y=coord[1];
      }
    }

    if(nchips) hmapclone->SetBinContent(binwithmax,0.);
    if(totalhits!=0) hhotmap->SetBinContent(x-1, y-1, 1);// to avoid a marker in 0,0 for empty histos
    delete hproj;
  }

  delete hmapclone;

  return fhrstave;

}

//
//Function to return the number of active (i.e. enabled) chips in a stave
//
int GetNchipsActive(THnSparse *hmap, int maxchip){
  int ix1=1;
  int ix2=1024;
  int activechips = maxchip;

  while(ix2<=hmap->GetAxis(0)->GetNbins()){
    hmap->GetAxis(0)->SetRange(ix1,ix2);
    TH2F *hproj = (TH2F*)hmap->Projection(1,0);
    if((hproj->Integral(1,1024,1,512)<1e-15))
      activechips--;
    ix1=ix2+1;
    ix2+=1024;
    delete hproj;
  }
  hmap->GetAxis(0)->SetRange(1,9216);//reset range
  return activechips;
}

//
//Function to return the number of runs in which a stave has no hits
//
int GetNrunsWOhits(TH2 *hFhrStv){
  int runswohits = 0;
  for(int iy=1; iy<=hFhrStv->GetNbinsY(); iy++){
    if(hFhrStv->GetBinContent(1,iy)<1e-15)
      runswohits++;
  }
  return runswohits;
}

//
//Set Style
//
void SetStyle(TH1 *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}
