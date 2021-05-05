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
#include <TStopwatch.h>
#include "inc/constants.h"

using namespace std;

//Functions
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap, bool HalfStave, bool IB,  bool isHotPixelMapDrawn);
int GetNchipsActive(THnSparse *hmap, int maxchip, int MaxRange, bool HalfStave, bool IB);
int GetNrunsWOhits(TH2 *hFhrStv);
void SetStyle(TH1 *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath_hit, string skipruns, int IBorOB, bool isHotPixelMapDrawn);
int nStavesInLay[7]= {0};

void MaskNoisyPixelsInRuns(){
  string fpath;

  cout<<"\n\nAvailable file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*FHRMAPS_HITMAPS* -Art | tail -n 500");
  cout<<"\nCopy file name: ";
  cin>>fpath;
  cout<<endl;
  int IBorOB;

  //IBorOB = 0 if I want to check all IB layers                                              
  //IBorOB = 1 if I want to check all OB layers                                             
  //IBorOB = 2 if I want to check all IB + OB layers or if I want to check a single layer         

  if(fpath.find("IB")!=string::npos){
    IBorOB = 0;
    nStavesInLay[0] = nStavesInLayAll[0];
    nStavesInLay[1] = nStavesInLayAll[1];
    nStavesInLay[2] = nStavesInLayAll[2];
  }
  else if (fpath.find("OB")!=string::npos){
    IBorOB = 1;
    nStavesInLay[0] = nStavesInLayAll[3];
    nStavesInLay[1] = nStavesInLayAll[4];
    nStavesInLay[2] = nStavesInLayAll[5];
    nStavesInLay[3] = nStavesInLayAll[6];
  }
  else if (fpath.find("all")!=string::npos){
    IBorOB = 2;
    for (Int_t i =0; i < 7; i ++){
      nStavesInLay[i] = nStavesInLayAll[i];
    }
  }
  else{
    string layernum = fpath.substr(fpath.find("Layer")+5, 1);
    IBorOB = 2;
    nStavesInLay[0] = nStavesInLayAll[stoi(layernum)];
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

  string drawpixelmap;
  Bool_t isHotPixelMapDrawn =0;
  cout << "Would you like to save the hot pixel map? [y/n] ";
  cin >> drawpixelmap;
  if (drawpixelmap=="y" || drawpixelmap=="Y"){
    cout << endl;
    isHotPixelMapDrawn =1;
  }

  DoAnalysis(fpath, skipruns, IBorOB, isHotPixelMapDrawn);
}

//
// Analysis
//
void DoAnalysis(string filepath_hit, string skipruns, int IBorOB, bool isHotPixelMapDrawn){

  //TStopwatch t;
  //TStopwatch t1;
  //t.Start();
  //t1.Start();
  gStyle->SetOptStat(0000);

  Int_t col[] = {TColor::GetColor("#ff3300"), 807, 799, 829, 409, 839, TColor::GetColor("#00ceff"), 867, 859, TColor::GetColor("#0033a1"), TColor::GetColor("#708090"), 877};

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
  Int_t MaxRange  = 0; 
  Int_t MaxRangeY = 0;
  Int_t NChipsPerHIC   = 7;  
  Int_t ChipRowsPerHIC = 2;
  Int_t HICsPerHS = 0; //Number of HICs in each half stave
  Int_t HalfStaveFound = 0;
  vector<double> ntrig;

  for(int ir=0; ir<nRuns; ir++){
    double fhr_run = h2_2[ir]->GetBinContent(1,1);
    int stavefound = 0;
    int chipfound = 0;
    int chipfoundeff = 0;
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

    MaxRange =  hmaps[stavefound*nRuns+ir] -> GetAxis(0)->GetXmax();
    MaxRangeY =  hmaps[stavefound*nRuns+ir] -> GetAxis(1)->GetXmax();

    chipfoundeff = chipfound;
    if (IBorOB == 0 || (IBorOB == 2 && stoi(laynums[0]) < 3)) {
      NChipsPerHIC=1;   //IB layers do not have HICs
      ChipRowsPerHIC=1; //IB layers do not have HICs
    }
    else {
      if (IBorOB ==1) HICsPerHS=4;
      else if (IBorOB == 2){
	if (stoi(laynums[0])<5) HICsPerHS = 4;
	else HICsPerHS = 7;
      }
      if (chipfound < HICsPerHS) HalfStaveFound=0;
      else {
	HalfStaveFound =1; 
	chipfoundeff = chipfound-HICsPerHS;
      }
    }

    hmaps[stavefound*nRuns+ir] -> GetAxis(0) ->SetRange(1+1024*chipfoundeff*NChipsPerHIC, NChipsPerHIC*1024*(chipfoundeff+1));
    hmaps[stavefound*nRuns+ir] -> GetAxis(1) ->SetRange(1+512*HalfStaveFound*ChipRowsPerHIC, ChipRowsPerHIC*512*(HalfStaveFound+1));

    TH2F *hprojsparse = (TH2F*)hmaps[stavefound*nRuns+ir]->Projection(1,0);
    hmaps[stavefound*nRuns+ir] -> GetAxis(0) ->SetRange(1, MaxRange);//reset the range
    hmaps[stavefound*nRuns+ir] -> GetAxis(1) ->SetRange(1, MaxRangeY);//reset the range

    double hits_chip = hprojsparse->Integral(1,1024*NChipsPerHIC,1,512*ChipRowsPerHIC);
    delete hprojsparse;

    if(hits_chip/(ChipRowsPerHIC*512.*NChipsPerHIC*1024.*fhr_run) < 1e-15){//to avoid bugs due to bad runs
      ntrig.push_back(-1.);
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers (ignored in the calculation of the average fhr)"<<endl;
    }
    else{
      ntrig.push_back(hits_chip/(ChipRowsPerHIC*512.*NChipsPerHIC*1024.*fhr_run));
      cout<<"Run "<<runlabel[ir]<<" has "<<ntrig[ntrig.size()-1]<<" triggers"<<endl;
    }
  }

  //Start masking hottest pixels for each stave in each run, Fill also the histo with the hot pixel maps for each stave
  cout<<endl;
  cout<<"... Analysing FHR with (hot) pixel masking (Making also hot pixel map)"<<endl;
  vector<array<float,nMasked+1>> fhrall;
  vector<array<float,nMasked+1>> fhrall1;

  TH2F *hHotMap[nLayers][48];
  TH2F *hHotMapCloned[nLayers][48];

  Int_t index =0;
  for(int ilay=0; ilay<nLayers; ilay++){
    for(int istave=0; istave<48; istave++){
      if (isHotPixelMapDrawn) hHotMap[ilay][istave] = new TH2F(Form("hHotMap_L%s_Stv%d",laynums[index].c_str(), istave), "; ; ", (int)hmaps[index] ->GetAxis(0)->GetXmax(),0,hmaps[index] ->GetAxis(0)->GetXmax(), (int)hmaps[index] ->GetAxis(1)->GetXmax(),0,hmaps[index] ->GetAxis(1)->GetXmax());
    }
    index += nRuns*nStavesInLay[ilay];
  }

  Int_t numStavePart = 2;
  int irun = nRuns-1;
  Int_t nChips =0;
  Bool_t IB = 1;
  Int_t ilayEff =0;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //start from the bottom in order to start with the oldest run
    if (stoi(laynums[ihist]) < 3)  { //IB layers                                                   
      nChips=9; 
      MaxRange = 9216; // 1024 x nChips=9
      IB = 1;
      numStavePart = 1;
    }
    else if (stoi(laynums[ihist]) ==3 || stoi(laynums[ihist]) ==4) { //L3 and L4 
      nChips=2*28.;  //#chips in a half stave: 7 chips in one HIC row x 2 chip rows in one HIC x 4 HICs in a half stave
      MaxRange = 28672; //1024 x 7 chips in one HIC row x 4 HICs in a half stave 
      IB = 0;
      numStavePart =2; //HS Lower and HS Upper
    }
    else if (stoi(laynums[ihist]) ==5 || stoi(laynums[ihist]) ==6)  { //L5 and L6 
      nChips=2*49.;  //#chips in a half stave: 7 chips in one HIC row x 2 chip rows in one HIC x 7 HICs in a half stave
      MaxRange = 50176; //1024 x 7 chips in one HIC row x 7 HICs in a half stave 
      IB = 0;
      numStavePart =2; //HS Lower and HS Upper
    }

    //cout<<"\n*************Layer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<<" Run: "<<runnumbers[ihist]<< endl;
    for (Int_t StavePart=0; StavePart< numStavePart; StavePart++){ //loop over the two Half Staves for OB      
      int nchipsactive = GetNchipsActive(hmaps[ihist],nChips, MaxRange, StavePart, IB);
      if (numStavePart==1){
      if(nchipsactive<nChips) cout<<"\nLayer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<<" Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nchipsactive<<endl;
      }
      else {
	if(nchipsactive<nChips) cout<<"\nLayer "<<laynums[ihist]<<" Stave "<<stavenums[ihist]<< " Half Stave " << StavePart << " Run: "<<runnumbers[ihist]<<" --> Chips active:"<<nchipsactive<<endl;
      }
      ilayEff = stoi(laynums[ihist]);
      if (numStavePart == 2) ilayEff = stoi(laynums[ihist]) - 3;

      if (StavePart==0)   {
	fhrall.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : ilayEff][stoi(stavenums[ihist])], StavePart, IB, isHotPixelMapDrawn));
	if (numStavePart==1) fhrall1.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : ilayEff][stoi(stavenums[ihist])], StavePart, IB,  isHotPixelMapDrawn));
      }
      else  {
	fhrall1.push_back(GetFHRwithMasking(hmaps[ihist],nchipsactive,ntrig[irun],hHotMap[nLayers==1 ? 0 : ilayEff][stoi(stavenums[ihist])], StavePart, IB,  isHotPixelMapDrawn));
      }
    }
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

  TH2F *hFhrStv[nLayers][100][2];
  Int_t indexEff = 0;
  Int_t indexLaynums = 0;
  for(int ilay=0; ilay<nLayers; ilay++){
    if (ilay==0) indexLaynums = 0;
    else {
      indexEff += nStavesInLay[ilay-1];
      indexLaynums = nRuns*indexEff;
    }
    for(int is=0; is<nStavesInLay[ilay]; is++){
      for (Int_t i=0; i<2; i++){
	hFhrStv[ilay][is][i] = new TH2F(Form("h2FhrStv_%s_%d_HS%i", laynums[indexLaynums].c_str(),is, i), Form("Layer-%s - Stave-%d; # Hot Pixel Masked;Run", laynums[indexLaynums].c_str(),is),nMasked+1, binsmasked, nRuns, 0.5, nRuns+0.5);
      }
    }
  }

  //Fill histogram
  int ilayer = nLayers-1;
  int istave = nStavesInLay[ilayer]-1;
  irun=0;
  for(int i=0; i<(int)fhrall.size(); i++){
    for(int ifhr=0; ifhr<(int)fhrall[i].size(); ifhr++){
      for (Int_t StavePart=0; StavePart< 2; StavePart++){ //loop over the two Half Staves for OB  
	if (StavePart==0) hFhrStv[ilayer][istave][StavePart]->SetBinContent(ifhr+1, irun+1, fhrall[i][ifhr]);
	else hFhrStv[ilayer][istave][StavePart]->SetBinContent(ifhr+1, irun+1, fhrall1[i][ifhr]);
      }
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

  cout << "\nFHR vs #masked pixels " << endl;
  //Make FHR (averaged on all runs) vs #masked pix for all staves in a layer
  indexEff = 0;
  indexLaynums = 0;
  for(int ilay=0; ilay<nLayers; ilay++){
    if (ilay==0) indexLaynums = 0;
    else {
      indexEff += nStavesInLay[ilay-1];
      indexLaynums = nRuns*indexEff;
    }
    //legend
    TLegend *leg = new TLegend(0.904, 0.127,0.997,0.898);
    if (nLayers==4 && (ilay==1 || ilay==2 || ilay==3)) leg->SetNColumns(2);
    else if (nLayers==7 && ilay>=3) leg->SetNColumns(2);
    else  if (stoi(laynums[0]) ==4 || stoi(laynums[0]) ==5 || stoi(laynums[0]) ==6) leg->SetNColumns(2);

    if (IBorOB==0) numStavePart=1;
    else if (IBorOB==1) numStavePart=2;
    else {
      if (stoi(laynums[ilay]) < 3) numStavePart=1;
      else numStavePart=2;
    }
    TString SStavePart[2] = {"HS Lower", "HS Upper"};
    if (numStavePart==1) {
      SStavePart[0] = "";
      SStavePart[1] = "";
    }

    for (Int_t StavePart=0; StavePart< numStavePart; StavePart++){ //loop over the two Half Staves for OB      
      TCanvas cnv(Form("cnv_%d_HS%i",ilay, StavePart), Form("cnv_%d_HS%i",ilay, StavePart));
      cnv.cd();
      cnv.SetLogy();
      cnv.SetLogx();
      cnv.SetTickx();
      cnv.SetTicky();
      cnv.SetMargin(0.0988,0.1,0.1,0.0993);

      TH1F *hframe = cnv.DrawFrame(0.1,7e-15,3*(nMasked),1e-3, Form("Layer %s %s - Average FHR %s; # Hot Pixels masked ; FHR (/event/pixel)",laynums[indexLaynums].c_str(), SStavePart[StavePart].Data(), filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str()));
      hframe->SetBins(nMasked+1,binsmasked);

      for(int is=0; is<nStavesInLay[ilay];is++){
	TH1F *proj = (TH1F*)hFhrStv[ilay][is][StavePart]->ProjectionX(Form("proj_%d%d_HS%i",ilay,is, StavePart));
	int runswohits = GetNrunsWOhits(hFhrStv[ilay][is][StavePart]);
	proj->Scale(1./(nRuns-runswohits)); //Divide by the number of runs minus the ones without hits

	if ((nLayers==4 && (ilay>=0 && ilay<=3)) || stoi(laynums[0]) < 3) {
	  SetStyle(proj, col[is<nStavesInLay[ilay]/2 ? is : is-nStavesInLay[ilay]/2],is<nStavesInLay[ilay]/2 ? 24:26);
	}
	else if ((nLayers==7 && ilay>=3 && ilay<=4 ) || (stoi(laynums[0])==3 || stoi(laynums[0])==4)){
	  if((is)<nStavesInLay[ilay]/3){
	    SetStyle(proj, col[is], 24);
	  }
	  else if ((is)<nStavesInLay[ilay]*2/3){
	    SetStyle(proj, col[is-nStavesInLay[ilay]/3], 26);
	  }
	  else{
	    SetStyle(proj, col[is-nStavesInLay[ilay]*2/3], 25);
	  }

	}
	else if ((nLayers==7 && ilay>=5 && ilay<=6 ) || (stoi(laynums[0])==5 || stoi(laynums[0])==6)){
	  if((is)<int(nStavesInLay[ilay]/4)) {
	    SetStyle(proj, col[is], 24);
	  }
	  else if ((is)<2*int(nStavesInLay[ilay]/4))
	    SetStyle(proj, col[is-int(nStavesInLay[ilay]/4)], 26);
	  else if ((is)<3*nStavesInLay[ilay]/4)
	    SetStyle(proj, col[is-2*int(nStavesInLay[ilay]/4)], 25);
	  else
	    SetStyle(proj, col[is-int(nStavesInLay[ilay]*3/4)], 30);
	}

	proj->Draw("PL same");

	if (StavePart==0) leg->AddEntry(proj, Form("Stv%d",is),"p");
      }
      leg->Draw("same");
      TString NameCnv = "";

      if (StavePart==0){
	if (numStavePart==1) NameCnv = Form("../Plots/Layer%s_FHRpixmask_%s", laynums[indexLaynums].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
	else NameCnv = Form("../Plots/Layer%s_HSLower_FHRpixmask_%s", laynums[indexLaynums].c_str(), filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
      }
      else NameCnv = Form("../Plots/Layer%s_HSUpper_FHRpixmask_%s", laynums[indexLaynums].c_str(), filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
      cnv.SaveAs(NameCnv + ".pdf");
      cnv.SaveAs(NameCnv + ".root");
    }
  }

  indexEff=0;
  indexLaynums =0;
  //  cout << "Running time up to here: " << endl;
  //  t1.Stop();
  //  t1.Print();
  if (isHotPixelMapDrawn){
    cout << "\nDrawing hot pixel map for each layer " << endl;
    for(int ilay=0; ilay<nLayers; ilay++){
      if (ilay==0) indexLaynums = 0;
      else {
	indexEff += nStavesInLay[ilay-1];
	indexLaynums = nRuns*indexEff;
      }
      TCanvas cnvT(Form("cnvT_%d",ilay), Form("cnv_%d",ilay),800,1200);
      TCanvas cnvB(Form("cnvB_%d",ilay), Form("cnvBOT_%d",ilay),800,1200);
      cnvT.SetTopMargin(0.4);
      if (numStavePart==1) cnvT.Divide(1,nStavesInLay[ilay],0,0);
      else cnvT.Divide(1,nStavesInLay[ilay]/2,0,0);
      cnvB.SetTopMargin(0.4);
      cnvB.Divide(1,nStavesInLay[ilay]/2,0,0);
      if (numStavePart==1) cnvT.SetTitle(Form("cnvTOP_%d",ilay));

      int istaveeff=0;
      for(int istave=0; istave<nStavesInLay[ilay]; istave++){
	istaveeff = istave;
	hHotMap[ilay][istave]->SetMarkerStyle(20);
	hHotMap[ilay][istave]->SetMarkerSize(0.6);
	hHotMap[ilay][istave]->SetMarkerColor(kRed);
	hHotMap[ilay][istave]->SetLineColor(kRed);

	if (numStavePart==1 || istave < nStavesInLay[ilay]/2) {
	  cnvT.cd(istave+1);
	  cnvT.GetPad(istave+1)->SetTickx();
	  cnvT.GetPad(istave+1)->SetTicky();
	  cnvT.GetPad(istave+1)->SetRightMargin(0.01);
	  if(!istave) cnvT.GetPad(istave+1)->SetTopMargin(0.3);
	}
	else {
	  istaveeff = istave - nStavesInLay[ilay]/2;
	  cnvB.cd(istaveeff+1);
	  cnvB.GetPad(istaveeff+1)->SetTickx();
	  cnvB.GetPad(istaveeff+1)->SetTicky();
	  cnvB.GetPad(istaveeff+1)->SetRightMargin(0.01);
	  if(!istaveeff) cnvB.GetPad(istaveeff+1)->SetTopMargin(0.3);
	}
	hHotMap[ilay][istave]->Rebin2D(4);
	hHotMapCloned[ilay][istave]= (TH2F*) hHotMap[ilay][istave]->Clone(Form("%s_", hHotMap[ilay][istave]->GetName()));
	hHotMap[ilay][istave]->Draw("P X+");
	hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.005);
	hHotMap[ilay][istave]->GetYaxis()->SetTickLength(0.005);
	hHotMap[ilay][istave]->GetYaxis()->SetLabelSize(0.13);
	hHotMap[ilay][istave]->GetXaxis()->SetLabelSize(0.13);
	if ( (numStavePart==1 && istave ==0) || (numStavePart ==2 && (istave == 0 || istave == nStavesInLay[ilay]/2)) ){
	  hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(0.003);
	  hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
	  hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
	}
	else{
	  hHotMap[ilay][istave]->GetXaxis()->SetLabelOffset(999);
	  hHotMap[ilay][istave]->GetXaxis()->SetTickLength(0.05);
	  hHotMap[ilay][istave]->GetXaxis()->SetNdivisions(530);
	}
	if (numStavePart == 2)   hHotMap[ilay][istave]->GetYaxis()->SetNdivisions(8);

	TLatex lat;
	lat.SetTextAngle(90);
	lat.SetNDC();
	if (numStavePart==1) {
	  lat.SetTextSize(0.15); 
	  lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
	}
	else {
	  if (istave < nStavesInLay[ilay]/2){
	    lat.SetTextSize(0.25);
	    lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
	  }
	  else {
	    lat.SetTextSize(0.25);
	    lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));
	  }
	}
      }

      cnvT.cd();
      TLatex latT;
      latT.SetNDC();
      latT.SetTextSize(0.03);
      if (numStavePart==1)    latT.DrawLatex(0.01,0.98,Form("L%s",laynums[indexLaynums].c_str()));
      else     latT.DrawLatex(0.01,0.98,Form("L%s - TOP",laynums[indexLaynums].c_str()));

      if (numStavePart==2){
	cnvB.cd();
	TLatex latB;
	latB.SetNDC();
	latB.SetTextSize(0.03);
	latB.DrawLatex(0.01,0.98,Form("L%s - BOTTOM",laynums[indexLaynums].c_str()));
      }
   
      TFile * fileOutputT; 
      TFile * fileOutputB; 
      TString NameCnvT="";
      TString NameCnvB="";
      if (numStavePart==1) {
	NameCnvT = Form("../Plots/Layer%s_Hotpixmap_%s", laynums[indexLaynums].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
      }
      else {
	NameCnvT = Form("../Plots/Layer%s-TOP_Hotpixmap_%s", laynums[indexLaynums].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
	NameCnvB = Form("../Plots/Layer%s-BOT_Hotpixmap_%s", laynums[indexLaynums].c_str(),filepath_hit.substr(filepath_hit.find("from"), filepath_hit.find(".root")-filepath_hit.find("from")).c_str());
	cnvB.SaveAs(NameCnvB + ".pdf");
	fileOutputB = new TFile(NameCnvB+ ".root", "RECREATE");
      }
      cnvT.SaveAs(NameCnvT + ".pdf");
      fileOutputT = new TFile(NameCnvT + ".root", "RECREATE");

      for(int istave=0; istave<nStavesInLay[ilay]; istave++){
	if (numStavePart==1 || istave < nStavesInLay[ilay]/2) {
	  fileOutputT->WriteTObject(hHotMapCloned[ilay][istave]);
	}
	else {
	  fileOutputB->WriteTObject(hHotMapCloned[ilay][istave]);
	}
      }
      fileOutputT->Close();
      if (numStavePart==2)    fileOutputB->Close();
      cout << "The following root files have been created:\n" << NameCnvT << ".root" << endl;
      if (numStavePart==2) cout << NameCnvB << ".root"<< endl;
    }
  }
  //t.Stop();
  //cout << "\nRunning time: " << endl;
  //t.Print();
}

//
// Function to compare two hitmaps --> returns an arrays with timestamp of run2, noisyPixInRefRun, noisyPixInRun2, noisyPixInCommon
//
std::array<float,nMasked+1> GetFHRwithMasking(THnSparse *hmap, const int nchips, double ntrig, TH2 *hhotmap, bool HS , bool IB, bool isHotPixelMapDrawn){

  array<float,nMasked+1> fhrstave;

  THnSparse *hmapclone = (THnSparse*)hmap->Clone(Form("%s_clone",hmap->GetName()));
  //cout<<"\n" <<hmap->GetName()<<" --> "<<hmapclone->GetNbins()<<" noisy pixels"<<endl;
  Int_t iyMin =1;
  Int_t iyMax =hmapclone->GetAxis(1)->GetNbins();
  if (IB==0) {
    if (HS==0) iyMax = hmap->GetAxis(1)->GetNbins()/2;
    else iyMin = hmap->GetAxis(1)->GetNbins()/2;
  }

  hmapclone->GetAxis(1)->SetRange(iyMin, iyMax);
  for(int iter=0; iter<nMasked+1; iter++){
    TH1F *hproj = (TH1F*)hmapclone->Projection(1); 
    long int totalhits = hproj->Integral();
    float fhr = nchips==0 ? 0. : (float)totalhits / (512.*1024.*nchips*ntrig); 

    if(ntrig<0) fhr=0.;
    fhrstave[iter] = fhr;

    int coord[2];
    double max = -1.;
    int x=0,y=0;
    long int binwithmax = 0;
    for(int ibin=0; ibin<hmapclone->GetNbins(); ibin++){
      double bincontent = hmapclone->GetBinContent(ibin, coord);
      if (IB==0){
	if (coord[1] < iyMin || coord[1] > iyMax) continue;
      }
      if(bincontent>max){
        max=bincontent;
        binwithmax = ibin;
        x=coord[0];
        y=coord[1];
      }
    }

    if(nchips) hmapclone->SetBinContent(binwithmax,0.);
    if (isHotPixelMapDrawn){
      if(totalhits!=0 && iter < 100) {
	hhotmap->SetBinContent(x-1, y-1, 1);// to avoid a marker in 0,0 for empty histos
      }
    }
    delete hproj;
  } 
  
  delete hmapclone;

  return fhrstave;

}

//
//Function to return the number of active (i.e. enabled) chips in a stave (or in a half stave when dealing with OB)
//
int GetNchipsActive(THnSparse *hmap, int maxchip, int MaxRange, bool HS, bool IB){
  int ix1=1;
  int ix2=1024;
  int iy1=1;
  int iy2=512;
  int iy2Max = 0;
  int iy2Min = 1;
  int activechips = maxchip;
  iy2Max =hmap->GetAxis(1)->GetNbins();
  if (IB==0) {
    if (HS==0) iy2Max = hmap->GetAxis(1)->GetNbins()/2;
    else {
      iy1 = 1025;
      iy2 = 1024 + 512;
      iy2Min = hmap->GetAxis(1)->GetNbins()/2;
    }
  }
  while (iy2 <= iy2Max && iy1 >= iy2Min){
    ix1 = 1;
    ix2=1024;
    while(ix2<=hmap->GetAxis(0)->GetNbins()){
      hmap->GetAxis(0)->SetRange(ix1,ix2);
      hmap->GetAxis(1)->SetRange(iy1,iy2);
      TH2F *hproj = (TH2F*)hmap->Projection(1,0);
      if((hproj->Integral(1,1024,1,512)<1e-15)){
	activechips--;
      }
      ix1=ix2+1;
      ix2+=1024;
      delete hproj;
    }
    iy1=iy2+1; 
    iy2+=512;
  }
  hmap->GetAxis(0)->SetRange(1,MaxRange);//reset range
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
