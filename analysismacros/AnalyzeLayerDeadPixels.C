#include <string>
#include <iostream>
#include <vector>
#include <TH2.h>
#include <TFile.h>
#include <TList.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMath.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TKey.h>

using namespace std;

void SetStyle(TGraph *h, Int_t col, Style_t mkr);
void DoAnalysis(string filepath, const int nChips, string skipruns, bool isIB);

//
// MAIN
//
void AnalyzeLayerDeadPixels(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
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


  //Call
  DoAnalysis(fpath, nchips, skipruns, isIB);

}

//
//Set Style
//
void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, string skipruns, bool isIB){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> hmaps;
  std::vector<string> timestamps, runnumbers, laynums;
  int nTimes=0, nRuns=1;
  Int_t col[] = {TColor::GetColor("#ff3300"), TColor::GetColor("#ec6e0a"), TColor::GetColor("#daaa14"), TColor::GetColor("#c7e51e"), TColor::GetColor("#85dd69"), TColor::GetColor("#42d6b4"), TColor::GetColor("#00ceff"), TColor::GetColor("#009adf"), TColor::GetColor("#0067c0"), TColor::GetColor("#0033a1")};

  //Read the file and the list of plots with entries
  TFile *infile=new TFile(filepath.c_str());
  TList *list = infile->GetListOfKeys();
  TKey *key;
  TObject *obj;
  TIter next(list);
  TH2 *h2;
  while((key = ((TKey*)next()))){
    obj = key->ReadObj();
    if ((strcmp(obj->IsA()->GetName(),"TProfile")!=0)
         && (!obj->InheritsFrom("TH2"))
	       && (!obj->InheritsFrom("TH1"))
         && (!obj->InheritsFrom("THnSparse"))
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    string objtitle = (string)obj->GetTitle();
    if(objname.find("Stv")!=string::npos) break;
    if(objtitle.find("Threshold")!=string::npos) continue;//take only maps of #dead pixels

    h2 = (TH2*)obj;
    cout<<"... Reading "<<obj->GetName()<<endl;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);
    if(skipruns.find(runnum)!=string::npos) continue;// eventually skip runs if specified
    hmaps.push_back(h2);
    //cout<<runnum<<"  "<<timestamp<<endl;
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);

    laynums.push_back(laynum);
    //cout<<"run: "<<runnum<<"   timestamp: "<<timestamp<<"    laynum: "<<laynum<<endl;
    nTimes++;
    if(nTimes>1 && laynum==laynums[laynums.size()-2])
      nRuns++;
    else nRuns=1;
  }

  const int nLayers = (int)hmaps.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1;

  //const int nRuns = (int)runnumbers.size() / nLayers;
  //cout<<"Lay: "<<nLayers<<"  Runs: "<<nRuns<<endl;
  TGraph *trend[nLayers][100];
  int ilayer=nLayers-1;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){
    for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){
      trend[ilayer][ibiny-1] = new TGraph();
    }
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        ilayer--;
      }
  }
  ilayer=nLayers-1;
  TH1F *hproj = new TH1F();
  string histname = hmaps[0]->GetName();
  int irun=0;
  double maxtot=-1.;
  int staveswithdead[nLayers];
  for(int ilay=0; ilay<nLayers; ilay++)
    staveswithdead[ilay] = 0;
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){//loop on y bins (staves)
      TH1D *hproj = hmaps[ihist]->ProjectionX("proj",ibiny,ibiny); //single stave
      trend[ilayer][ibiny-1]->SetName(Form("gr_L%s_stave%d",laynums[ihist].c_str(),ibiny-1));
      int chipswithdeadpix = 0;
      for(int ibinx=1; ibinx<=hmaps[ihist]->GetNbinsX(); ibinx++){//evaluate the chips with dead pixels
        if(hmaps[ihist]->GetBinContent(ibinx,ibiny)>0)
          chipswithdeadpix++;
      }
      if(chipswithdeadpix>0){
        cout<<"Layer "<<laynums[ihist]<<" Stave "<<ibiny-1<<" Run: "<<runnumbers[ihist]<<" --> # Chips with dead pix: "<<chipswithdeadpix<<endl;
        staveswithdead[ilayer]++;
      }

      trend[ilayer][ibiny-1]->SetPoint(irun, irun, hproj->Integral());//total number of dead pix for this stave, if 0 pixels put 0.1 to allow log scale
      if(hproj->Integral()>maxtot)
        maxtot=hproj->Integral();

      if((ibiny-1)<hmaps[ihist]->GetNbinsY()/2)
        SetStyle(trend[ilayer][ibiny-1], col[ibiny-1], 24);
      else
        SetStyle(trend[ilayer][ibiny-1], col[ibiny-1-hmaps[ihist]->GetNbinsY()/2], 26);
    }
    irun++;
    if(ihist>0)
      if(laynums[ihist-1]!=laynums[ihist]){
        irun=0;
        ilayer--;
      }
  }

  int npoints = trend[0][0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; # Dead Pixels", npoints, -0.5, (double)npoints-0.5);

  for(int ir=0; ir<(int)runnumbers.size()/nLayers; ir++)
      hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[(int)runnumbers.size()/nLayers-1-ir])));

  //find max for plots (same for all staves)
  double max = -1.;
  for(int ihist=0; ihist<(int)hmaps.size(); ihist++){
    if(hmaps[ihist]->GetMaximum()>max)
      max=hmaps[ihist]->GetMaximum();
  }

  //Draw
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogy();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.0988,0.1,0.194,0.0993);
    TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
    for(int istave=0; istave<hmaps[ilay*nRuns]->GetNbinsY(); istave++)
      leg->AddEntry(trend[ilay][istave], Form("Stv%d",istave), "p");
    hfake->GetYaxis()->SetRangeUser(8e-1, maxtot+0.5*maxtot);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetTitle(Form("Layer-%s (%d Staves with dead pix), %s",laynums[ilay*nRuns].c_str(), (int)ceil(staveswithdead[ilay]/nRuns), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    hfake->Draw();
    for(int istave=0; istave<hmaps[ilay*nRuns]->GetNbinsY(); istave++)
      trend[ilay][istave]->Draw("P same");
    leg->Draw("same");
    canvas->SaveAs(Form("../Plots/Layer%s_deadpixels_%s.pdf", laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas->SaveAs(Form("../Plots/Layer%s_deadpixels_%s.root", laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    delete canvas;
    delete leg;
  }

  //Make GIF with TH2 for each run and for each layer
  irun=1;
  ilayer=nLayers-1;
  gStyle->SetPalette(1);
  for(int ilay=0; ilay<nLayers; ilay++)//remove images if they exist already
    gSystem->Unlink(Form("../Plots/Layer%s_deadpixelmap_%s.gif", laynums[ilay*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){// start from the last in order to have the runs from the oldest to the newest
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetRightMargin(0.15);
    hmaps[ihist]->Draw("colz");
    hmaps[ihist]->SetMinimum(0);
    hmaps[ihist]->SetMaximum(max);
    hmaps[ihist]->GetZaxis()->SetTitle("# Dead Pixels");
    hmaps[ihist]->SetTitle(Form("Layer-%s, Run %06d (%d/%d)", laynums[ilayer*nRuns].c_str(), stoi(runnumbers[ihist]), irun,nRuns));
    canvas->Print(Form("../Plots/Layer%s_deadpixelmap_%s.gif+40", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    irun++;
    if(!ihist && nLayers==1){
      canvas->Print(Form("../Plots/Layer%s_deadpixelmap_%s.gif++40++", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      break;
    }
    if(nLayers>1){
      if(ihist>0){
        if(laynums[ihist-1]!=laynums[ihist]){
          canvas->Print(Form("../Plots/Layer%s_deadpixelmap_%s.gif++40++", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
          ilayer--;
          irun=1;
        }
      }
      else if(!ihist && !ilayer){
        canvas->Print(Form("../Plots/Layer%s_deadpixelmap_%s.gif++40++", laynums[ilayer*nRuns].c_str(), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
      }
    }

    delete canvas;
  }

}
