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
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns);

//
// MAIN
//
void AnalyzeTrgFlgFHR(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*w_error_and_trig* -Art | tail -n 500");
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


  //Call
  DoAnalysis(fpath, nchips, isIB, skipruns);

}

//
//Set Style
//
void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.2);
  h->SetMarkerColor(col);
  //h->SetFillStyle(0);
  //h->SetFillColorAlpha(col,0.8);
}

//
// Analyse data
//
void DoAnalysis(string filepath, const int nChips, bool isIB, string skipruns){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> herr;
  std::vector<string> timestamps, runnumbers, laynums;
  int nRuns = 1;
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
       ) {
            cout<<"<W> Object "<<obj->GetName()<<" is not 1D or 2D histogram : will not be converted"<<endl;
       }
    string objname = (string)obj->GetName();
    if(objname.find("trg")==string::npos) continue;
    h2 = (TH2*)obj;
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    string laynum = objname.substr(objname.find("L")+1,1);

    if(skipruns.find(runnum)!=string::npos) continue; //eventually skip runs specified by the user

    cout<<"... Reading "<<obj->GetName()<<endl;
    herr.push_back(h2);
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    laynums.push_back(laynum);

    if(laynums.size()>1 && laynum==laynums[laynums.size()-2])
      nRuns++;
    else nRuns=1;
  }

  const int nLayers = (int)herr.size()==nRuns ? 1 : stoi(laynums[laynums.size()-1])+1;

  //sum all the histos in a single histogram (for summary plot) for each layer
  TH2D *hSummary[nLayers];
  for(int ilay=0; ilay<nLayers; ilay++)
    hSummary[ilay] = (TH2D*)herr[stoi(laynums[ilay*nRuns])]->Clone(Form("hSummary_L%d",nLayers>1 ? ilay:stoi(laynums[0])));
  for(int iplot=1; iplot<(int)herr.size(); iplot++){
    int layidx = stoi(laynums[iplot]);
    if(iplot == layidx*nRuns) continue;
    hSummary[layidx]->Add(herr[iplot]);
  }

  //Make plots with Error IDs vs Run for each layer
  TGraph *trend[nLayers][hSummary[0]->GetNbinsY()];

  int ir = 0;
  double max[nLayers];
  for(int ilay=0; ilay<nLayers; ilay++)
    max[ilay] = -1.;
  for(int iplot=0; iplot<(int)herr.size(); iplot++){
    TH1D *hproj = (TH1D*)herr[iplot]->ProjectionY(Form("herr_%d",iplot));
    for(int ibin=1; ibin<=hSummary[0]->GetNbinsY(); ibin++){
      if(ir==0){
        trend[stoi(laynums[iplot])][ibin-1] = new TGraph();
        trend[stoi(laynums[iplot])][ibin-1]->SetName(Form("gr_L%s_trgID%d",laynums[iplot].c_str(),ibin));
        SetStyle(trend[stoi(laynums[iplot])][ibin-1], col[ibin<=10?ibin-1:ibin-11], ibin<=10?24:25);
      }
      trend[stoi(laynums[iplot])][ibin-1]->SetPoint(ir,ir, hproj->GetBinContent(ibin));
      if(hproj->GetBinContent(ibin)>max[stoi(laynums[iplot])])
        max[stoi(laynums[iplot])]=hproj->GetBinContent(ibin);
    }
    delete hproj;
    ir++;
    if(ir==nRuns) ir=0;
  }

  //Draw summary plot
  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas canvas;
    canvas.cd();
    canvas.SetTickx();
    canvas.SetTicky();
    canvas.SetLogz();
    canvas.SetMargin(0.18,0.2,0.194,0.0993);
    canvas.SetRightMargin(0.15);
    hSummary[ilay]->SetTitle(Form("Trigger & Flags L%d, %s", nLayers>1 ? ilay:stoi(laynums[0]),filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
    hSummary[ilay]->Draw("colz");
    //hSummary[ilay]->GetXaxis()->SetNdivisions(530);
    //hSummary[ilay]->GetYaxis()->SetNdivisions(516);
    hSummary[ilay]->GetXaxis()->SetLabelSize(0.045);
    hSummary[ilay]->GetYaxis()->SetLabelSize(0.045);
    hSummary[ilay]->GetZaxis()->SetLabelSize(0.045);
    hSummary[ilay]->GetXaxis()->SetTitleSize(0.05);
    hSummary[ilay]->GetYaxis()->SetTitleSize(0.05);
    hSummary[ilay]->GetYaxis()->SetTitleOffset(0.7);
    hSummary[ilay]->GetZaxis()->SetTitleSize(0.05);
    hSummary[ilay]->GetZaxis()->SetTitleOffset(0.9);

    if(!ilay) canvas.SaveAs(Form("../Plots/IB_FHRTrgFlgPlotSummary_%s.pdf[", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    canvas.SaveAs(Form("../Plots/IB_FHRTrgFlgPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

  //Draw trends
  int npoints = trend[0][0]->GetN();
  TH1F *hfake = new TH1F("hfake", "; Run; # Errors", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<(int)runnumbers.size()/nLayers; ir++)
      hfake->GetXaxis()->SetBinLabel(ir+1, Form("run%06d", stoi(runnumbers[(int)runnumbers.size()/nLayers-1-ir])));

  TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
  leg->SetHeader("IDs");
  for(int iid=1; iid<=hSummary[0]->GetNbinsY();iid++)
    leg->AddEntry(trend[0][iid-1], Form("%s",hSummary[0]->GetYaxis()->GetBinLabel(iid)), "p");

  for(int ilay=0; ilay<nLayers; ilay++){
    TCanvas canvas;
    canvas.cd();
    canvas.SetTickx();
    canvas.SetTicky();
    canvas.SetLogy();
    canvas.SetMargin(0.0988,0.1,0.194,0.0993);

    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetTitle(Form("Layer-%d, Trigger & Flag trends %s",nLayers>1 ? ilay:stoi(laynums[0]), filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    hfake->GetYaxis()->SetRangeUser(1, 10*max[ilay]);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->Draw();
    for(int iid=1; iid<=hSummary[ilay]->GetNbinsY();iid++){
      trend[ilay][iid-1]->Draw("P same");
    }
    leg->Draw("same");
    canvas.SaveAs(Form("../Plots/IB_FHRTrgFlgPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
    if(ilay==nLayers-1) canvas.SaveAs(Form("../Plots/IB_FHRTrgFlgPlotSummary_%s.pdf]", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  }

}
