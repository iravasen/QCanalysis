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
void DoAnalysis(string filepath, const int nChips, bool isIB);

//
// MAIN
//
void AnalyzeErrorsFHR(){
  string fpath;
  int nchips=9;
  cout<<"\n\n=> Available file(s) for the analysis (the last should be the file you want!): \n"<<endl;
  gSystem->Exec("ls ../Data/*w_error_data* -Art | tail -n 500");
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


  //Call
  DoAnalysis(fpath, nchips, isIB);

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
void DoAnalysis(string filepath, const int nChips, bool isIB){

  gStyle->SetOptStat(0000);

  std::vector<TH2*> herr;
  std::vector<string> timestamps, runnumbers;
  int nTimes=0, nRuns=1;

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
    if(objname.find("err")==string::npos) continue;
    h2 = (TH2*)obj;
    if(!h2->GetEntries()) continue;
    cout<<"... Reading "<<obj->GetName()<<endl;
    herr.push_back(h2);
    string timestamp = objname.find("run")==string::npos ? objname.substr(objname.find("_",2)+1, 13) : objname.substr(objname.find("_",6)+1, 13);
    string runnum =  objname.find("run")==string::npos ? "norun":objname.substr(objname.find("run")+3, 6);
    timestamps.push_back(timestamp);
    runnumbers.push_back(runnum);
    nTimes++;
    nRuns++;
  }

  //sum all the histos in a single histogram (for summary plot)
  TH2D *hSummary = (TH2D*)herr[0]->Clone("hSummary");
  for(int iplot=1; iplot<(int)herr.size(); iplot++)
    hSummary->Add(herr[iplot]);

  //Draw summary plot
  TCanvas *canvas = new TCanvas();
  canvas->Divide(1,2);
  canvas->cd(1);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->SetMargin(0.0988,0.2,0.194,0.0993);
  canvas->GetPad(1)->SetRightMargin(0.15);
  hSummary->SetTitle(Form("Errors IB, %s", filepath.substr(filepath.find("from"), filepath.find("_w_")-filepath.find("from")).c_str()));
  hSummary->Draw("colz");
  hSummary->GetXaxis()->SetNdivisions(530);
  hSummary->GetYaxis()->SetNdivisions(516);
  hSummary->GetXaxis()->SetLabelSize(0.045);
  hSummary->GetYaxis()->SetLabelSize(0.045);
  hSummary->GetZaxis()->SetLabelSize(0.045);
  hSummary->GetXaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleSize(0.05);
  hSummary->GetYaxis()->SetTitleOffset(0.7);
  hSummary->GetZaxis()->SetTitleSize(0.05);
  hSummary->GetZaxis()->SetTitleOffset(0.7);
  canvas->cd(2);
  TLegend *legend = new TLegend(0.1,0.1,0.9,0.9);
  legend->SetHeader("Error ID legend:");
  legend->AddEntry("", "0: RDH page counters for the same RU/trigger are not continuous","");
  legend->AddEntry("", "1: RDH ang GBT header page counters are not consistent","");
  legend->AddEntry("", "2: GBT payload header was expected but not found","");
  legend->AddEntry("", "3: GBT payload trailer was expected but not found","");
  legend->AddEntry("", "4: All lanes were stopped but the page counter in not 0","");
  legend->AddEntry("", "5: End of FEE data reached while not all lanes received stop","");
  legend->AddEntry("", "6: Data was received for stopped lane","");
  legend->AddEntry("", "7: No data was seen for lane (which was not in timeout)","");
  legend->AddEntry("", "8: ChipID (on module) was different from the lane ID on the IB stave","");
  legend->AddEntry("", "9: Cable data does not start with chip header or empty chip","");
  legend->AddEntry("", "10: Active lanes pattern conflicts with expected for given RU type","");
  legend->AddEntry("", "11: Jump in RDH_packetCounters","");
  legend->Draw();

  canvas->SaveAs(Form("../Plots/IB_FHRErrorPlotSummary_%s.pdf", filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  canvas->SaveAs(Form("../Plots/IB_FHRErrorPlotSummary_%s.root",filepath.substr(filepath.find("from"), filepath.find(".root")-filepath.find("from")).c_str()));
  delete canvas;
}
