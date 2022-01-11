#include <string>
#include <iostream>
#include <vector>
#include <array>
#include "inc/itsAnalysis.hh"

void MakeDeadPixelMap(){
  itsAnalysis myAnalysis("Dead pixel Hits"); // Change to "Hits on Layer" if using FHR as fake data

  auto nLayers      = myAnalysis.nLayers();     // int of number of layers
  auto laynums      = myAnalysis.Layers();      //vec of layers
  auto runNumbers   = myAnalysis.Runs();        //vec of run numbers
  auto hmaps        = myAnalysis.loadedHists(); // all histograms for layers and runs needed
  auto nRuns        = myAnalysis.nRuns();

  std::unique_ptr<TFile> myFile( TFile::Open(Form("DeadPixMapResults_run%s_to_run%s.root",runNumbers[nRuns-1].c_str(),runNumbers[0].c_str()), "RECREATE") );

  for (string layer : laynums){ // loop over layers
    int ilay=stoi(layer);
    int nStavesInLay = myAnalysis.stavesInLayer(ilay);
    
    TCanvas cnv(Form("cnv_%d",ilay), Form("cnv_%d",ilay),800,1200);

    cnv.SetTopMargin(0.4);
    if (ilay>=2) cnv.Divide(1,20,0,0);
    else cnv.Divide(1,nStavesInLay,0,0);

    int iHists = -1;

    for(int istave=0; istave<nStavesInLay; istave++){
      iHists ++;
      TH2F *hHotMap;
      if (ilay==3 || ilay==4) TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%i",layer.c_str(),istave), "; ; ", 28672/4,0.5,28672.5,2048/4,0.5,2048.5);
      else if (ilay>=5) TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%i",layer.c_str(),istave), "; ; ", 50176/4,0.5,50176.5,2048/4,0.5,2048.5);
      else TH2F *hHotMap = new TH2F(Form("hHotMap_L%s_Stv%i",layer.c_str(),istave), "; ; ", 9216,0.5,9216.5,512,0.5,512.5);

      int cnt = 0;

      for(auto hist : myAnalysis.loadLayerSparse(stoi(layer))){ //loop hist for layer
        if(hist->GetEntries()>1e4){
          if(istave==0)cout<<hist->GetName()<<" skipped because has more than 10000 entries (probably a bad run)."<<endl;
          continue;
        }
        int stave = myAnalysis.getStaveNumber(hist);
        if(stave!=istave) continue;
        if(!cnt) {hHotMap =(TH2F*)hist->Projection(1,0);
                  if(ilay>=3) hHotMap->Rebin2D(4,4); 
        }
        else {
          TH2F *htemp = (TH2F*)hist->Projection(1,0);
          if(ilay>=3) htemp->Rebin2D(4,4);
          hHotMap->Add(htemp);
          delete htemp;
        }
        cnt++;
      }

      hHotMap->SetTitle(Form("hHotMap_lay%i_stv%i",ilay,istave));
      hHotMap->SetName(Form("hHotMap_lay%i_stv%i",ilay,istave));
      hHotMap->SetMarkerStyle(20);
      hHotMap->SetMarkerSize(0.6);
      hHotMap->SetMarkerColor(kRed);
      hHotMap->SetLineColor(kRed);
      hHotMap->Write();

      cnv.cd(iHists+1);
      cnv.GetPad(iHists+1)->SetTickx();
      cnv.GetPad(iHists+1)->SetTicky();
      cnv.GetPad(iHists+1)->SetRightMargin(0.01);
      if(!iHists) cnv.GetPad(iHists+1)->SetTopMargin(0.1);
      hHotMap->SetTitle(" ");

      hHotMap->GetXaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetTickLength(0.005);
      hHotMap->GetYaxis()->SetLabelSize(0.13);
      hHotMap->GetXaxis()->SetLabelSize(0.13);

      if(iHists>0){
        hHotMap->GetXaxis()->SetLabelOffset(999);
        hHotMap->GetXaxis()->SetTickLength(0.05);
        hHotMap->GetXaxis()->SetNdivisions(530);
      }
      else{
        hHotMap->GetXaxis()->SetLabelOffset(0.003);
        hHotMap->GetXaxis()->SetTickLength(0.05);
        hHotMap->GetXaxis()->SetNdivisions(530);
      }
      hHotMap->SetStats(0);
      hHotMap->DrawCopy("P X+");

      TLatex lat;
      lat.SetTextAngle(90);
      lat.SetNDC();
      lat.SetTextSize(0.15);
      lat.DrawLatex(0.04,0.3,Form("Stv%d",istave));

      delete hHotMap;

      cout<<"Layer: "<<layer<<" stv: "<<istave<<" / "<<nStavesInLay-1<<endl;

      if(((istave+1)%20==0 && istave!=0) || (istave+1)==nStavesInLay){
          cnv.cd();
          TLatex lat;
          lat.SetNDC();
          lat.SetTextSize(0.03);
          lat.DrawLatex(0.01,0.98,Form("L%s",layer.c_str()));
          cnv.SaveAs(Form("../Plots/Layer%i_Deadpixmap_run%s_to_run%s_upto_stv%i.pdf",ilay,runNumbers[nRuns-1].c_str(),runNumbers[0].c_str(),istave));

          cnt=0;
          iHists = -1;
      }

    } // End loop on staves
  } // End loop on layers
} // end of qMakeDeadPixelMap()