#include "inc/itsAnalysis.hh"

void SetStyle(TGraph *h, Int_t col, Style_t mkr){
  h->SetLineColor(col);
  h->SetMarkerStyle(mkr);
  h->SetMarkerSize(1.4);
  h->SetMarkerColor(col);
}

Int_t col[] = {810, 807, 797, 827, 417, 841, 868, 867, 860, 602, 921, 874};

void AnalyzeLayerDeadPixels(){
  itsAnalysis myAnalysis("DeadPixel");

  auto nLayers      = myAnalysis.nLayers();     // int of number of layers
  auto laynums      = myAnalysis.Layers();      //vec of layers
  auto runNumbers   = myAnalysis.Runs();        //vec of run numbers
  auto nRuns        = myAnalysis.nRuns();
  auto hmaps        = myAnalysis.loadedHists(); // all histograms for layers and runs needed

  TGraph *trend[nLayers][100]; //list of trends
  for (int ilayer = 0; ilayer < 6; ++ilayer){ //loop over layers
    for(int ihist=(int)hmaps.size()-1; ihist>=0; ihist--){ //loop over number of histograms
      for(int ibiny=1; ibiny<=hmaps[ihist]->GetNbinsY(); ibiny++){ // Loop over y-bins (staves)
        trend[ilayer][ibiny-1] = new TGraph();
      }
    }
  }

  int staveswithdead[nLayers];
  double maxtot=-1.;

  TH1F *hproj = new TH1F();
  for (string layer : laynums){ // loop over layers
    staveswithdead[stoi(layer)] = 0;
    int irun = myAnalysis.nRuns()-1; // Count back to have right order in histogram
    for(auto hist : myAnalysis.loadLayer(stoi(layer))){ //loop hist for layer
      for(int ibiny=1; ibiny<=hist->GetNbinsY(); ibiny++){//loop on y bins (staves)
        TH1D *hproj = hist->ProjectionX("proj",ibiny,ibiny); //single stave
        trend[stoi(layer)][ibiny-1]->SetName(Form("gr_L%s_stave%d",layer.c_str(),ibiny-1));
        int nChips = myAnalysis.nChips(stoi(layer));
        int chipswithdeadpix = 0;
        for(int ibinx=1; ibinx<=hist->GetNbinsX(); ibinx++){//evaluate the number of disabled chips
          if(hist->GetBinContent(ibinx,ibiny)>0) chipswithdeadpix++;
        }
        if(chipswithdeadpix>0){
          cout<<"Layer: "<<layer<<" Stave: "<<ibiny-1<<" Run: "<<myAnalysis.getRunNumber(hist)<<" --> # Chips with dead pix: "<<chipswithdeadpix<<endl;
          staveswithdead[stoi(layer)]++;
        }

        trend[stoi(layer)][ibiny-1]->SetPoint(irun, irun, hproj->Integral());
        if(hproj->Integral()>maxtot)
          maxtot=hproj->Integral();

        // Set style IB
        if (stoi(layer) < 3){
          if((ibiny-1)<hist->GetNbinsY()/2)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1], 24);
          else
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-hist->GetNbinsY()/2], 26);
        }
        //Set Style OB
        else if (stoi(layer) ==3 || stoi(layer)==4){
          if((ibiny-1)<hist->GetNbinsY()/6)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1], 24);
          else if ((ibiny-1)<hist->GetNbinsY()*2/6)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()/6)], 26);
          else if ((ibiny-1)<hist->GetNbinsY()*3/6)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*2/6)], 25);
          else if ((ibiny-1)<hist->GetNbinsY()*4/6)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*3/6)], 24);
          else if ((ibiny-1)<hist->GetNbinsY()*5/6)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*4/6)], 26);
          else
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*5/6)], 25);
        }
        else if (stoi(layer) ==5 || stoi(layer)==6){
          if((ibiny-1)<hist->GetNbinsY()/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1], 24);
          else if ((ibiny-1)<hist->GetNbinsY()*2/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()/8)], 26);
          else if ((ibiny-1)<hist->GetNbinsY()*3/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*2/8)], 25);
          else if ((ibiny-1)<hist->GetNbinsY()*4/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*3/8)], 30);
          else if ((ibiny-1)<hist->GetNbinsY()*5/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*4/8)], 24);
          else if ((ibiny-1)<hist->GetNbinsY()*6/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*5/8)], 26);
          else if ((ibiny-1)<hist->GetNbinsY()*7/8)
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*6/8)], 25);
          else
            SetStyle(trend[stoi(layer)][ibiny-1], col[ibiny-1-(hist->GetNbinsY()*7/8)], 30);
        }
      }
      irun--;
    }
  }

  int npoints = myAnalysis.nRuns();
  TH1F *hfake = new TH1F("hfake", "; # Dead Pixels", npoints, -0.5, (double)npoints-0.5);
  for(int ir=0; ir<npoints; ir++) // Set bin labels to run numbers
      hfake->GetXaxis()->SetBinLabel(npoints-(ir), Form("run%06d",stoi(myAnalysis.Runs()[ir])));

  for (string layer : laynums){ // loop over layers
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetLogy();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetMargin(0.0988,0.1,0.194,0.0993);
    TLegend *leg = new TLegend(0.904, 0.197,0.997,0.898);
    hfake->GetYaxis()->SetRangeUser(8e-1, maxtot+0.5*maxtot);
    hfake->GetXaxis()->SetTitleOffset(2.8);
    hfake->SetStats(0);
    hfake->SetTitle(Form("Layer-%s (%i Staves with DeadPixels), from run%s to run%s",layer.c_str(),(int)ceil(staveswithdead[stoi(layer)]/nRuns),runNumbers.back().c_str(),runNumbers[0].c_str()));
    hfake->Draw();
    if(stoi(layer)<=2){
      for(int istave=0; istave<myAnalysis.stavesInLayer(stoi(layer)); istave++){
        leg->AddEntry(trend[stoi(layer)][istave], Form("Stv%d",istave), "p");
        trend[stoi(layer)][istave]->Draw("P same");
      }
      leg->Draw("same");
      canvas->SaveAs(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      canvas->SaveAs(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      delete canvas;
      delete leg;
    }
    else{
      leg->SetNColumns(2);
      hfake->SetTitle(Form("Layer-%s Upper (%i Staves with DeadPixels), from run%s to run%s",layer.c_str(),(int)ceil(staveswithdead[stoi(layer)]/nRuns),runNumbers.back().c_str(),runNumbers[0].c_str()));
      for(int istave=0; istave<myAnalysis.stavesInLayer(stoi(layer))/2; istave++){
        leg->AddEntry(trend[stoi(layer)][istave], Form("Stv%d",istave), "p");
        trend[stoi(layer)][istave]->Draw("P same");
      }
      leg->Draw("same");
      canvas->SaveAs(Form("../Plots/Layer%sU_DeadPixels_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      canvas->SaveAs(Form("../Plots/Layer%sU_DeadPixels_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));

      TLegend *leg2 = new TLegend(0.904, 0.197,0.997,0.898);
      leg2->SetNColumns(2);
      auto *hfake2 = (TH1F*)hfake->Clone();
      hfake2->SetTitle(Form("Layer-%s Lower (%i Staves with DeadPixels), from run%s to run%s",layer.c_str(),(int)ceil(staveswithdead[stoi(layer)]/nRuns),runNumbers.back().c_str(),runNumbers[0].c_str()));
      hfake2->Draw();
      for(int istave=myAnalysis.stavesInLayer(stoi(layer))/2; istave<myAnalysis.stavesInLayer(stoi(layer)); istave++){
        leg2->AddEntry(trend[stoi(layer)][istave], Form("Stv%d",istave), "p");
        trend[stoi(layer)][istave]->Draw("P same");
      }
      leg2->Draw("same");
      canvas->SaveAs(Form("../Plots/Layer%sL_DeadPixels_run%s-run%s.pdf", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      canvas->SaveAs(Form("../Plots/Layer%sL_DeadPixels_run%s-run%s.root", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      delete canvas;
      delete leg;
      delete leg2;
    }
  }

  gStyle->SetPalette(1);
  for (string layer : laynums){ // loop over layers
    TCanvas *canvas = new TCanvas();
    canvas->cd();
    canvas->SetTickx();
    canvas->SetTicky();
    canvas->SetLogz();
    canvas->SetRightMargin(0.15);
    auto histos = myAnalysis.loadLayer(stoi(layer));
    for(int i = histos.size()-1; i >= 0; i--){ //loop over number of histograms
      auto hist = histos[i];
      hist->Draw("colz");
      hist->SetMinimum(8e-1);
      hist->SetMaximum(maxtot+0.5*maxtot);
      hist->SetTitle(Form("Layer%s Run%s (%i/%i);Chip Number; Stave Number",layer.c_str(),myAnalysis.getRunNumber(hist).c_str(),(int)histos.size()-i,(int)histos.size()));
      hist->GetZaxis()->SetTitle("# Dead Pixels");
      // Save frames of GIF
      canvas->Print(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.gif+40", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
      // Save last frame as gif++ so gif loops
      if (i==0)canvas->Print(Form("../Plots/Layer%s_DeadPixels_run%s-run%s.gif++", layer.c_str(),runNumbers.back().c_str(),runNumbers[0].c_str()));
    }
  }

} // end of analyzeLayerDeadPixels()