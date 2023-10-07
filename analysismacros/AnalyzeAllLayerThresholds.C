//#include "inc/itsAnalysis.hh"
//#include "QualityControl/PostProcessingInterface.h"
//#include "QualityControl/Reductor.h"
//#include "QualityControl/DatabaseFactory.h"
//#include "QualityControl/RootClassFactory.h"
//#include "QualityControl/DatabaseInterface.h"
//#include "QualityControl/MonitorObject.h"
//#include "QualityControl/QcInfoLogger.h"
//#include "QualityControl/CcdbDatabase.h"
//#include "inc/ccdb.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "vector"
#include "string"
#include "set"
//include "stdio.h"



//using namespace o2::quality_control::repository;
//using namespace o2::quality_control::core;
//using namespace std;




Double_t GetAvgThresholdsForStave(TH2F* h, Int_t stvNum, Int_t numChip) {
	TH1D* hproj = h->ProjectionX("proj", stvNum, stvNum);
	Double_t sum = 0;
	Int_t cNumChip = numChip;
	for (int i = 1; i <= numChip; i++) {
		Double_t forChip = hproj->GetBinContent(i);
		if (forChip == 0) {
			cNumChip--;
		}
		sum += forChip;
	}
	Double_t avg = 0;
	if (cNumChip > 0) {
		avg = sum / cNumChip;
	}
	return avg;
}

Double_t GetNumberOfDisabledNChip(TH2F* h, Int_t stvNum, Int_t numChip) {
	TH1D* hproj = h->ProjectionX("proj", stvNum, stvNum);
	Int_t disabledNumChip = 0;
	for (int i = 1; i <= numChip; i++) {
		if (hproj->GetBinContent(i) == 0) {
			disabledNumChip++;
		}
	}
	return disabledNumChip;
}


set <string> RemoveEmptyRuns(set <string> allRunNumbers, set <string> toDelete) {
	for (auto run : toDelete) {
		allRunNumbers.erase(run);
	}
	return allRunNumbers;
}


TGraph* AddPointAndStyle(TGraph* result, Double_t x, Double_t y, Int_t color, Int_t style){
	result->AddPoint(x, y);
	result->SetMarkerColor(color);
	result->SetMarkerStyle(style);
	result->SetMarkerSize(1.5);
	return result;
}




set <string> GetRunNumbers(TFile* analyzefile)
{
	set <string> result;
	TList* listKeys = analyzefile->GetListOfKeys();
	for (TObject* key : *listKeys) {
		string fileName = key->GetName();
		string sRunNumber = fileName.substr(fileName.find("run") + 3, 6);
		result.insert(sRunNumber);
	}
	return result;
}

void PrintGraphToPDF(int numLayer, TCanvas* canvas, TMultiGraph* graph, TLegend* legend, string nameFile) {
	canvas->SetBottomMargin(0.135);
	if (numLayer == 0) {
		canvas->Print(Form("%s[", nameFile.c_str()));
	}
	graph->Draw("AP");
	legend->Draw();
	canvas->Print(Form("%s", nameFile.c_str()));
	if (numLayer == 10) {
		canvas->Print(Form("%s]", nameFile.c_str()));
	}
}



void AnalyzeAllLayerThresholds()
{
	vector <string> l0 = { "IB", "1", "12", "9", "1" };
	vector <string> l1 = { "IB", "13", "28", "9", "13" };
	vector <string> l2 = { "IB", "29", "48", "9", "29" };
	vector <string> l3Top = { "ML", "1", "12", "112", "1" };
	vector <string> l3Bottom = { "ML", "13", "24", "112", "1" };
	vector <string> l4Top = { "ML", "25", "39", "112", "25" };
	vector <string> l4Bottom = { "ML", "40", "54", "112", "25" };
	vector <string> l5Top = { "OL", "1", "21", "196", "1" };
	vector <string> l5Bottom = { "OL", "22", "42", "196", "1" };
	vector <string> l6Top = { "OL", "43", "66", "196", "43" };
	vector <string> l6Bottom = { "OL", "67", "90", "196", "43" };
	vector <vector<string>> config = { l0, l1, l2, l3Top, l3Bottom, l4Top, l4Bottom, l5Top, l5Bottom, l6Top, l6Bottom };
	vector <Int_t> colorConfig = { kRed, kGreen, kBlue, kMagenta, kOrange };
	vector <Int_t> styleConfig = { 24,25,26,30,28,32 };
	TFile* analyzefile = new TFile("/mnt/c/Users/User/source/repos/QCanalysis/Data/Output_all-layers_THRMAPS_DEADPIXMAPS_from_run542392_to_run543954.root");
	set <string> allRunNumbers = GetRunNumbers(analyzefile);
	for (int numLayer = 0; numLayer <= 10; numLayer++) {
		vector<string> configToUse = config[numLayer];
		TMultiGraph* avgLayerThreshold = new TMultiGraph();
		TMultiGraph* disabledChipsForLayer = new TMultiGraph();
		TLegend* legend = new TLegend(1, 0.135, 0.9, 0.9);
		TLegend* legendForChips = new TLegend(1, 0.135, 0.9, 0.9);
		int colorPointer = 0;
		int stylePointer = 0;
		set <string> emptyRunNumbers = {};
		for (int currentStv = stoi(configToUse[1]); currentStv <= stoi(configToUse[2]); currentStv++) {
			TGraph* stvGraph = new TGraph();
			TGraph* disabledChips = new TGraph();
			Double_t runSequenceNumber = 1;
			for (auto runNumber : allRunNumbers) {
				TList* listKeys = analyzefile->GetListOfKeys();
				for (TObject* key : *listKeys) {
					string fileName = key->GetName();
					if (fileName.find(configToUse[0]) != string::npos && fileName.find(runNumber) != string::npos) {
						TH2F* h = (TH2F*)analyzefile->Get(key->GetName());
						if (h->Integral() == 0) {
							emptyRunNumbers.insert(runNumber);
							continue;
						}
						Double_t avgThrForStv = GetAvgThresholdsForStave(h, currentStv, stoi(configToUse[3]));
						Int_t disabledNumChip = GetNumberOfDisabledNChip(h, currentStv, stoi(configToUse[3]));
						stvGraph = AddPointAndStyle(stvGraph, runSequenceNumber, avgThrForStv, colorConfig[colorPointer], styleConfig[stylePointer]);
						//if (disabledNumChip != 0) {
						disabledChips = AddPointAndStyle(disabledChips, runSequenceNumber, disabledNumChip, colorConfig[colorPointer], styleConfig[stylePointer]);
						//}
						runSequenceNumber++;
					}
				}
			}
			if (colorPointer == colorConfig.size() - 1) {
				colorPointer = 0;
				stylePointer = stylePointer + 1;
			}
			else {
				colorPointer = colorPointer + 1;
			}
			legend->AddEntry(stvGraph, Form("Stv%d", currentStv - stoi(configToUse[4])), "p");
			avgLayerThreshold->Add(stvGraph);
			legendForChips->AddEntry(disabledChips, Form("Stv%d", currentStv - stoi(configToUse[4])), "p");
			disabledChipsForLayer->Add(disabledChips);
		}
		string layerName = (numLayer < 3) ? to_string(numLayer) : to_string((numLayer + 1) / 2 + 1) + ((numLayer % 2 == 0) ? "Bottom" : "Top");

		avgLayerThreshold->SetTitle(Form("Average thresholds for stave at layer-%s", layerName.c_str()));
		TAxis* X = avgLayerThreshold->GetXaxis();

		disabledChipsForLayer->SetTitle(Form("Disabled chips for stave at layer-%s", layerName.c_str()));
		TAxis* XChips = disabledChipsForLayer->GetXaxis();
		set <string> notEmptyRuns = RemoveEmptyRuns(allRunNumbers, emptyRunNumbers);
		XChips->SetNdivisions((notEmptyRuns.size() * 2));
		X->SetNdivisions((notEmptyRuns.size() * 2));
		Double_t runSequenceNumber = 2;
		X->ChangeLabel(1, 325, 0, 10, -1, -1, "");
		X->ChangeLabel(notEmptyRuns.size() + 2, 325, 0, 10, -1, -1, "");
		XChips->ChangeLabel(1, 325, 0, 10, -1, -1, "");
		XChips->ChangeLabel(notEmptyRuns.size() + 2, 325, 0, 10, -1, -1, "");
		for (auto runNumber : notEmptyRuns) {
			X->ChangeLabel(runSequenceNumber, 325, -1, 10, -1, -1, Form("run%i", stoi(runNumber)));
			XChips->ChangeLabel(runSequenceNumber, 325, -1, 10, -1, -1, Form("run%i", stoi(runNumber)));
			runSequenceNumber++;
		}
		avgLayerThreshold->GetYaxis()->SetTitle("Average thresholds for stave");
		disabledChipsForLayer->GetYaxis()->SetTitle("Disabled chips for stave");
		string nameFileAvg = "avgThreshold.pdf";
		string nameFileChips = "disabled_chips.pdf";
		TCanvas* canvas = new TCanvas(Form("c_L%s", layerName.c_str()), Form("c_L%s", layerName.c_str()), 0, 0, 640, 480);
		PrintGraphToPDF(numLayer, canvas, avgLayerThreshold, legend, nameFileAvg);
		canvas->Close();
		TCanvas* canvasForChips = new TCanvas(Form("c2_L%s", layerName.c_str()), Form("c2_L%s", layerName.c_str()), 0, 0, 640, 480);
		PrintGraphToPDF(numLayer, canvasForChips, disabledChipsForLayer, legendForChips, nameFileChips);
		canvasForChips->Close();
	}
}