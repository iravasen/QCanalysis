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



struct LayerParameters
{
	string name;
	int nFirstStave;
	int nLastStave;
	int nChipsPerStave;
	int nStavesForLayer;
};

//struct LayerParameters l0 = { "IB", 1, 12, 9, 1 };
//struct LayerParameters l1 = { "IB", 13, 28, 9, 13 };
//struct LayerParameters l2 = { "IB", 29, 48, 9, 29 };
//struct LayerParameters l3Top = { "ML", 1, 12, 112, 1 };
//struct LayerParameters l3Bottom = { "ML", 13, 24, 112, 1 };
//struct LayerParameters l4Top = { "ML", 25, 39, 112, 25 };
//struct LayerParameters l4Bottom = { "ML", 40, 54, 112, 25 };
//struct LayerParameters l5Top = { "OL", 1, 21, 196, 1 };
//struct LayerParameters l5Bottom = { "OL", 22, 42, 196, 1 };
//struct LayerParameters l6Top = { "OL", 43, 66, 196, 43 };
//struct LayerParameters l6Bottom = { "OL", 67, 90, 196, 43 };



std::pair <Double_t, Double_t> GetAvgThresholdsForStaveAndNumberOfFailedNChip(TH2F* h, Int_t stvNum, Int_t numChip) {
	std::pair <Double_t, Double_t> avgFailed; 
	TH1D* hproj = h->ProjectionX("proj", stvNum, stvNum);
	Double_t sum = 0;
	Int_t cNumChip = numChip;
	Int_t failedNumChip = 0;
	for (int i = 1; i <= numChip; i++) {
		Double_t forChip = hproj->GetBinContent(i);
		if (forChip == 0) {
			cNumChip--;
			failedNumChip++;
		}
		sum += forChip;
	}
	Double_t avg = 0;
	if (cNumChip > 0) {
		avg = sum / cNumChip;
	}
	avgFailed.first = avg;
	avgFailed.second = failedNumChip;
	return avgFailed;
}

//Double_t GetNumberOfFailedNChip(TH2F* h, Int_t stvNum, Int_t numChip) {
//	TH1D* hproj = h->ProjectionX("proj", stvNum, stvNum);
//	Int_t failedNumChip = 0;
//	for (int i = 1; i <= numChip; i++) {
//		if (hproj->GetBinContent(i) == 0) {
//			failedNumChip++;
//		}
//	}
//	return failedNumChip;
//}


set <string> RemoveEmptyRuns(set <string> allRunNumbers, set <string> toDelete) {
	for (auto run : toDelete) {
		allRunNumbers.erase(run);
	}
	return allRunNumbers;
}



TGraph* AddStyle(TGraph* result, Int_t color, Int_t style) {
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

TGraph* AddLine(int numLayer, Double_t beginX, Double_t endX, Double_t beginY, Double_t endY) {
	if (numLayer >= 3) {
		Double_t plotX2[2] = { beginX, endX };
		Double_t plotY2[2] = { beginY, endY };
		TGraph* graph = new TGraph(2, plotX2, plotY2);
		return graph;
	}
	return new TGraph(); 
}


void PrintGraphToPDF(int numLayer, TCanvas* canvas, TMultiGraph* graph, TLegend* legend, string nameFile) {
	canvas->SetBottomMargin(0.19);
	canvas->SetTopMargin(0.075);
	if (numLayer == 0) {
		canvas->Print(Form("%s.pdf[", nameFile.c_str()));
	}
	graph->Draw("AP");
	legend->Draw();
	canvas->Print(Form("%s.pdf", nameFile.c_str()));
	canvas->Print(Form("./analysismacros/plot/%s_Layer%i.png", nameFile.c_str(), numLayer));
	if (numLayer == 10) {
		canvas->Print(Form("%s.pdf]", nameFile.c_str()));
	}
}

vector <LayerParameters> init() {
	LayerParameters l0 = { "IB", 1, 12, 9, 1 };
	LayerParameters l1 = { "IB", 13, 28, 9, 13 };
	LayerParameters l2 = { "IB", 29, 48, 9, 29 };
	LayerParameters l3Top = { "ML", 1, 12, 112, 1 };
	LayerParameters l3Bottom = { "ML", 13, 24, 112, 1 };
	LayerParameters l4Top = { "ML", 25, 39, 112, 25 };
	LayerParameters l4Bottom = { "ML", 40, 54, 112, 25 };
	LayerParameters l5Top = { "OL", 1, 21, 196, 1 };
	LayerParameters l5Bottom = { "OL", 22, 42, 196, 1 };
	LayerParameters l6Top = { "OL", 43, 66, 196, 43 };
	LayerParameters l6Bottom = { "OL", 67, 90, 196, 43 };

	vector <LayerParameters> config = { l0, l1, l2, l3Top, l3Bottom, l4Top, l4Bottom, l5Top, l5Bottom, l6Top, l6Bottom };

	return config;
}



void AnalyzeAllLayerThresholds() {
	
	/*vector <string> l0 = { "IB", "1", "12", "9", "1" };
	vector <string> l1 = { "IB", "13", "28", "9", "13" };
	vector <string> l2 = { "IB", "29", "48", "9", "29" };
	vector <string> l3Top = { "ML", "1", "12", "112", "1" };
	vector <string> l3Bottom = { "ML", "13", "24", "112", "1" };
	vector <string> l4Top = { "ML", "25", "39", "112", "25" };
	vector <string> l4Bottom = { "ML", "40", "54", "112", "25" };
	vector <string> l5Top = { "OL", "1", "21", "196", "1" };
	vector <string> l5Bottom = { "OL", "22", "42", "196", "1" };
	vector <string> l6Top = { "OL", "43", "66", "196", "43" };
	vector <string> l6Bottom = { "OL", "67", "90", "196", "43" };*/
	vector <LayerParameters> config = init();

	vector <Int_t> colorConfig = { kRed, kGreen, kBlue, kMagenta, kOrange };
	vector <Int_t> styleConfig = { 24,25,26,30,28,32 };
	string filepath;
	cout << endl << "Please enter the full path to the file to analyze" << endl;
	gSystem->Exec("ls ../Data/*all-layers* -Art | tail -n 500");
	cin >> filepath;
	//"/mnt/c/Users/User/source/repos/QCanalysis/Data/Output_all-layers_THRMAPS_DEADPIXMAPS_from_run543470_to_run545333.root"
	TFile* analyzefile = new TFile(filepath.c_str());
	set <string> allRunNumbers = GetRunNumbers(analyzefile);
	for (int numLayer = 0; numLayer <= 10; numLayer++) {
		LayerParameters configToUse = config[numLayer];
		TMultiGraph* avgLayerThreshold = new TMultiGraph();
		TMultiGraph* failedChipsForLayer = new TMultiGraph();
		TLegend* legend = new TLegend(1, 0.19, 0.9, 0.925);
		TLegend* legendForChips = new TLegend(1, 0.19, 0.9, 0.925);
		int colorPointer = 0;
		int stylePointer = 0;
		set <string> emptyRunNumbers = {};
		for (int currentStv = configToUse.nFirstStave; currentStv <= configToUse.nLastStave; currentStv++) {
			TGraph* stvGraph = new TGraph();
			TGraph* failedChips = new TGraph();
			Double_t runSequenceNumber = 1;
			for (auto runNumber : allRunNumbers) {
				TList* listKeys = analyzefile->GetListOfKeys();
				for (TObject* key : *listKeys) {
					string fileName = key->GetName();
					if (fileName.find(configToUse.name) != string::npos && fileName.find(runNumber) != string::npos) {
						TH2F* h = (TH2F*)analyzefile->Get(key->GetName());
						if (h->Integral() == 0) {
							emptyRunNumbers.insert(runNumber);
							continue;
						}
						std::pair <Double_t, Double_t> avgFailed = GetAvgThresholdsForStaveAndNumberOfFailedNChip(h, currentStv, configToUse.nChipsPerStave);
						Double_t avgThrForStv = avgFailed.first;
						Int_t failedNumChip = avgFailed.second;
						if (avgThrForStv != 0) {
							stvGraph->AddPoint(runSequenceNumber, avgThrForStv);
						}
						failedChips->AddPoint(runSequenceNumber, failedNumChip); 
						runSequenceNumber++;
					}
				}
			}
			stvGraph = AddStyle(stvGraph, colorConfig[colorPointer], styleConfig[stylePointer]);
			failedChips = AddStyle(failedChips, colorConfig[colorPointer], styleConfig[stylePointer]);
			if (colorPointer == colorConfig.size() - 1) {
				colorPointer = 0;
				stylePointer = stylePointer + 1;
			}
			else {
				colorPointer = colorPointer + 1;
			}
			legend->AddEntry(stvGraph, Form("Stv%d", currentStv - configToUse.nStavesForLayer), "p");
			avgLayerThreshold->Add(stvGraph);
			legendForChips->AddEntry(failedChips, Form("Stv%d", currentStv - configToUse.nStavesForLayer), "p");
			failedChipsForLayer->Add(failedChips);
		}
		string layerName = (numLayer < 3) ? to_string(numLayer) : to_string((numLayer + 1) / 2 + 1) + ((numLayer % 2 == 0) ? "Bottom" : "Top");

		avgLayerThreshold->SetTitle(Form("Average thresholds for stave at layer-%s", layerName.c_str()));
		TAxis* X = avgLayerThreshold->GetXaxis();

		failedChipsForLayer->SetTitle(Form("Failed chips for stave at layer-%s", layerName.c_str()));
		TAxis* XChips = failedChipsForLayer->GetXaxis();
		set <string> notEmptyRuns = RemoveEmptyRuns(allRunNumbers, emptyRunNumbers);
		XChips->SetNdivisions((notEmptyRuns.size() * 2));
		X->SetNdivisions((notEmptyRuns.size() * 2));
		Double_t runSequenceNumber = 2;
		X->ChangeLabel(1, 315, 0, 10, -1, -1, "");
		X->ChangeLabel(notEmptyRuns.size() + 2, 315, 0, 10, -1, -1, "");
		XChips->ChangeLabel(1, 315, 0, 10, -1, -1, "");
		XChips->ChangeLabel(notEmptyRuns.size() + 2, 315, 0, 10, -1, -1, "");
		for (auto runNumber : notEmptyRuns) {
			X->ChangeLabel(runSequenceNumber, 300, -1, 10, -1, -1, Form("run%i", stoi(runNumber)));
			XChips->ChangeLabel(runSequenceNumber, 300, -1, 10, -1, -1, Form("run%i", stoi(runNumber)));
			runSequenceNumber++;
		}
		avgLayerThreshold->GetYaxis()->SetTitle("Average thresholds for stave");
		failedChipsForLayer->GetYaxis()->SetTitle("Failed chips for stave");
		string nameFileAvg = "avgThreshold";
		string nameFileChips = "Failed_chips";
		TCanvas* canvas = new TCanvas(Form("c_L%s", layerName.c_str()), Form("c_L%s", layerName.c_str()), 0, 0, 640, 480);
		PrintGraphToPDF(numLayer, canvas, avgLayerThreshold, legend, nameFileAvg);
		canvas->Close();
		TCanvas* canvasForChips = new TCanvas(Form("c2_L%s", layerName.c_str()), Form("c2_L%s", layerName.c_str()), 0, 0, 640, 480);
		TGraph* line = AddLine(numLayer, -1.0, notEmptyRuns.size() * 2, configToUse.nChipsPerStave / 4.0, configToUse.nChipsPerStave / 4.0);
		failedChipsForLayer->Add(line, "L");
		PrintGraphToPDF(numLayer, canvasForChips, failedChipsForLayer, legendForChips, nameFileChips);
		canvasForChips->Close();
	}
}