#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TMath.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLine.h"
#include "TGraph.h"

#include <math.h>
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <vector>

//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
//#include "classes/DelphesClasses.h"
//#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#include "external/ExRootAnalysis/ExRootResult.h"
//#endif

#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip> 

std::string path="13_03_18/";
bool m_DEBUG=true;

float HTCutValue = -1.0;
std::vector<float> HTCut;

Color_t ci[] = {kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack,kWhite,kBlack};
Color_t cf[] = {kMagenta+3,kGreen-3,kGreen+5,kCyan,kMagenta-3,kGreen,kBlue,kRed,kYellow,kBlack,kBlue-3, kBlue+3, kMagenta+3,kGreen-3,kGreen+5,kCyan,kMagenta-3,kGreen,kBlue,kRed,kYellow,kBlack,kBlue-3, kBlue+3, kMagenta+3};

/*
root -l PlotMET_BKG.C\(\"delphes_output.root\"\)
*/

TFile *f;

std::vector<float> fillHTCutFlowMap(std::vector<std::string> dataset, std::string hist_name, std::vector<std::string> process, float HTCutVal)
{
   TH1F* h[50];
   std::stringstream temp_name;
   std::vector<float> counts;

   for(Int_t i = 0; i< dataset.size(); i++){
   //for(Int_t i = dataset.size()-1; i >= 0; --i){
     if(m_DEBUG)std::cout << "i: " << i << std::endl;
     temp_name << hist_name << "_" << process.at(i);
     /*if(m_DEBUG)*/std::cout << "temp_name: " << temp_name.str() << std::endl;
     TFile *f = TFile::Open(dataset.at(i).c_str(),"READ");
     TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
     int thisBin=htemp->FindBin(HTCutVal);
     float integral=htemp->Integral(thisBin,9999);
     std::cout << "This bin: " << thisBin << " HTCutValue: " << HTCutVal << " Integral: " << htemp->Integral(thisBin,9999) << std::endl;
     counts.push_back(integral);
     temp_name.str("");
   }

   return counts;
}

void writeCutFlow(std::vector<std::string> dataset, std::string hist_name, std::vector<std::string> process,std::map<std::string,std::string> label_map, std::string channel, std::map<int,std::vector<float> > counts)
{

   std::cout << "On write cut flow function" << std::endl;

   float TotBkG[20]={0};

   Double_t x1[9],y1[HTCut.size()][9],x2[9],y2[HTCut.size()][9];
   Int_t n = 9, k1 = 0, k2 = 0;

   TCanvas *c1 = new TCanvas("c", "", 800, 600);
   c1->SetLogy();

   std::ofstream myfile;
   std::string filename=path+"weighted_cutflow_"+channel+".tex";
   if(hist_name.find("Weighted") != std::string::npos) myfile.open (filename.c_str());

   std::string latex_name;
   if(channel=="emu"){
     latex_name="$e\\mu$";
   }else if(channel=="mumu"){
     latex_name="$\\mu\\mu$";
   }else if(channel=="ee"){
     latex_name="ee";
   }


   //if(channel=="emu") myfile << "\\hspace{-3cm} \n";

   myfile << "{\\tiny \n";
   myfile << "\\begin{tabular}{| c | c | c | c | c";
   if(channel=="ee") myfile << " | c";
   for(UInt_t icut=0;icut<HTCut.size();icut++) myfile << " | c ";
   myfile << "|} \n";
   myfile << "\\hline \n";
   myfile << "$m_{\\Phi}$ TeV (mixing) & $\\sigma_{MG}\\times \\mathcal{L}_{int}$ & " << latex_name << " & $\\ell^{\\pm}\\ell^{\\pm}$ & $\\geq 4~$ jets";
   if(channel=="ee") myfile << "& $m_{ee}>110$ GeV";
   for(UInt_t icut=0;icut<HTCut.size();icut++) myfile << " & " << HTCut.at(icut)/1000. << "TeV ";
   myfile << "\\\\ \n";
   myfile << "\\hline \n";

   /*if(channel!="ee") myfile << "\\begin{tabular}{| c | c | c | c | c | c |} \n";
   else myfile << "\\begin{tabular}{| c | c | c | c | c | c | c |} \n";
   myfile << "\\hline \n";
   if(channel!="ee") myfile << "\\multicolumn{6}{|c|}{"<< latex_name << "} \\\\ \n";
   else myfile << "\\multicolumn{7}{|c|}{"<< latex_name << "} \\\\ \n";
   myfile << "\\hline \n";
   if(channel!="ee") myfile << "$m_{\\Phi}$ TeV (channel) & no req. & " << latex_name << " & SS & $>3~$ jets & $H_{T}>" << HTCutValue << "$~GeV \\\\ \n";
   else myfile << "$m_{\\Phi}$ TeV (channel) & no req. & $ee$ & SS & $>3~$ jets & $m_{ee}>110 GeV$ & $H_{T}>" << HTCutValue << "$~GeV \\\\ \n";
   myfile << "\\hline \n";*/

   std::string postfix="", postfix2="";

   int writeIt=0;
   for(Int_t i = dataset.size()-1; i >= 0; --i){
     //std::cout << "dataset: " << dataset.at(i) << "counts: " << counts.at(i) << std::endl;
     if(process.at(i).find("signal1.5") == std::string::npos &&
         process.at(i).find("signal2.0") == std::string::npos &&
         process.at(i).find("signal2.4PSI2.5") == std::string::npos && process.at(i).find("signal") != std::string::npos) continue;
     TFile *f = TFile::Open(dataset.at(i).c_str(),"READ");
     TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
     if(writeIt==0 && (process.at(i).find("signal") == std::string::npos)){ std::cout << "\\hline" << std::endl; myfile << "\\hline \n";}
     std::cout << label_map[process.at(i)] << " & "
               << htemp->GetBinContent(1) << " & ";
     float normFactor=2.4/2.02, normFactor2=0.5;
     if(process.at(i).find("signal") == std::string::npos) {normFactor=1.0;normFactor2=1.0;}
     myfile << label_map[process.at(i)] << " & "
               << std::fixed << std::setprecision(2) << normFactor*htemp->GetBinContent(1) << " & ";
     if(channel=="ee"){ std::cout << normFactor*htemp->GetBinContent(2) << " & "; myfile << normFactor*htemp->GetBinContent(2) << " & ";}
     if(channel=="emu"){ std::cout << normFactor*htemp->GetBinContent(3) << " & "; myfile << normFactor*htemp->GetBinContent(3) << " & ";}
     if(channel=="mumu"){ std::cout << normFactor*htemp->GetBinContent(4) << " & "; myfile << normFactor*htemp->GetBinContent(4) << " & ";}
     std::cout << htemp->GetBinContent(5);
     myfile << normFactor*normFactor2*htemp->GetBinContent(5);
     if(channel == "ee"){
       std::cout << " & " << htemp->GetBinContent(6) << " & " << htemp->GetBinContent(7);
       myfile << " & " << normFactor*normFactor2*htemp->GetBinContent(6) << " & " << normFactor*normFactor2*htemp->GetBinContent(7);
     }else{
       std::cout << " & "<< htemp->GetBinContent(6);
       myfile << " & "<< normFactor*normFactor2*htemp->GetBinContent(6);
     }
     for(UInt_t icut=0;icut<HTCut.size();icut++){
       std::cout << " & " << counts[icut].at(i);
       myfile << " & " << normFactor*normFactor2*counts[icut].at(i);
     }
     std::cout << " \\\\" << std::endl;
     myfile << " \\\\ \n";
     if(process.at(i).find("signal") != std::string::npos ){
        if(process.at(i).find("emu") != std::string::npos){
          postfix="emu";
          for(UInt_t icut=0;icut<HTCut.size();icut++){
             if(m_DEBUG) std::cout << "counts1: " << counts[icut].at(i) << " k1: " << k1 << std::endl;
             y1[icut][k1] = 3./counts[icut].at(i);
          }
          std::string label = label_map[process.at(i)];
          //std::string label = "1.5 blah";
          //std::cout << label << std::endl; 
          std::string masa = label.substr(0,3);
          //std::string masa = "10.0";
          float mass = (float)std::atof (masa.c_str());
          std::cout << mass << std::endl;
          x1[k1] = mass;
          k1++;
        } else{
          postfix2=channel;
          for(UInt_t icut=0;icut<HTCut.size();icut++){
             if(m_DEBUG) std::cout << "counts2: "<< counts[icut].at(i) << " k2: " << k2 << std::endl;
             y2[icut][k2] = 3./counts[icut].at(i);
          }
          std::string label = label_map[process.at(i)];
          //std::string label = "1.5 blah";
          //std::cout << label << std::endl; 
          std::string masa = label.substr(0,3);
          //std::string masa = "10.0";
          float mass = (float)std::atof (masa.c_str());
          std::cout << mass << std::endl;
          x2[k2] = mass;
          k2++;
        }
        //continue;
     }else{
      writeIt++;
      TotBkG[0]+=htemp->GetBinContent(1);
      TotBkG[1]+=htemp->GetBinContent(2);
      TotBkG[2]+=htemp->GetBinContent(3);
      TotBkG[3]+=htemp->GetBinContent(4);
      TotBkG[4]+=htemp->GetBinContent(5);
      TotBkG[5]+=htemp->GetBinContent(6);
      if(channel == "ee") {
        TotBkG[6]+=htemp->GetBinContent(7);
        for(UInt_t icut=0;icut<HTCut.size();icut++){
           TotBkG[7+icut]+=counts[icut].at(i);
        }
      }else{
        for(UInt_t icut=0;icut<HTCut.size();icut++){
           TotBkG[6+icut]+=counts[icut].at(i);
           std::cout << "Checkthis: " << counts[icut].at(i) << " " << TotBkG[6+icut] << std::endl;
        }
      }
     }
   }
   std::cout << "\\hline" << std::endl;
   std::cout << "Total Background " << " & "
   << TotBkG[0]  << " & ";
   myfile << "\\hline \n";
   myfile << "Total Background " << " & "
   << TotBkG[0]  << " & ";
   if(channel=="ee"){ std::cout << TotBkG[1]  << " & "; myfile << TotBkG[1]  << " & ";}
   if(channel=="emu"){ std::cout << TotBkG[2]  << " & "; myfile << TotBkG[2]  << " & ";}
   if(channel=="mumu"){ std::cout << TotBkG[3]  << " & "; myfile << TotBkG[3]  << " & ";}
   std::cout << TotBkG[4]  << " & " << TotBkG[5];
   myfile << TotBkG[4]  << " & " << TotBkG[5];

   if(channel == "ee"){
     std::cout << " & " << TotBkG[6];
     myfile << " & " << TotBkG[6];
   }


   for(UInt_t icut=0;icut<HTCut.size();icut++){
      if(channel == "ee"){
        std::cout << " & " << TotBkG[7+icut] << std::endl;
        myfile << " & " << TotBkG[7+icut];
      }else{
        std::cout << " & " << TotBkG[6+icut] << std::endl;
        myfile << " & " << TotBkG[6+icut];;
      }
   }

   myfile << " \\\\ \n";
   myfile << "\\hline \n";
   myfile << "\\end{tabular} \n";
   myfile << "} \n";

   myfile.close();

   std::stringstream graph_name;

   for(UInt_t icut=0;icut<HTCut.size();icut++){

      double y3[9]={0},y4[9]={0};
      int ncounts=sizeof y1[0]/sizeof y1[0][0];
      int ncounts2=sizeof y2[0]/sizeof y2[0][0];
      std::cout << "Array columns size is: " << ncounts << " and is: " << ncounts2 <<std::endl;
      for(UInt_t icounts=0;icounts<ncounts;icounts++){
         y3[icounts]=y1[icut][icounts];
         y4[icounts]=y2[icut][icounts];
      }

      TGraph* gr1 = new TGraph(n,x1,y3);

      std::cout << "x1[0]: " << x1[0] << " y3[0] " << y3[0] << std::endl;
      std::cout << "x1[1]: " << x1[1] << " y3[1] " << y3[1] << std::endl;
      std::cout << "x1[2]: " << x1[2] << " y3[2] " << y3[2] << std::endl;
      std::cout << "x1[3]: " << x1[3] << " y3[3] " << y3[3] << std::endl;
      std::cout << "x1[4]: " << x1[4] << " y3[4] " << y3[4] << std::endl;
      std::cout << "x1[5]: " << x1[5] << " y3[5] " << y3[5] << std::endl;
      std::cout << "x1[6]: " << x1[6] << " y3[6] " << y3[6] << std::endl;
      std::cout << "x1[7]: " << x1[7] << " y3[7] " << y3[7] << std::endl;
      std::cout << "x1[8]: " << x1[8] << " y3[8] " << y3[8] << std::endl;
      //std::cout << "x1[9]: " << x1[9] << " y3[9] " << y3[9] << std::endl;

      std::cout << "x2[0]: " << x2[0] << " y4[0] " << y4[0] << std::endl;
      std::cout << "x2[1]: " << x2[1] << " y4[1] " << y4[1] << std::endl;
      std::cout << "x2[2]: " << x2[2] << " y4[2] " << y4[2] << std::endl;
      std::cout << "x2[3]: " << x2[3] << " y4[3] " << y4[3] << std::endl;
      std::cout << "x2[4]: " << x2[4] << " y4[4] " << y4[4] << std::endl;
      std::cout << "x2[5]: " << x2[5] << " y4[5] " << y4[5] << std::endl;
      std::cout << "x2[6]: " << x2[6] << " y4[6] " << y4[6] << std::endl;
      std::cout << "x2[7]: " << x2[7] << " y4[7] " << y4[7] << std::endl;
      std::cout << "x2[8]: " << x2[8] << " y4[8] " << y4[8] << std::endl;
      //std::cout << "x2[9]: " << x2[9] << " y4[9] " << y4[9] << std::endl;
      //std::cout << "x[10]: " << x[10] << " y[10] " << y[10] << std::endl;
      //std::cout << "x[11]: " << x[11] << " y[11] " << y[11] << std::endl;

      graph_name.str("");
      graph_name << "HT" << HTCut.at(icut);

      gr1->Draw("AC*");
      c1->SaveAs((path+"LimitPlot_"+channel+"_"+postfix+"_"+graph_name.str()+".png").c_str());

      if(channel!="emu"){
        c1->Modified();
        c1->Clear();
        TGraph* gr2 = new TGraph(n,x2,y4);
        gr2->Draw("AC*");
        c1->SaveAs((path+"LimitPlot_"+channel+"_"+postfix2+"_"+graph_name.str()+".png").c_str());
      }

   }   //std::cout << "HTCutValue 2: " << HTCutValue << std::endl;

}

void doStackPlot(std::vector<std::string> dataset, std::string hist_name, std::vector<std::string> process,float scaleF,std::string xLabel, std::map<std::string,std::string> label_map, std::string channel, bool isLog)
{
  if(m_DEBUG) std::cout << "Channel: " << channel << std::endl;

  std::map<std::string,Color_t> color_map;
  color_map["ttbarlep"]=kRed;
  color_map["zjets_mllcut_all"]=kBlue;
  //color_map["ttw"]=kCyan;
  color_map["ttbarw-v2"]=kCyan;
  color_map["ttz"]=kBlue-7;
  color_map["twlep"]=kYellow;
  color_map["ww"]=kGreen+2;
  color_map["wwjj"]=kOrange+3;
  //color_map["zw"]=kViolet;
  color_map["wz-3j-v2"]=kViolet;
  color_map["zz"]=kTeal+9;
  color_map["signal1.5PSI2.5LQee"]=kBlue-10;
  color_map["signal1.6PSI2.5LQee"]=kBlue-9;
  color_map["signal1.7PSI2.5LQee"]=kBlue-8;
  color_map["signal1.8PSI2.5LQee"]=kBlue-7;
  color_map["signal1.9PSI2.5LQee"]=kBlue-6;
  color_map["signal2.0PSI2.5LQee"]=kBlue-5;
  color_map["signal2.1PSI2.5LQee"]=kBlue-4;
  color_map["signal2.2PSI2.5LQee"]=kBlue-3;
  color_map["signal2.3PSI2.5LQee"]=kBlue-2;
  color_map["signal2.4PSI2.5LQee"]=kBlue-1;
  color_map["signal2.5PSI3.0LQee"]=kGreen-4;
  color_map["signal2.6PSI3.0LQee"]=kGreen-3;
  color_map["signal2.7PSI3.0LQee"]=kGreen-2;
  color_map["signal2.8PSI3.0LQee"]=kGreen-1;
  color_map["signal1.5PSI2.5LQmumu"]=kGreen-10;
  color_map["signal1.6PSI2.5LQmumu"]=kGreen-9;
  color_map["signal1.7PSI2.5LQmumu"]=kGreen-8;
  color_map["signal1.8PSI2.5LQmumu"]=kGreen-7;
  color_map["signal1.9PSI2.5LQmumu"]=kGreen-6;
  color_map["signal2.0PSI2.5LQmumu"]=kGreen-5;
  color_map["signal2.1PSI2.5LQmumu"]=kGreen-4;
  color_map["signal2.2PSI2.5LQmumu"]=kGreen-3;
  color_map["signal2.3PSI2.5LQmumu"]=kGreen-2;
  color_map["signal2.4PSI2.5LQmumu"]=kGreen-1;
  color_map["signal2.4PSI3.0LQmumu"]=kBlue-5;
  color_map["signal2.5PSI3.0LQmumu"]=kBlue-4;
  color_map["signal2.6PSI3.0LQmumu"]=kBlue-3;
  color_map["signal2.7PSI3.0LQmumu"]=kBlue-2;
  color_map["signal2.8PSI3.0LQmumu"]=kBlue-1;
  color_map["signal1.5PSI2.5LQemu"]=kMagenta-10;
  color_map["signal1.6PSI2.5LQemu"]=kMagenta-9;
  color_map["signal1.7PSI2.5LQemu"]=kMagenta-8;
  color_map["signal1.8PSI2.5LQemu"]=kMagenta-7;
  color_map["signal1.9PSI2.5LQemu"]=kMagenta-6;
  color_map["signal2.0PSI2.5LQemu"]=kMagenta-5;
  color_map["signal2.1PSI2.5LQemu"]=kMagenta-4;
  color_map["signal2.2PSI2.5LQemu"]=kMagenta-3;
  color_map["signal2.3PSI2.5LQemu"]=kMagenta-2;
  color_map["signal2.4PSI2.5LQemu"]=kMagenta-1;
  color_map["signal2.5PSI3.0LQemu"]=kRed-4;
  color_map["signal2.6PSI3.0LQemu"]=kRed-3;
  color_map["signal2.7PSI3.0LQemu"]=kRed-2;
  color_map["signal2.8PSI3.0LQemu"]=kRed-1;

  THStack *hs = new THStack("hs","");

  TH1F* h[100];
  std::vector<std::string> newname;

  std::stringstream temp_name;
  float maxi = 0;

  gStyle->SetPadRightMargin(0.1);
  TCanvas *c = new TCanvas("c", "", 800, 600);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  TLegend *leg = new TLegend(0.11, 0.70, 0.9, 0.89);
  leg->SetNColumns(3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  int index=0;
  int index1=0;

  std::vector<std::string> labels;
  labels.clear();

  for(Int_t i = dataset.size()-1; i >= 0; --i){
    std::cout << "Now running on dataset: " << dataset.at(i) << " and process: " << process.at(i) << std::endl;
    //for(Int_t i = 0; i < dataset.size(); i++){
    if(m_DEBUG)std::cout << "i: " << i << std::endl;
    // Only stack background histograms
    if(process.at(i).find("signal") != std::string::npos ) continue;
    index1++;
    temp_name << hist_name << "_" << process.at(i);
    if(m_DEBUG)std::cout << "temp_name: " << temp_name.str() << std::endl;
    TFile *f = TFile::Open(dataset.at(i).c_str(),"READ");
    //std::cout << "Seg fault is up here: -1"<< std::endl;
    TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
    //std::cout << "Seg fault is up here: 0"<< std::endl;
    if(htemp->Integral()==0) {temp_name.str("");if(m_DEBUG)std::cout << "Histogram have zero integral" << std::endl; continue;}
    //std::cout << "Seg fault is up here: 1"<< std::endl;
    labels.push_back(label_map[process.at(i)]);
    h[index]=(TH1F*)htemp->Clone(temp_name.str().c_str());
    h[index]->SetLineWidth(2);
    h[index]->SetLineColor(kBlack);
    h[index]->SetFillColor(color_map[process.at(i)]/*rainbow[ki]*/);
    //std::cout << "Seg fault is up here: 2"<< std::endl;
    if(m_DEBUG)std::cout << "Now running on: " << i << " counter " << index << " " << process.at(i) << " Integral: " << htemp->Integral() << std::endl;
    h[index]->SetFillStyle(1001);
    h[index]->SetLineStyle(1/*ls[ki]*/);
    h[index]->SetStats(0);
    //h[i]->SetYTitle("# of events/100 fb^{-1}");
    //h[i]->SetXTitle(xLabel.c_str());
    h[index]->SetTitle("");
    //h[i]->Scale(1./6);
    maxi += h[index]->GetMaximum();
    if(m_DEBUG)std::cout << "This is the maximum: " << maxi << " i " << i << std::endl;
    hs->Add(h[index]);
    if(m_DEBUG)std::cout << "After adding the histogram" << std::endl;
    temp_name.str("");
    index++;
  }

  if(maxi==0) return;

  if(m_DEBUG)std::cout << "After the first loop " << index1 << std::endl;
  int index0=index;

  for(UInt_t k=index1;k<dataset.size();k++){
   if(process.at(k).find("signal") == std::string::npos ) continue;
   temp_name << hist_name << "_" << process.at(k);
   TFile *f = TFile::Open(dataset.at(k).c_str(),"READ");
   if(m_DEBUG)std::cout << "dataset: " << dataset.at(k) << std::endl;
   TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
   if(htemp->Integral()==0) {temp_name.str("");continue;}
   labels.push_back(label_map[process.at(k)]);
   if(m_DEBUG)std::cout << "Now running on: " << k << " counter " << index0 << " " << process.at(k) << " Integral: " << htemp->Integral() << " Name: " << temp_name.str() << std::endl;
   h[index0]=(TH1F*)htemp->Clone(temp_name.str().c_str());
   h[index0]->SetLineStyle(1);
   h[index0]->SetLineColor(color_map[process.at(k)]);
   h[index0]->SetLineWidth(4);
   h[index0]->SetStats(0);
   index0++;
   temp_name.str("");
  }

  if(m_DEBUG)std::cout << "After the second loop on signal samples: " << index << " index0: " << index0 << " labels size " << labels.size() << std::endl;
  for(UInt_t l=0;l<labels.size();l++) if(m_DEBUG)std::cout << "label: " << labels.at(l) << std::endl;

   for (UInt_t ki = 0; ki < index; ++ki) {
    if(m_DEBUG)std::cout << "Adding legends: " << ki << " " << labels.at(ki) << std::endl;
    if(m_DEBUG)std::cout << "Check this ! " << h[ki]->GetName() << " process @" << process.at(ki) << std::endl;
    if(h[ki]->Integral()==0) continue;
    leg->AddEntry(h[ki], labels.at(ki).c_str(), "F");
   }

   for (UInt_t ki = index; ki < index0; ++ki) {
    if(m_DEBUG)std::cout << "Adding legends: " << ki << " " << process.at(ki) << std::endl;
    if(m_DEBUG)std::cout << "Check this ! " << h[ki]->GetName() << std::endl;
    if(h[ki]->Integral()==0) continue;
    if((process.at(ki).find("signal1.5") != std::string::npos ||
       process.at(ki).find("signal2.0") != std::string::npos ||
       process.at(ki).find("signal2.4PSI2.5") != std::string::npos) && channel!="ee"){ leg->AddEntry(h[ki], labels.at(ki).c_str(), "F");}
    else if((process.at(ki).find("signal1.5") != std::string::npos ||
       process.at(ki).find("signal2.0") != std::string::npos ||
       process.at(ki).find("signal2.4") != std::string::npos) && channel=="ee") {
    std::cout << "Found process: " << process.at(ki) << " hist_name " << h[ki-1]->GetName() << " labels: " << labels.at(ki-1) << std::endl;
    leg->AddEntry(h[ki-1], labels.at(ki-1).c_str(), "F");
    }

   }

   TCanvas *cs = new TCanvas("cs","",600,800);
   if(isLog) cs->SetLogy();
   hs->SetMaximum(scaleF*maxi);
   hs->SetMinimum(1e-1);
   //hs->GetYaxis()->SetRangeUser(1e-5,1e+4);
   //hs->SetMinimum(1.0);
   hs->Draw("HIST");
   leg->Draw();

   for(UInt_t j=index;j<index0;j++) {
      /*if(j==index) h[j]->Draw("HIST");
      else*/ 
      h[j]->Scale(0.5*2.4/2.02);
      if((process.at(j).find("signal1.5") != std::string::npos ||
         process.at(j).find("signal2.0") != std::string::npos ||
         process.at(j).find("signal2.4PSI2.5") != std::string::npos) && channel!="ee") {h[j]->Draw("HIST:same");std::cout<<h[j]->GetName()<<std::endl;}
      else if((process.at(j).find("signal1.5") != std::string::npos ||
       process.at(j).find("signal2.0") != std::string::npos ||
       process.at(j).find("signal2.4") != std::string::npos) && channel=="ee") {h[j-1]->Draw("HIST:same");std::cout<<h[j-1]->GetName()<<std::endl;}
      if(m_DEBUG)std::cout << "hist_name: " << h[j]->GetName() << " Integral: " << h[j]->Integral() << " index: " << j << std::endl;
   }

   TPaveText *pt = new TPaveText(0.6,0.6,0.8,0.7,"NDC");
   if(channel=="mumu") pt->AddText("SS #mu#mu");
   else if(channel=="emu") pt->AddText("SS e#mu");
   else if(channel=="ee") pt->AddText("SS ee");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->Draw();

   hs->GetXaxis()->SetTitle(xLabel.c_str());
   hs->GetYaxis()->SetTitleOffset(1.4);
   hs->GetYaxis()->SetTitle("# of events/300 fb^{-1}");
   std::stringstream histname;
   histname << hist_name << "_" << channel << ".png";
   //histname << hist_name << "_nosamescut_bkg.png";

   cs->SaveAs((path+histname.str()).c_str());
   cs->Modified();
   cs->Clear();

   if(hist_name.find("VR_nJets") != std::string::npos){
     f->cd();
     if(channel == "ee"){
       if(m_DEBUG)std::cout << "Here ee" << hist_name << std::endl;
       h[40] =(TH1F*)(hs->GetStack()->Last());
       if(m_DEBUG)std::cout << h[40]->GetName() << std::endl;
       h[40]->Write((hist_name+"_"+channel).c_str());
     }else if(channel == "mumu"){
       if(m_DEBUG)std::cout << "Here mm" << std::endl;
       h[41] =(TH1F*)(hs->GetStack()->Last());
       if(m_DEBUG)std::cout << h[41]->GetName() << std::endl;
       h[41]->Write((hist_name+"_"+channel).c_str());
     }
   }

   std::vector<std::string> testVars;
   testVars.push_back("Heejjjj");testVars.push_back("Hmumujjjj");testVars.push_back("Hemujjjj");
   testVars.push_back("eejjjjMass");testVars.push_back("emujjjjMass");testVars.push_back("mumujjjjMass");
   testVars.push_back("LQMassE");testVars.push_back("LQMassMU");

   bool getSB=false;
   for(UInt_t k=0;k<testVars.size();k++){
      if(hist_name.find(testVars.at(k)) != std::string::npos){
        getSB=true;
        break;
      }
   }

   if(getSB){
     //std::vector<float> HTCut;
     std::vector<float> ZVector;
     //Get background histogram
     h[42] = (TH1F*)(hs->GetStack()->Last());
     TH1F *eeS = (TH1F*)h[42]->Clone("eeS");
     eeS->Reset();
     for(unsigned int j=1;j<eeS->GetXaxis()->GetNbins()+1; j++){
        float content=sqrt(h[42]->GetBinContent(j));
        float bincenter=h[42]->GetBinCenter(j);
        eeS->Fill(bincenter,content);
     }
     //Loop on signals
     if(m_DEBUG)std::cout << "index: " << index << " index0: " << index0 << std::endl;
     for(UInt_t j=index;j<index0;j++) {
        ZVector.clear();
        bool thisFirst=true;
        TH1F *hTemp = (TH1F*)eeS->Clone("hTemp");
        TH1F *hTemp2 = (TH1F*)eeS->Clone("hTemp2");
        TH1F *hTemp3 = (TH1F*)eeS->Clone("hTemp3");
        std::string histname(h[j]->GetName());
        std::string processName;
        for(unsigned int k=0;k<process.size();k++){
          processName=process.at(k);
          if(histname.find(process.at(k))!=std::string::npos) break;
        }
        /*if(m_DEBUG)*/std::cout << "Process: " << processName << " Hist Name: " << h[j]->GetName() << std::endl;
        hTemp3->Reset();
        hTemp2->Reset();
        hTemp->Reset();
        hTemp->Add(h[j]);
        hTemp->Divide(eeS);
        //Z=sqrt(2*(S+B)*ln(1+S/B)-S)
        for(unsigned int k=1;k<hTemp->GetXaxis()->GetNbins()+1; k++){
           float IntS = h[j]->Integral(k,h[j]->GetXaxis()->GetNbins()+1);
           if(IntS<3.0) break;
           float SqIntB = h[42]->Integral(k,h[42]->GetXaxis()->GetNbins()+1);
           if(m_DEBUG)std::cout << "IntS= " << IntS << " BinCenter: " << hTemp->GetBinCenter(k) << std::endl;
           if(m_DEBUG)std::cout << "SqIntB= " << SqIntB << std::endl;
           if(m_DEBUG)std::cout << "IntS/sqrt(SqIntB)= " << IntS/sqrt(SqIntB) << std::endl;
           if(m_DEBUG)std::cout << "intSB: " << IntS/sqrt(SqIntB) << std::endl;
           float Z = sqrt(2*((IntS+SqIntB)*log(1+IntS/SqIntB)-IntS));
           float HTCutVar=hTemp2->GetBinLowEdge(k)+hTemp2->GetBinWidth(k);
           std::cout << "HT Bin: " << HTCutVar << std::endl;
           if(SqIntB==0){
              ZVector.push_back(0);
              continue;
           }
           else ZVector.push_back(Z);
           if(m_DEBUG)std::cout << "full Z: " << Z << std::endl;
           float bincenter = hTemp->GetBinCenter(k);
           hTemp2->Fill(bincenter,IntS/sqrt(SqIntB));
           hTemp3->Fill(bincenter,Z);
           //if(HTCutVar<1200 || HTCutVar>2100) continue;
           //if((histname.find("LQ"+channel)==std::string::npos || channel == "emu") && histname.find("H")!=std::string::npos ) HTCut.push_back(HTCutVar);
           /*if(IntS/sqrt(SqIntB)>5 && thisFirst){
             std::cout << "Found bin : " << hTemp2->FindBin(bincenter) << " S/Sqrt(B): " << IntS/sqrt(SqIntB) << " k: " << k << " Cut: " << hTemp2->GetBinLowEdge(k)+hTemp2->GetBinWidth(k) << " bin width: " << hTemp2->GetBinWidth(k) << std::endl;
             thisFirst = false;
             HTCutVar=hTemp2->GetBinLowEdge(k)+hTemp2->GetBinWidth(k);
             if((histname.find("LQ"+channel)==std::string::npos || channel == "emu") && histname.find("H")!=std::string::npos) HTCut.push_back(HTCutVar);
           }*/
        }
        std::vector<float>::iterator result;
        result = std::max_element(ZVector.begin(), ZVector.end());
        int max_position = std::distance(ZVector.begin(), result);
        float max_z=-1;
        if(ZVector.size()!=0){ max_z=ZVector.at(max_position);
        std::cout << "max element at: " << max_position << " max value " << max_z << std::endl;
        std::cout << "ZVector size: " << ZVector.size() << " Hist size: " << hTemp->GetXaxis()->GetNbins()+1 <<std::endl;
        float HTCutVar2=hTemp3->GetBinLowEdge(max_position+1)+hTemp2->GetBinWidth(max_position+1);
        //if((histname.find("LQ"+channel)==std::string::npos || channel == "emu") && histname.find("H")!=std::string::npos ){
        //  HTCut.push_back(HTCutVar2);
        //}
        std::cout << "HTCutVar2: " << HTCutVar2 << std::endl;
        }
        //std::cout << "max element is: " << result << std::endl;
        //std::cout << "" << ZVector
        hTemp->GetYaxis()->SetTitle("S/#sqrt{B}");
        hTemp->SetTitle(processName.c_str());
        hTemp->Draw();
        //TLine *line = new TLine(300,5,5000,5);
        TLine *line = new TLine(300,max_z,5000,max_z);
        line->SetLineColor(kRed);
        line->SetLineStyle(1);
        line->SetLineWidth(2);
        line->Draw();
        cs->SaveAs((path+string(h[j]->GetName())+"_"+channel+".png").c_str());
        cs->Modified();
        cs->Clear();
        hTemp2->GetYaxis()->SetTitle("S/#sqrt{B}");
        hTemp2->SetTitle(processName.c_str());
        hTemp2->Draw();
        //line->Draw();
        cs->SaveAs((path+string(h[j]->GetName())+"_SBInt_"+channel+".png").c_str());
        cs->Modified();
        cs->Clear();
        hTemp3->GetYaxis()->SetTitle("Full Significance Z");
        hTemp3->SetTitle(processName.c_str());
        hTemp3->Draw();
        line->Draw();
        cs->SaveAs((path+string(h[j]->GetName())+"_FullZ_"+channel+".png").c_str());
        cs->Modified();
        cs->Clear();
     }
     /*if(channel=="emu")*/ HTCut.push_back(2100.0);HTCut.push_back(3200.0);
     if(HTCut.size()!=0){
       sort( HTCut.begin(), HTCut.end() );
       HTCut.erase( unique( HTCut.begin(), HTCut.end() ), HTCut.end() );
       for(UInt_t icut=0;icut<HTCut.size();icut++){
          std::cout << "This is the HT cut: " << HTCut.at(icut) << std::endl;
       }
       std::cout << "HTCut Value: " << HTCut.at(HTCut.size()-1) << std::endl;
       HTCutValue = HTCut.at(HTCut.size()-1);
       UInt_t bin=h[42]->FindBin(HTCutValue);
       std::cout << "Bin " << h[42]->FindBin(HTCutValue) << " Integral: " << h[42]->Integral(bin,h[42]->GetXaxis()->GetNbins()+1) << std::endl;
     }
   }

}

void stackThemOut(std::string channel)
{

  std::cout << "Channel: " << channel << std::endl;

  std::vector<std::string> process;
  if(channel!="mumu")process.push_back("ttbarlep");
  if(channel=="ee") process.push_back("zjets_mllcut_all");
  //process.push_back("ttw");
  process.push_back("ttbarw-v2");
  process.push_back("ttz");
  if(channel!="mumu")process.push_back("twlep");
  if(channel!="mumu")process.push_back("ww");
  process.push_back("wwjj");
  //process.push_back("zw");
  process.push_back("wz-3j-v2");
  process.push_back("zz");
  if(channel=="ee"){
  process.push_back("signal1.5PSI2.5LQee");
  process.push_back("signal1.6PSI2.5LQee");
  process.push_back("signal1.7PSI2.5LQee");
  process.push_back("signal1.8PSI2.5LQee");
  process.push_back("signal1.9PSI2.5LQee");
  process.push_back("signal2.0PSI2.5LQee");
  process.push_back("signal2.1PSI2.5LQee");
  process.push_back("signal2.2PSI2.5LQee");
  process.push_back("signal2.3PSI2.5LQee");
  process.push_back("signal2.4PSI2.5LQee");
  process.push_back("signal2.5PSI3.0LQee");
  process.push_back("signal2.6PSI3.0LQee");
  process.push_back("signal2.7PSI3.0LQee");
  process.push_back("signal2.8PSI3.0LQee");}
  process.push_back("signal1.5PSI2.5LQemu");
  process.push_back("signal1.6PSI2.5LQemu");
  process.push_back("signal1.7PSI2.5LQemu");
  process.push_back("signal1.8PSI2.5LQemu");
  process.push_back("signal1.9PSI2.5LQemu");
  process.push_back("signal2.0PSI2.5LQemu");
  process.push_back("signal2.1PSI2.5LQemu");
  process.push_back("signal2.2PSI2.5LQemu");
  process.push_back("signal2.3PSI2.5LQemu");
  process.push_back("signal2.4PSI2.5LQemu");
  process.push_back("signal2.5PSI3.0LQemu");
  process.push_back("signal2.6PSI3.0LQemu");
  process.push_back("signal2.7PSI3.0LQemu");
  process.push_back("signal2.8PSI3.0LQemu");
  if(channel=="mumu"){
  process.push_back("signal1.5PSI2.5LQmumu");
  process.push_back("signal1.6PSI2.5LQmumu");
  process.push_back("signal1.7PSI2.5LQmumu");
  process.push_back("signal1.8PSI2.5LQmumu");
  process.push_back("signal1.9PSI2.5LQmumu");
  process.push_back("signal2.0PSI2.5LQmumu");
  process.push_back("signal2.1PSI2.5LQmumu");
  process.push_back("signal2.2PSI2.5LQmumu");
  process.push_back("signal2.3PSI2.5LQmumu");
  process.push_back("signal2.4PSI2.5LQmumu");
  process.push_back("signal2.4PSI3.0LQmumu");
  process.push_back("signal2.5PSI3.0LQmumu");
  process.push_back("signal2.6PSI3.0LQmumu");
  process.push_back("signal2.7PSI3.0LQmumu");
  process.push_back("signal2.8PSI3.0LQmumu");}

  std::vector<std::string> dataset;
  for(UInt_t iproc=0;iproc<process.size();iproc++){
   dataset.push_back(path+"hists_"+channel+"_"+process.at(iproc)+".root");
  }

  std::map<std::string,std::string> label_map;
  label_map["ttbarlep"]="t#bar{t}";
  label_map["zjets_mllcut_all"]="Z+jets";
  //label_map["ttw"]="t#bar{t}W";
  label_map["ttbarw-v2"]="t#bar{t}W";
  label_map["ttz"]="t#bar{t}Z";
  label_map["twlep"]="tW";
  label_map["ww"]="WW";
  label_map["wwjj"]="WWjj";
  //label_map["zw"]="ZW";
  label_map["wz-3j-v2"]="ZW";
  label_map["zz"]="ZZ";
  label_map["signal1.5PSI2.5LQee"]="M_{#Psi}=1.5 TeV"/*,  M_{LQ}=2.5 ee*/;
  label_map["signal1.6PSI2.5LQee"]="M_{#Psi}=1.6 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal1.7PSI2.5LQee"]="M_{#Psi}=1.7 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal1.8PSI2.5LQee"]="M_{#Psi}=1.8 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal1.9PSI2.5LQee"]="M_{#Psi}=1.9 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal2.0PSI2.5LQee"]="M_{#Psi}=2.0 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal2.1PSI2.5LQee"]="M_{#Psi}=2.1 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal2.2PSI2.5LQee"]="M_{#Psi}=2.2 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal2.3PSI2.5LQee"]="M_{#Psi}=2.3 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal2.4PSI2.5LQee"]="M_{#Psi}=2.4 TeV"/*, M_{LQ}=2.5 ee*/;
  label_map["signal2.5PSI3.0LQee"]="M_{#Psi}=2.5 TeV"/*, M_{LQ}=3.0 ee*/;
  label_map["signal2.6PSI3.0LQee"]="M_{#Psi}=2.6 TeV"/*, M_{LQ}=3.0 ee*/;
  label_map["signal2.7PSI3.0LQee"]="M_{#Psi}=2.7 TeV"/*, M_{LQ}=3.0 ee*/;
  label_map["signal2.8PSI3.0LQee"]="M_{#Psi}=2.8 TeV"/*, M_{LQ}=3.0 ee*/;

  label_map["signal1.5PSI2.5LQemu"]="M_{#Psi}^{e#mu}=1.5 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal1.6PSI2.5LQemu"]="M_{#Psi}^{e#mu}=1.6 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal1.7PSI2.5LQemu"]="M_{#Psi}^{e#mu}=1.7 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal1.8PSI2.5LQemu"]="M_{#Psi}^{e#mu}=1.8 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal1.9PSI2.5LQemu"]="M_{#Psi}^{e#mu}=1.9 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal2.0PSI2.5LQemu"]="M_{#Psi}^{e#mu}=2.0 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal2.1PSI2.5LQemu"]="M_{#Psi}^{e#mu}=2.1 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal2.2PSI2.5LQemu"]="M_{#Psi}^{e#mu}=2.2 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal2.3PSI2.5LQemu"]="M_{#Psi}^{e#mu}=2.3 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal2.4PSI2.5LQemu"]="M_{#Psi}^{e#mu}=2.4 TeV"/*, M_{LQ}=2.5 e#mu*/;
  label_map["signal2.5PSI3.0LQemu"]="M_{#Psi}^{e#mu}=2.5 TeV"/*, M_{LQ}=3.0 e#mu*/;
  label_map["signal2.6PSI3.0LQemu"]="M_{#Psi}^{e#mu}=2.6 TeV"/*, M_{LQ}=3.0 e#mu*/;
  label_map["signal2.7PSI3.0LQemu"]="M_{#Psi}^{e#mu}=2.7 TeV"/*, M_{LQ}=3.0 e#mu*/;
  label_map["signal2.8PSI3.0LQemu"]="M_{#Psi}^{e#mu}=2.8 TeV"/*, M_{LQ}=3.0 e#mu*/;

  label_map["signal1.5PSI2.5LQmumu"]="M_{#Psi}=1.5 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal1.6PSI2.5LQmumu"]="M_{#Psi}=1.6 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal1.7PSI2.5LQmumu"]="M_{#Psi}=1.7 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal1.8PSI2.5LQmumu"]="M_{#Psi}=1.8 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal1.9PSI2.5LQmumu"]="M_{#Psi}=1.9 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal2.0PSI2.5LQmumu"]="M_{#Psi}=2.0 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal2.1PSI2.5LQmumu"]="M_{#Psi}=2.1 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal2.2PSI2.5LQmumu"]="M_{#Psi}=2.2 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal2.3PSI2.5LQmumu"]="M_{#Psi}=2.3 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal2.4PSI2.5LQmumu"]="M_{#Psi}=2.4 TeV"/*, M_{LQ}=2.5 #mu#mu*/;
  label_map["signal2.4PSI3.0LQmumu"]="M_{#Psi}=2.4 TeV"/*, M_{LQ}=3.0 #mu#mu*/;
  label_map["signal2.5PSI3.0LQmumu"]="M_{#Psi}=2.5 TeV"/*, M_{LQ}=3.0 #mu#mu*/;
  label_map["signal2.6PSI3.0LQmumu"]="M_{#Psi}=2.6 TeV"/*, M_{LQ}=3.0 #mu#mu*/;
  label_map["signal2.7PSI3.0LQmumu"]="M_{#Psi}=2.7 TeV"/*, M_{LQ}=3.0 #mu#mu*/;
  label_map["signal2.8PSI3.0LQmumu"]="M_{#Psi}=2.8 TeV"/*, M_{LQ}=3.0 #mu#mu*/;

  doStackPlot(dataset,"hist_SR_LQMassE",process,1000,"m_{ejj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_LQMassMU",process,1000,"m_{#mu jj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_emujjjjMass",process,1000, "m_{e#mujjjj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_eejjjjMass",process,1000, "m_{eejjjj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_mumujjjjMass",process,1000, "m_{#mu#mujjjj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_ElePT",process,1000, "p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_MuPT",process,1000, "p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_JetPT",process,1000, "p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Heejjjj",process,80, "H_{T}=p_{T}^{e1}+p_{T}^{e2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Hemujjjj",process,80, "H_{T}=p_{T}^{e}+p_{T}^{#mu}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Hmumujjjj",process,80, "H_{T}=p_{T}^{#mu1}+p_{T}^{#mu2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_MET",process,1000, "E_{T}^{miss} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_ElePT_chflip",process,1000,"p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_JetPT_chflip",process,1000,"p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_MuPT_chflip",process,1000,"p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_nJets",process,1000,"N Jets",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Mee",process,1000,"M_{ee}",label_map,channel,1);

  doStackPlot(dataset,"hist_CR_ElePT",process,1000, "p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_MuPT",process,1000, "p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_MET",process,1000, "E_{T}^{miss} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_ElePT_chflip",process,1000,"p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_MuPT_chflip",process,1000,"p_{T}^{#mu} [GeV]",label_map,channel,1);

  doStackPlot(dataset,"hist_VR_ElePT",process,1000, "p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_MuPT",process,1000, "p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_JetPT",process,1000, "p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_MET",process,1000, "E_{T}^{miss} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_ElePT_chflip",process,1000,"p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_JetPT_chflip",process,1000,"p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_MuPT_chflip",process,1000,"p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_nJets",process,1.2,"N Jets",label_map,channel,0);

  std::map<std::string,std::string> label_map2;
  label_map2["ttbarlep"]="$t\\bar{t}$";
  label_map2["zjets_mllcut_all"]="Z+jets";
  //label_map2["ttw"]="$t\\bar{t}W$";
  label_map2["ttbarw-v2"]="$t\\bar{t}W$";
  label_map2["ttz"]="$t\\bar{t}Z$";
  label_map2["twlep"]="tW";
  label_map2["ww"]="WW";
  label_map2["wwjj"]="WWjj";
  //label_map2["zw"]="ZW";
  label_map2["wz-3j-v2"]="ZW";
  label_map2["zz"]="ZZ";
  label_map2["signal1.5PSI2.5LQee"]="1.5 ($ee$)";
  label_map2["signal1.6PSI2.5LQee"]="1.6 ($ee$)";
  label_map2["signal1.7PSI2.5LQee"]="1.7 ($ee$)";
  label_map2["signal1.8PSI2.5LQee"]="1.8 ($ee$)";
  label_map2["signal1.9PSI2.5LQee"]="1.9 ($ee$)";
  label_map2["signal2.0PSI2.5LQee"]="2.0 ($ee$)";
  label_map2["signal2.1PSI2.5LQee"]="2.1 ($ee$)";
  label_map2["signal2.2PSI2.5LQee"]="2.2 ($ee$)";
  label_map2["signal2.3PSI2.5LQee"]="2.3 ($ee$)";
  label_map2["signal2.4PSI2.5LQee"]="2.4 ($ee$)";
  label_map2["signal2.5PSI3.0LQee"]="2.5 ($ee$)";
  label_map2["signal2.6PSI3.0LQee"]="2.6 ($ee$)";
  label_map2["signal2.7PSI3.0LQee"]="2.7 ($ee$)";
  label_map2["signal2.8PSI3.0LQee"]="2.8 ($ee$)";

  label_map2["signal1.5PSI2.5LQemu"]="1.5 ($e\\mu$)";
  label_map2["signal1.6PSI2.5LQemu"]="1.6 ($e\\mu$)";
  label_map2["signal1.7PSI2.5LQemu"]="1.7 ($e\\mu$)";
  label_map2["signal1.8PSI2.5LQemu"]="1.8 ($e\\mu$)";
  label_map2["signal1.9PSI2.5LQemu"]="1.9 ($e\\mu$)";
  label_map2["signal2.0PSI2.5LQemu"]="2.0 ($e\\mu$)";
  label_map2["signal2.1PSI2.5LQemu"]="2.1 ($e\\mu$)";
  label_map2["signal2.2PSI2.5LQemu"]="2.2 ($e\\mu$)";
  label_map2["signal2.3PSI2.5LQemu"]="2.3 ($e\\mu$)";
  label_map2["signal2.4PSI2.5LQemu"]="2.4 ($e\\mu$)";
  label_map2["signal2.5PSI3.0LQemu"]="2.5 ($e\\mu$)";
  label_map2["signal2.6PSI3.0LQemu"]="2.6 ($e\\mu$)";
  label_map2["signal2.7PSI3.0LQemu"]="2.7 ($e\\mu$)";
  label_map2["signal2.8PSI3.0LQemu"]="2.8 ($e\\mu$)";

  label_map2["signal1.5PSI2.5LQmumu"]="1.5 ($\\mu\\mu$)";
  label_map2["signal1.6PSI2.5LQmumu"]="1.6 ($\\mu\\mu$)";
  label_map2["signal1.7PSI2.5LQmumu"]="1.7 ($\\mu\\mu$)";
  label_map2["signal1.8PSI2.5LQmumu"]="1.8 ($\\mu\\mu$)";
  label_map2["signal1.9PSI2.5LQmumu"]="1.9 ($\\mu\\mu$)";
  label_map2["signal2.0PSI2.5LQmumu"]="2.0 ($\\mu\\mu$)";
  label_map2["signal2.1PSI2.5LQmumu"]="2.1 ($\\mu\\mu$)";
  label_map2["signal2.2PSI2.5LQmumu"]="2.2 ($\\mu\\mu$)";
  label_map2["signal2.3PSI2.5LQmumu"]="2.3 ($\\mu\\mu$)";
  label_map2["signal2.4PSI2.5LQmumu"]="2.4 ($\\mu\\mu$)";
  label_map2["signal2.4PSI3.0LQmumu"]="2.4 ($\\mu\\mu$)";
  label_map2["signal2.5PSI3.0LQmumu"]="2.5 ($\\mu\\mu$)";
  label_map2["signal2.6PSI3.0LQmumu"]="2.6 ($\\mu\\mu$)";
  label_map2["signal2.7PSI3.0LQmumu"]="2.7 ($\\mu\\mu$)";
  label_map2["signal2.8PSI3.0LQmumu"]="2.8 ($\\mu\\mu$)";

  std::cout << "HT Cut Vector Size: " << HTCut.size() << std::endl;

  std::map<int,std::vector<float> > cutmap;

  std::vector<float> weightsHT;

  for(UInt_t icut=0;icut<HTCut.size();icut++){
     weightsHT.clear();
     HTCutValue=HTCut.at(icut);
     std::string histo_name="hist_SR_H"+channel+"jjjj";
     weightsHT=fillHTCutFlowMap(dataset,histo_name,process,HTCutValue);
     //else if(channel=="emu") weightsHT=fillHTCutFlowMap(dataset,"hist_SR_Hemujjjj",process,HTCutValue);
     //else if(channel=="mumu") weightsHT=fillHTCutFlowMap(dataset,"hist_SR_Hmumujjjj",process,HTCutValue);
     cutmap[icut]=weightsHT;
  }

  writeCutFlow(dataset,"hist_CutFlow_Weighted",process,label_map2,channel,cutmap);
  writeCutFlow(dataset,"hist_CutFlow",process,label_map2,channel,cutmap);

  const char* command;
  command=("tar -zcvf "+path+channel+"Plots.tar.gz "+path+"/*"+channel+".png ").c_str();
  gSystem->Exec(command);
  command=("rm "+path+"/*.png").c_str();
  //gSystem->Exec(command);

}

void StackPlotsNew()
{
  gROOT->SetBatch(1);
  //gSystem->Load("libDelphes");

  f = TFile::Open((path+"GetScaleFactor.root").c_str(),"RECREATE");

  stackThemOut("ee");
  HTCut.clear();
  stackThemOut("emu");
  HTCut.clear();
  stackThemOut("mumu");
  HTCut.clear();

  f->Write();
  f->Close();
}
