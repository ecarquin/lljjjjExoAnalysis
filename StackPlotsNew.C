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

#include <math.h>
#include <cmath>
#include <cstdlib>
#include <sstream>


//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
//#include "classes/DelphesClasses.h"
//#include "external/ExRootAnalysis/ExRootTreeReader.h"
//#include "external/ExRootAnalysis/ExRootResult.h"
//#endif

#include <iostream>

std::string path="09_03_17/";

Color_t ci[] = {kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack,kWhite,kBlack};
Color_t cf[] = {kMagenta+3,kGreen-3,kGreen+5,kCyan,kMagenta-3,kGreen,kBlue,kRed,kYellow,kBlack,kBlue-3, kBlue+3, kMagenta+3,kGreen-3,kGreen+5,kCyan,kMagenta-3,kGreen,kBlue,kRed,kYellow,kBlack,kBlue-3, kBlue+3, kMagenta+3};

/*
root -l PlotMET_BKG.C\(\"delphes_output.root\"\)
*/

TFile *f;

void doStackPlot(std::vector<std::string> dataset, std::string hist_name, std::vector<std::string> process,float scaleF,std::string xLabel, std::map<std::string,std::string> label_map, std::string channel, bool isLog)
{
  std::cout << "Channel: " << channel << std::endl;
  THStack *hs = new THStack("hs","");

  TH1F* h[50];
  std::vector<std::string> newname;

  std::stringstream temp_name;
  float maxi = 0;

  gStyle->SetPadRightMargin(0.1);
  TCanvas *c = new TCanvas("c", "", 800, 600);
  gROOT->SetStyle("Plain");
  gStyle->SetTitleBorderSize(0);
  TLegend *leg = new TLegend(0.11, 0.70, 0.89, 0.89);
  leg->SetNColumns(3);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);

  int index=0;
  int index1=0;

  std::vector<std::string> labels;
  labels.clear();

  for(Int_t i = dataset.size()-1; i >= 0; --i){
    //for(Int_t i = 0; i < dataset.size(); i++){
    std::cout << "i: " << i << std::endl;
    // Only stack background histograms
    if(process.at(i).find("signal") != std::string::npos ) continue;
    index1++;
    temp_name << hist_name << "_" << process.at(i);
    std::cout << "temp_name: " << temp_name.str() << std::endl;
    TFile *f = TFile::Open(dataset.at(i).c_str(),"READ");
    TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
    if(htemp->Integral()==0) {temp_name.str("");std::cout << "Histogram have zero integral" << std::endl; continue;}
    labels.push_back(label_map[process.at(i)]);
    h[index]=(TH1F*)htemp->Clone(temp_name.str().c_str());
    h[index]->SetLineWidth(2);
    h[index]->SetLineColor(ci[i]);
    h[index]->SetFillColor(cf[i]/*rainbow[ki]*/);
    std::cout << "Now running on: " << i << " counter " << index << " " << process.at(i) << " Integral: " << htemp->Integral() << std::endl;
    h[index]->SetFillStyle(1001);
    h[index]->SetLineStyle(1/*ls[ki]*/);
    h[index]->SetStats(0);
    //h[i]->SetYTitle("# of events/100 fb^{-1}");
    //h[i]->SetXTitle(xLabel.c_str());
    h[index]->SetTitle("");
    //h[i]->Scale(1./6);
    maxi += h[index]->GetMaximum();
    std::cout << "This is the maximum: " << maxi << " i " << i << std::endl;
    hs->Add(h[index]);
    std::cout << "After adding the histogram" << std::endl;
    temp_name.str("");
    index++;
  }

  if(maxi==0) return;

  std::cout << "After the first loop " << index1 << std::endl;
  int index0=index;

  for(UInt_t k=index1;k<dataset.size();k++){
   if(process.at(k).find("signal") == std::string::npos ) continue;
   temp_name << hist_name << "_" << process.at(k);
   TFile *f = TFile::Open(dataset.at(k).c_str(),"READ");
   std::cout << "dataset: " << dataset.at(k) << std::endl;
   TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
   if(htemp->Integral()==0) {temp_name.str("");continue;}
   labels.push_back(label_map[process.at(k)]);
   std::cout << "Now running on: " << k << " counter " << index0 << " " << process.at(k) << " Integral: " << htemp->Integral() << " Name: " << temp_name.str() << std::endl;
   h[index0]=(TH1F*)htemp->Clone(temp_name.str().c_str());
   h[index0]->SetLineStyle(1);
   h[index0]->SetLineColor(cf[k]);
   h[index0]->SetLineWidth(4);
   h[index0]->SetStats(0);
   index0++;
   temp_name.str("");
  }

  std::cout << "After the second loop on signal samples: " << index << " index0: " << index0 << " labels size " << labels.size() << std::endl;
  for(UInt_t l=0;l<labels.size();l++) std::cout << "label: " << labels.at(l) << std::endl;

   for (UInt_t ki = 0; ki < index; ++ki) {
    std::cout << "Adding legends: " << ki << " " << labels.at(ki) << std::endl;
    std::cout << "Check this ! " << h[ki]->GetName() << std::endl;
    if(h[ki]->Integral()==0) continue;
    leg->AddEntry(h[ki], labels.at(ki).c_str(), "F");
   }

   for (UInt_t ki = index; ki < index0; ++ki) {
    std::cout << "Adding legends: " << ki << " " << process.at(ki) << std::endl;
    std::cout << "Check this ! " << h[ki]->GetName() << std::endl;
    if(h[ki]->Integral()==0) continue;
    leg->AddEntry(h[ki], labels.at(ki).c_str(), "F");
   }

   TCanvas *cs = new TCanvas("cs","",600,800);
   if(isLog) cs->SetLogy();
   hs->SetMaximum(scaleF*maxi);
   hs->SetMinimum(1e-3);
   
   hs->SetMinimum(1.0);
   hs->Draw("HIST");
   leg->Draw();

   for(UInt_t j=index;j<index0;j++) {
      /*if(j==index) h[j]->Draw("HIST");
      else*/ 
      h[j]->Draw("HIST:same");
      std::cout << "hist_name: " << h[j]->GetName() << " Integral: " << h[j]->Integral() << std::endl;
   }

   TPaveText *pt = new TPaveText(0.6,0.6,0.7,0.7,"NDC");
   pt->AddText(("SS "+channel).c_str());
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->Draw();

   hs->GetXaxis()->SetTitle(xLabel.c_str());
   hs->GetYaxis()->SetTitleOffset(1.2);
   hs->GetYaxis()->SetTitle("# of events/100 fb^{-1}");
   std::stringstream histname;
   histname << hist_name << "_" << channel << ".png";
   //histname << hist_name << "_nosamescut_bkg.png";
   

   cs->SaveAs((path+histname.str()).c_str());

   if(hist_name.find("VR_nJets") != std::string::npos){
     f->cd();
     if(channel == "ee"){
       std::cout << "Here ee" << hist_name << std::endl;
       h[40] =(TH1F*)(hs->GetStack()->Last());
       std::cout << h[40]->GetName() << std::endl;
       h[40]->Write((hist_name+"_"+channel).c_str());
     }else if(channel == "mumu"){
       std::cout << "Here mm" << std::endl;
       h[41] =(TH1F*)(hs->GetStack()->Last());
       std::cout << h[41]->GetName() << std::endl;
       h[41]->Write((hist_name+"_"+channel).c_str());
     }
   }

}

void stackThemOut(std::string channel)
{

  std::vector<std::string> process;
  process.push_back("ttbarlep");
  //process.push_back("zjjhpT");
  process.push_back("zjjBig");
  process.push_back("ttw");
  process.push_back("ttz");
  process.push_back("twlep");
  process.push_back("ww");
  process.push_back("wwjj");
  process.push_back("zw");
  process.push_back("zz");
  process.push_back("signal1.5PSI2.5LQee");
  process.push_back("signal1.5PSI2.5LQemu");
  process.push_back("signal1.6PSI2.5LQemu");
  process.push_back("signal1.7PSI2.5LQemu");
  process.push_back("signal1.8PSI2.5LQemu");
  process.push_back("signal1.9PSI2.5LQemu");
  process.push_back("signal2.0PSI1.5LQemu");
  process.push_back("signal2.0PSI2.5LQee");
  process.push_back("signal2.1PSI2.5LQemu");
  process.push_back("signal2.2PSI2.5LQemu");
  process.push_back("signal2.3PSI2.5LQemu");
  process.push_back("signal2.4PSI2.5LQemu");
  process.push_back("signal2.5PSI2.6LQemu");
  process.push_back("signal1.5PSI2.5LQmumu");
  process.push_back("signal2.0PSI2.5LQmumu");

  std::vector<std::string> dataset;
  for(UInt_t iproc=0;iproc<process.size();iproc++)
   dataset.push_back(path+"hists_"+channel+"_"+process.at(iproc)+".root");

  std::map<std::string,std::string> label_map;
  label_map["ttbarlep"]="t#bar{t}";
  //label_map["zjjhpT"]="Z+jets high pT";
  label_map["zjjBig"]="Z+jets";
  label_map["ttw"]="t#bar{t}W";
  label_map["ttz"]="t#bar{t}Z";
  label_map["twlep"]="tW";
  label_map["ww"]="WW";
  label_map["wwjj"]="WWjj";
  label_map["zw"]="ZW";
  label_map["zz"]="ZZ";
  label_map["signal1.5PSI2.5LQee"]="M_{#Psi}=1.5 TeV, M_{LQ}=2.5 ee";
  label_map["signal1.5PSI2.5LQemu"]="M_{#Psi}=1.5 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal1.6PSI2.5LQemu"]="M_{#Psi}=1.6 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal1.7PSI2.5LQemu"]="M_{#Psi}=1.7 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal1.8PSI2.5LQemu"]="M_{#Psi}=1.8 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal1.9PSI2.5LQemu"]="M_{#Psi}=1.9 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal2.0PSI1.5LQemu"]="M_{#Psi}=2.0 TeV, M_{LQ}=1.5 e#mu";
  label_map["signal2.0PSI2.5LQee"]="M_{#Psi}=2.0 TeV, M_{LQ}=2.5 ee";
  label_map["signal2.1PSI2.5LQemu"]="M_{#Psi}=2.1 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal2.2PSI2.5LQemu"]="M_{#Psi}=2.2 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal2.3PSI2.5LQemu"]="M_{#Psi}=2.3 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal2.4PSI2.5LQemu"]="M_{#Psi}=2.4 TeV, M_{LQ}=2.5 e#mu";
  label_map["signal2.5PSI2.6LQemu"]="M_{#Psi}=2.5 TeV, M_{LQ}=2.6 e#mu";
  label_map["signal1.5PSI2.5LQmumu"]="M_{#Psi}=1.5 TeV, M_{LQ}=2.5 #mu#mu";
  label_map["signal2.0PSI2.5LQmumu"]="M_{#Psi}=2.0 TeV, M_{LQ}=2.5 #mu#mu";

  doStackPlot(dataset,"hist_SR_LQMassE",process,10,"m_{ejj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_LQMassMU",process,10,"m_{#mu jj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_emujjjjMass",process,10, "m_{e#mujjjj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_eejjjjMass",process,10, "m_{eejjjj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_mumujjjjMass",process,10, "m_{#mu#mujjjj} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_ElePT",process,10, "p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_MuPT",process,10, "p_{T}^{#mu] [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_JetPT",process,10, "p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Heejjjj",process,10, "H=p_{T}^{e1}+p_{T}^{e2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Hemujjjj",process,10, "H=p_{T}^{e}+p_{T}^{#mu}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Hmumujjjj",process,10, "H=p_{T}^{#mu1}+p_{T}^{#mu2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_MET",process,10, "E_{T}^{miss} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_ElePT_chflip",process,10,"p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_JetPT_chflip",process,10,"p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_MuPT_chflip",process,10,"p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_nJets",process,10,"N Jets",label_map,channel,1);
  doStackPlot(dataset,"hist_SR_Mee",process,10,"M_{ee}",label_map,channel,1);

  doStackPlot(dataset,"hist_CR_ElePT",process,10, "p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_MuPT",process,10, "p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_MET",process,10, "E_{T}^{miss} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_ElePT_chflip",process,10,"p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_CR_MuPT_chflip",process,10,"p_{T}^{#mu} [GeV]",label_map,channel,1);

  doStackPlot(dataset,"hist_VR_ElePT",process,10, "p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_MuPT",process,10, "p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_JetPT",process,10, "p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_MET",process,10, "E_{T}^{miss} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_ElePT_chflip",process,10,"p_{T}^{e} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_JetPT_chflip",process,10,"p_{T}^{jet} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_MuPT_chflip",process,10,"p_{T}^{#mu} [GeV]",label_map,channel,1);
  doStackPlot(dataset,"hist_VR_nJets",process,1.2,"N Jets",label_map,channel,0);

  const char* command;
  command=("tar -zcvf "+path+channel+"Plots.tar.gz "+path+"/*"+channel+".png ").c_str();
  gSystem->Exec(command);
  command=("rm "+path+"/*.png").c_str();
  gSystem->Exec(command);

}

void StackPlotsNew()
{
  gROOT->SetBatch(1);
  gSystem->Load("libDelphes");

  f = TFile::Open((path+"GetScaleFactor.root").c_str(),"RECREATE");

  stackThemOut("ee");
  stackThemOut("emu");
  stackThemOut("mumu");

  f->Write();
  f->Close();
}
