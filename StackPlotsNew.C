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

std::string path="22_10_16/";

Color_t ci[] = {kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack, kBlack,kWhite,kBlack};
Color_t cf[] = {kMagenta+3,kGreen-3,kGreen+5,kCyan,kMagenta-3,kGreen,kBlue,kRed,kYellow,kBlack,kBlue-3, kBlue+3, kMagenta+3,kGreen-3,kGreen+5,kCyan,kMagenta-3,kGreen,kBlue,kRed,kYellow,kBlack,kBlue-3, kBlue+3};

/*
root -l PlotMET_BKG.C\(\"delphes_output.root\"\)
*/

TFile *f;

void doStackPlot(std::vector<std::string> dataset, std::string hist_name, std::vector<std::string> process,float scaleF,std::string xLabel, std::vector<std::string> labels, std::string channel, bool isLog)
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

  for(Int_t i = dataset.size()-1; i >= 0; --i){
    //std::cout << "i: " << i << std::endl;
    // Only stack background histograms
    if(process.at(i).find("signal") != std::string::npos ) continue;
    temp_name << hist_name << "_" << process.at(i);
    TFile *f = TFile::Open(dataset.at(i).c_str(),"READ");
    TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
    h[index]=(TH1F*)htemp->Clone(temp_name.str().c_str());
    h[index]->SetLineWidth(2);
    h[index]->SetLineColor(ci[i]);
    h[index]->SetFillColor(cf[i]/*rainbow[ki]*/);
    std::cout << "Now running on: " << i << " counter " << i << " " << process.at(i) << std::endl;
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

  std::cout << "After the first loop " << index << std::endl;
  int index0=index;

  for(UInt_t k=0;k<dataset.size();k++){
   if(process.at(k).find("signal") == std::string::npos ) continue;
   temp_name << hist_name << "_" << process.at(k);
   TFile *f = TFile::Open(dataset.at(index).c_str(),"READ");
   TH1F *htemp = (TH1F*)f->Get(hist_name.c_str());
   h[index0]=(TH1F*)htemp->Clone(temp_name.str().c_str());
   h[index0]->SetLineStyle(1);
   h[index0]->SetLineColor(cf[k]);
   h[index0]->SetLineWidth(4);
   h[index0]->SetStats(0);
   index0++;
  }

  std::cout << "After the second loop on signal samples: " << index << " index0: " << index0 << std::endl;

   for (UInt_t ki = 0; ki < index0-1; ++ki) {
    std::cout << "Adding legends: " << ki << std::endl;
    leg->AddEntry(h[ki], labels.at(ki).c_str(), "F");
   }

   TCanvas *cs = new TCanvas("cs","",600,800);
   //hs->GetYaxis()->SetRangeUser(0,200);
   if(isLog) cs->SetLogy();
   hs->SetMaximum(scaleF*maxi);
   
   hs->SetMinimum(1.0);
   hs->Draw("HIST");
   leg->Draw();

   for(UInt_t j=index;j<index0-1;j++) h[j]->Draw("HIST:same");

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
     h[40] =(TH1F*)(hs->GetStack()->Last());
     h[40]->Write((hist_name+"_"+channel).c_str());
   }

}

void stackThemOut(std::string channel)
{

  f = TFile::Open("GetScaleFactor.root","RECREATE");

  std::vector<std::string> process;
  process.push_back("ttbarlep");
  process.push_back("zjj");
  process.push_back("zjjhpT");
  process.push_back("ttw");
  process.push_back("ttz");
  process.push_back("twlep");
  process.push_back("ww");
  /*if(channel != "emu")*/ process.push_back("wwjj");
  process.push_back("zw");
  process.push_back("zz");
  if(channel != "emu" && channel != "mumu") process.push_back("signal1.5PSI2.5LQee");
  if(channel != "mumu") process.push_back("signal1.5PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal1.6PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal1.7PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal1.8PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal1.9PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal2.0PSI1.5LQemu");
  if(channel != "emu" && channel != "mumu") process.push_back("signal2.0PSI2.5LQee");
  if(channel != "mumu") process.push_back("signal2.1PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal2.2PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal2.3PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal2.4PSI2.5LQemu");
  if(channel != "mumu") process.push_back("signal2.5PSI2.6LQemu");
  /*if(channel != "mumu")*/ process.push_back("signal1.5PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal1.6PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal1.7PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal1.8PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal1.9PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal2.0PSI1.5LQmumu");
  /*if(channel != "mumu")*/ process.push_back("signal2.0PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal2.1PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal2.2PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal2.3PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal2.4PSI2.5LQmumu");
  //if(channel != "ee" && channel != "emu") process.push_back("signal2.5PSI2.6LQmumu");

  std::vector<std::string> dataset;
  for(UInt_t iproc=0;iproc<process.size();iproc++)
   dataset.push_back(path+"hists_"+channel+"_"+process.at(iproc)+".root");

  std::vector<std::string> labels;
  labels.push_back("t#bar{t}");
  //labels.push_back("ttbarhad");
  labels.push_back("Z+jets");
  labels.push_back("Z+jets high pT");
  labels.push_back("t#bar{t}W");
  labels.push_back("t#bar{t}Z");
  labels.push_back("tW");
  //labels.push_back("twhad");
  labels.push_back("WW");
  /*if(channel != "emu")*/ labels.push_back("WWjj");
  labels.push_back("ZW");
  labels.push_back("ZZ");
  if(channel != "emu" && channel != "mumu") labels.push_back("M_{#Psi}=1.5 TeV, M_{LQ}=2.5 ee");
  if(channel != "mumu") labels.push_back("M_{#Psi}=1.5 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=1.6 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=1.7 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=1.8 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=1.9 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=2.0 TeV, M_{LQ}=1.5 e#mu");
  if(channel != "emu" && channel != "mumu") labels.push_back("M_{#Psi}=2.0 TeV, M_{LQ}=2.5 ee");
  if(channel != "mumu") labels.push_back("M_{#Psi}=2.1 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=2.2 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=2.3 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=2.4 TeV, M_{LQ}=2.5 e#mu");
  if(channel != "mumu") labels.push_back("M_{#Psi}=2.5 TeV, M_{LQ}=2.6 e#mu");

  if(channel != "ee" && channel != "emu") labels.push_back("M_{#Psi}=1.5 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*///if(channel != "ee" && channel != "emu") labels.push_back("M_{#Psi}=1.6 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*///if(channel != "ee" && channel != "emu") labels.push_back("M_{#Psi}=1.7 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*///if(channel != "ee" && channel != "emu") labels.push_back("M_{#Psi}=1.8 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*///if(channel != "ee" && channel != "emu") labels.push_back("M_{#Psi}=1.9 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*///if(channel != "ee" && channel != "emu") labels.push_back("M_{#Psi}=2.0 TeV, M_{LQ}=1.5 #mu#mu");
  /*if(channel != "mumu")*/ labels.push_back("M_{#Psi}=2.0 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*/ //labels.push_back("M_{#Psi}=2.1 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*/ //labels.push_back("M_{#Psi}=2.2 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*/ //labels.push_back("M_{#Psi}=2.3 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*/ //labels.push_back("M_{#Psi}=2.4 TeV, M_{LQ}=2.5 #mu#mu");
  /*if(channel != "mumu")*/ //labels.push_back("M_{#Psi}=2.5 TeV, M_{LQ}=2.6 #mu#mu");

  /*doStackPlot(process,"hist_LQMass",labels,2,"m_{ejj} [GeV]");
  doStackPlot(process,"hist_eejjjjMass",labels,10, "m_{eejjjj} [GeV]");
  doStackPlot(process,"hist_ElePT_lead",labels,10, "p_{T} [GeV]");
  doStackPlot(process,"hist_ElePT_subl",labels,10, "p_{T} [GeV]");
  doStackPlot(process,"hist_JetPT",labels,10, "p_{T} [GeV]");
  doStackPlot(process,"hist_Heejjjj",labels,10, "H=p_{T}^{e1}+p_{T}^{e2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]");*/

  doStackPlot(dataset,"hist_SR_LQMassE",process,10,"m_{ejj} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_LQMassMU",process,10,"m_{#mu jj} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_emujjjjMass",process,10, "m_{e#mujjjj} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_eejjjjMass",process,10, "m_{eejjjj} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_mumujjjjMass",process,10, "m_{#mu#mujjjj} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_ElePT",process,10, "p_{T}^{e} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_MuPT",process,10, "p_{T}^{#mu] [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_JetPT",process,10, "p_{T}^{jet} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_Heejjjj",process,10, "H=p_{T}^{e1}+p_{T}^{e2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_Hemujjjj",process,10, "H=p_{T}^{e}+p_{T}^{#mu}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_Hmumujjjj",process,10, "H=p_{T}^{#mu1}+p_{T}^{#mu2}+p_{T}^{j1}+p_{T}^{j2}+p_{T}^{j3}+p_{T}^{j4} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_MET",process,10, "E_{T}^{miss} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_ElePT_chflip",process,10,"p_{T}^{e} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_JetPT_chflip",process,10,"p_{T}^{jet} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_MuPT_chflip",process,10,"p_{T}^{#mu} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_SR_nJets",process,10,"N Jets",labels,channel,1);

  doStackPlot(dataset,"hist_CR_ElePT",process,10, "p_{T}^{e} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_CR_MuPT",process,10, "p_{T}^{#mu} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_CR_MET",process,10, "E_{T}^{miss} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_CR_ElePT_chflip",process,10,"p_{T}^{e} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_CR_MuPT_chflip",process,10,"p_{T}^{#mu} [GeV]",labels,channel,1);

  doStackPlot(dataset,"hist_VR_ElePT",process,10, "p_{T}^{e} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_MuPT",process,10, "p_{T}^{#mu} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_JetPT",process,10, "p_{T}^{jet} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_MET",process,10, "E_{T}^{miss} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_ElePT_chflip",process,10,"p_{T}^{e} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_JetPT_chflip",process,10,"p_{T}^{jet} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_MuPT_chflip",process,10,"p_{T}^{#mu} [GeV]",labels,channel,1);
  doStackPlot(dataset,"hist_VR_nJets",process,1.2,"N Jets",labels,channel,0);

  f->Write();
  f->Close();

}

void StackPlotsNew()
{
  gROOT->SetBatch(1);
  gSystem->Load("libDelphes");

  stackThemOut("ee");
  stackThemOut("emu");
  stackThemOut("mumu");

}
