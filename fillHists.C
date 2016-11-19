#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TMath.h"
#include <math.h>
#include <cmath>
#include <cstdlib>
#include "TLorentzVector.h"
#include "TRandom.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <sstream>

std::string outDir = "22_10_16/";


void processChain(std::string inputFile, std::string tag,std::string channel,std::string path, bool hasMerge, int nSplit, int fileN,std::map<std::string,float>& yields, std::map<std::string,float>& yieldsw)
{

  // Variables
  float lum = 100; //in fb^{-1}
  bool apply_reweight=false;
  float cs;

  int Nevents = 0;
  std::string cs_name = path+tag+"/cs.txt";
  std::cout << "cs_name: " << cs_name << std::endl;
  if(!hasMerge) cs = readCS(cs_name,Nevents);
  else cs = readAvMergedCS(cs_name,Nevents);

  TFile *f = TFile::Open("chargeFlipAtlas.root","READ");
  TH1F* etaBinning = (TH1F*)f->Get("etaBinning");
  TH1F* lowPt = (TH1F*)f->Get("lowPt");
  TH1F* mediumPt = (TH1F*)f->Get("mediumPt");
  TH1F* highPt = (TH1F*)f->Get("highPt");

  TRandom *random = new TRandom();
 
  // Create chain of root trees
  //Add the case in which the input is split in several input Files
  TChain chain("Delphes");
  if(fileN==-1)chain.Add(inputFile.c_str());
  else{
   for(unsigned int i=0;i<nSplit;i++){
      std::stringstream run;
      run <<i;
      std::string file=path+tag+"/Events/run_01_"+run.str()+"/tag_1_delphes_events.root";
      std::cout << file << std::endl;
      if(i==fileN)chain.Add(file.c_str());
      run.str("");
   }
  }
  
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  float xBins[14] = {24.0,30.0,40.0,50.0,60.0,70.0,90.0,110.0,130.0,150.0,200.0,250.0,500.0,800.0};

  TFile *SF = TFile::Open("SF.root","READ");
  TH1F* sf_ee=(TH1F*)SF->Get("sf_ee");
  TH1F* sf_mumu=(TH1F*)SF->Get("sf_mumu");

  // Book histograms
  std::map<std::string, TH1F> TH1Hists;

  if(tag.find("signal") != std::string::npos){
    TH1Hists["hist_SR_ElePT"] = TH1F("hist_SR_ElePT", "Electron P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_SR_ElePT_chflip"] = TH1F("hist_SR_ElePT_chflip", "Electron P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_SR_MuPT"] = TH1F("hist_SR_MuPT", "Muon P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_SR_MuPT_chflip"] = TH1F("hist_SR_MuPT_chflip", "Muon P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_SR_JetPT"] = TH1F("hist_SR_JetPT", "jet P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_SR_JetPT_chflip"] = TH1F("hist_SR_JetPT_chflip", "jet P_{T}", 20, 0.0, 800.0);
  }else{
    TH1Hists["hist_SR_ElePT"] = TH1F("hist_SR_ElePT", "Electron P_{T}", 13, xBins);
    TH1Hists["hist_SR_ElePT_chflip"] = TH1F("hist_SR_ElePT_chflip", "Electron P_{T}", 13, xBins);
    TH1Hists["hist_SR_MuPT"] = TH1F("hist_SR_MuPT", "Muon P_{T}", 13, xBins);
    TH1Hists["hist_SR_MuPT_chflip"] = TH1F("hist_SR_MuPT_chflip", "Muon P_{T}", 13, xBins);
    TH1Hists["hist_SR_JetPT"] = TH1F("hist_SR_JetPT", "jet P_{T}", 13, xBins);
    TH1Hists["hist_SR_JetPT_chflip"] = TH1F("hist_SR_JetPT_chflip", "jet P_{T}", 13, xBins);
  }

  TH1Hists["hist_SR_MET"] = TH1F("hist_SR_MET", "MET", 40, 0.0, 400.0);
  TH1Hists["hist_SR_LQMassE"] = TH1F("hist_SR_LQMassE","LQ mass electron",34,100.0,3500.0);
  TH1Hists["hist_SR_LQMassMU"] = TH1F("hist_SR_LQMassMU","LQ mass muon",34,100.0,3500.0);

  TH1Hists["hist_SR_emujjjjMass"] = TH1F("hist_SR_emujjjjMass","emujjjj Mass",38,200.0,4000.0);
  TH1Hists["hist_SR_Hemujjjj"] = TH1F("hist_SR_Hemujjjj","Hemujjjj ",39,100.0,4000.0);

  TH1Hists["hist_SR_mumujjjjMass"] = TH1F("hist_SR_mumujjjjMass","emujjjj Mass",38,200.0,4000.0);
  TH1Hists["hist_SR_Hmumujjjj"] = TH1F("hist_SR_Hmumujjjj","Hemujjjj ",39,100.0,4000.0);

  TH1Hists["hist_SR_eejjjjMass"] = TH1F("hist_SR_eejjjjMass","eejjjj Mass",38,200.0,4000.0);
  TH1Hists["hist_SR_Heejjjj"] = TH1F("hist_SR_Heejjjj","Heejjjj ",39,100.0,4000.0);

  TH1Hists["hist_SR_nJets"] = TH1F("hist_SR_nJets","N-Jets",6,-0.5,5.5);

  // Book histograms
  TH1Hists["hist_CR_ElePT"] = TH1F("hist_CR_ElePT", "Electron P_{T}", 20, 0.0, 200.0);
  TH1Hists["hist_CR_ElePT_chflip"] = TH1F("hist_CR_ElePT_chflip", "Electron P_{T}", 20, 0.0, 200.0);
  TH1Hists["hist_CR_MuPT"] = TH1F("hist_CR_MuPT", "Muon P_{T}", 20, 0.0, 200.0);
  TH1Hists["hist_CR_MuPT_chflip"] = TH1F("hist_CR_MuPT_chflip", "Muon P_{T}", 20, 0.0, 200.0);
  TH1Hists["hist_CR_MET"] = TH1F("hist_CR_MET", "MET", 20, 0.0, 400.0);

  // Book histograms
  if(tag.find("signal") != std::string::npos){
    TH1Hists["hist_VR_ElePT"] = TH1F("hist_VR_ElePT", "Electron P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_VR_ElePT_chflip"] = TH1F("hist_VR_ElePT_chflip", "Electron P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_VR_MuPT"] = TH1F("hist_VR_MuPT", "Muon P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_VR_MuPT_chflip"] = TH1F("hist_VR_MuPT_chflip", "Muon P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_VR_JetPT"] = TH1F("hist_VR_JetPT", "jet P_{T}", 20, 0.0, 800.0);
    TH1Hists["hist_VR_JetPT_chflip"] = TH1F("hist_VR_JetPT_chflip", "jet P_{T}", 20, 0.0, 800.0);
  }else{
    TH1Hists["hist_VR_ElePT"] = TH1F("hist_VR_ElePT", "Electron P_{T}", 13, xBins);
    TH1Hists["hist_VR_ElePT_chflip"] = TH1F("hist_VR_ElePT_chflip", "Electron P_{T}", 13, xBins);
    TH1Hists["hist_VR_MuPT"] = TH1F("hist_VR_MuPT", "Muon P_{T}", 13, xBins);
    TH1Hists["hist_VR_MuPT_chflip"] = TH1F("hist_VR_MuPT_chflip", "Muon P_{T}", 13, xBins);
    TH1Hists["hist_VR_JetPT"] = TH1F("hist_VR_JetPT", "jet P_{T}", 13, xBins);
    TH1Hists["hist_VR_JetPT_chflip"] = TH1F("hist_VR_JetPT_chflip", "jet P_{T}", 13, xBins);
  }

  TH1Hists["hist_VR_MuEta"] = TH1F("hist_VR_MuEta", "Muon Eta",16,-5,5);
  TH1Hists["hist_VR_JetEta"] = TH1F("hist_VR_JetEta", "jet Eta",16,-5,5);
  TH1Hists["hist_VR_EleEta"] = TH1F("hist_VR_EleEta", "Electron Eta",16,-5,5);
  TH1Hists["hist_VR_MET"] = TH1F("hist_VR_MET", "MET", 20, 0.0, 400.0);
  TH1Hists["hist_VR_nJets"] = TH1F("hist_VR_nJets","N-Jets",6,-0.5,5.5);

  TH1Hists["hist_SR_Mee"] = TH1F("hist_SR_Mee","Mee",100,0,200);

  TH1Hists["hist_SR_LQMassE"].SetYTitle("Events");
  TH1Hists["hist_SR_LQMassE"].SetXTitle("m_{ejj} [GeV]");

  TH1Hists["hist_SR_LQMassMU"].SetYTitle("Events");
  TH1Hists["hist_SR_LQMassMU"].SetXTitle("m_{#mu jj} [GeV]");

  TH1Hists["hist_SR_eejjjjMass"].SetYTitle("Events");
  TH1Hists["hist_SR_eejjjjMass"].SetXTitle("m_{eejjjj} [GeV]");

  TH1Hists["hist_SR_Heejjjj"].SetYTitle("Events");
  TH1Hists["hist_SR_Heejjjj"].SetXTitle("H [GeV]");

  TH1Hists["hist_SR_emujjjjMass"].SetYTitle("Events");
  TH1Hists["hist_SR_emujjjjMass"].SetXTitle("m_{eejjjj} [GeV]");

  TH1Hists["hist_SR_Hemujjjj"].SetYTitle("Events");
  TH1Hists["hist_SR_Hemujjjj"].SetXTitle("H [GeV]");

  TH1Hists["hist_SR_mumujjjjMass"].SetYTitle("Events");
  TH1Hists["hist_SR_mumujjjjMass"].SetXTitle("m_{eejjjj} [GeV]");

  TH1Hists["hist_SR_Hmumujjjj"].SetYTitle("Events");
  TH1Hists["hist_SR_Hmumujjjj"].SetXTitle("H [GeV]");

  TH1Hists["hist_SR_Mee"].SetYTitle("Events");
  TH1Hists["hist_SR_Mee"].SetXTitle("M [GeV]");

  std::map<std::string, TH1F>::iterator it0 = TH1Hists.begin();
        for( ; it0 != TH1Hists.end(); ++it0 ) it0->second.Sumw2();

  std::string fname;

  if(fileN!=-1){
     std::stringstream ss;
     ss << fileN;
     std::string N= ss.str();
     fname = "hists_"+channel+"_"+tag+"_"+N+".root";
  }else  
    fname = "hists_"+channel+"_"+tag+".root";

  std::cout << "After that" << std::endl;
  TFile *f0 = TFile::Open((outDir+fname).c_str(),"RECREATE");

  if( (tag.find("signal") != std::string::npos) || (tag.find("wwjj") != std::string::npos))
    numberOfEntries=20000;

  // Variables
  //float lum = 100; //in fb^{-1}
  if(fileN==-1) Nevents  = numberOfEntries;
  float weight = cs * lum * 1e3 / Nevents;
  float rand = 0.0;
  int njets=0;

  std::cout << "Will run over: " << numberOfEntries << " number of entries" << " cross section is: " << cs << " Nevents of total sample is: " << Nevents << std::endl; 

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    if(entry%1000==0){
      std::cout << "Entry: " << entry << std::endl;
      const char* command="date \"+%H:%M:%S   %d/%m/%y\"";
      gSystem->Exec(command);
    }
    //if(entry==100) break;

    //Get electrons
    std::vector<Electron*> electrons;
    for(unsigned int i=0;i<branchElectron->GetEntries();i++){
       Electron *el=(Electron*)branchElectron->At(i);
       if(fabs(el->Eta) < 2.47 && el->PT > 20) electrons.push_back(el);
    }

    //Get muons
    std::vector<Muon*> muons;
    for(unsigned int i=0;i<branchMuon->GetEntries();i++){
       Muon *mu=(Muon*)branchMuon->At(i);
       if(fabs(mu->Eta) < 2.5 && mu->PT > 20) muons.push_back(mu);
    }

    //Get jets
    std::vector<Jet*> jetz;
    for(unsigned int i=0;i<branchJet->GetEntries();i++){
       Jet *jet=(Jet*)branchJet->At(i);
       if(fabs(jet->Eta) < 2.8 && jet->PT > 20) jetz.push_back(jet);
    }

    int nElectrons = electrons.size();
    int nMuons     = muons.size();

    bool passSel=false;

    yields["NoSelection"]+=1.0;
    yieldsw["NoSelection"]+=weight;

    if(channel == "ee" && nElectrons == 2) passSel=true;
    else if(channel == "emu" && (nElectrons == 1 && nMuons == 1) ) passSel=true;
    else if(channel == "mumu" && nMuons ==2) passSel=true;

    if(!passSel) continue;

    if(channel == "ee"){ yields["2ELectrons"] += 1.0; yieldsw["2ELectrons"] += weight;}
    if(channel == "emu"){ yields["1ELec1Muon"] += 1.0; yieldsw["1Elec1Muon"] += weight;}
    if(channel == "mumu"){ yields["2Muons"] += 1.0; yieldsw["2Muons"] += weight;}

    //std::cout << "nElectrons= " << nElectrons << "nMuons= " << nMuons << std::endl;

    //std::cout << "nJets: " << branchJet->GetEntries() << std::endl;

    Electron *el1, *el2;
    Muon *mu1, *mu2;

    Int_t charge1 = -999;
    Int_t charge2 = -999;

    // If event contains at least 2 electrons
    if(channel == "ee"){
      el1 = (Electron *) electrons.at(0);
      el2 = (Electron *) electrons.at(1);

      charge1 = el1->Charge;
      charge2 = el2->Charge;
    }else if(channel == "emu"){
      el1 = (Electron *) electrons.at(0);
      mu1 = (Muon *) muons.at(0);

      charge1 = el1->Charge;
      charge2 = mu1->Charge;
    }else if(channel == "mumu"){
      mu1 = (Muon*) muons.at(0);
      mu2 = (Muon*) muons.at(1);

      charge1 = mu1->Charge;
      charge2 = mu2->Charge;
    }

    std::vector<int> charge;
    charge.push_back(charge1);
    charge.push_back(charge2);

    njets = jetz.size();
    //std::cout << "nJets: " << njets << std::endl;

    bool isChFlip = false;
    int initsign=charge.at(0)*charge.at(1);

    //int nElectrons = branchElectron->GetEntries();

    //Only do charge flip for electrons
    for(UInt_t iel=0;iel<nElectrons;iel++){
      rand = random->Rndm();
      Electron *el = (Electron *) electrons.at(iel);
      if(el->PT > 15 && el->PT<100){
        Int_t binx=lowPt->GetXaxis()->FindBin(el->Eta);
        float cfeff=lowPt->GetBinContent(binx);
        if(rand<cfeff) {
          std::cout << "Initial charge: " << charge[iel] << std::endl;
          charge[iel] *=-1;
          std::cout << "Charge flipped: " << charge[iel] << " elePt: " << el->PT << " eleEta: " << el->Eta << std::endl;
        }
      }else if (el->PT > 100 && el->PT < 200){
        Int_t binx=mediumPt->GetXaxis()->FindBin(el->Eta);
        float cfeff=mediumPt->GetBinContent(binx);
        if(rand<cfeff) {
          std::cout << "Initial charge: " << charge[iel] << std::endl;
          charge[iel] *=-1;
          std::cout << "Charge flipped: " << charge[iel] << " elePt: " << el->PT << " eleEta: " << el->Eta << std::endl;
        }
      }else if (el->PT > 200 && el->PT < 1000){
        Int_t binx=highPt->GetXaxis()->FindBin(el->Eta);
        float cfeff=highPt->GetBinContent(binx);
        if(rand<cfeff) {
          std::cout << "Initial charge: " << charge[iel] << std::endl;
          charge[iel] *=-1;
          std::cout << "Charge flipped: " << charge[iel] << " elePt: " << el->PT << " eleEta: " << el->Eta << std::endl;
        }
      }
    }

    int finalsign=charge.at(0)*charge.at(1);

    if(initsign != finalsign) isChFlip=true;

    if(charge.at(0)*charge.at(1)<0) continue;

    yields["SSLeptons"]  += 1.0;
    yieldsw["SSLeptons"] += weight;

    std::cout << "SS event selected " << weight << std::endl;

    //Get the leading jet
    if(njets!=0) jet=jetz.at(0);
    if(njets > 5) njets=5;

    float sfactor = 1.0;

    //Applying weight factor according to number of jets and channel
    if(channel == "ee" && apply_reweight == true) {
      sfactor *= sf_ee->GetBinContent(njets+1);
    }else if(channel == "emu" && apply_reweight == true){
      sfactor *= sf_ee->GetBinContent(njets+1)*sf_mumu->GetBinContent(njets+1)/2;
    }else if(channel == "mumu" && apply_reweight == true){
      sfactor *= sf_mumu->GetBinContent(njets+1);
    }

    //Only apply the factor for events with 2 or more jets
    if(njets >= 2 && (tag.find("signal") == std::string::npos) && apply_reweight == true) {
      weight *= sfactor;
      std::cout << "SF: " << sfactor << std::endl;
    }

    TH1Hists["hist_VR_nJets"].Fill(njets,weight);

    if(!isChFlip){
        if(njets>=4){
          if(channel=="ee"){
            TH1Hists["hist_SR_ElePT"].Fill(el1->PT,weight);
            TH1Hists["hist_SR_ElePT"].Fill(el2->PT,weight);
          }else if(channel=="emu"){
            TH1Hists["hist_SR_ElePT"].Fill(el1->PT,weight);
            TH1Hists["hist_SR_MuPT"].Fill(mu1->PT,weight);
          }else if(channel=="mumu"){
            TH1Hists["hist_SR_MuPT"].Fill(mu1->PT,weight);
            TH1Hists["hist_SR_MuPT"].Fill(mu2->PT,weight);
          }
          TH1Hists["hist_SR_JetPT"].Fill(jet->PT,weight);
        } else if(njets==0){
          if(channel=="ee"){
            TH1Hists["hist_CR_ElePT"].Fill(el1->PT,weight);
            TH1Hists["hist_CR_ElePT"].Fill(el2->PT,weight);
          }else if(channel=="emu"){
            TH1Hists["hist_CR_ElePT"].Fill(el1->PT,weight);
            TH1Hists["hist_CR_MuPT"].Fill(mu1->PT,weight);
          }else if(channel=="mumu"){
            TH1Hists["hist_CR_MuPT"].Fill(mu1->PT,weight);
            TH1Hists["hist_CR_MuPT"].Fill(mu2->PT,weight);
          }
        }
        if(channel=="ee"){
          TH1Hists["hist_VR_ElePT"].Fill(el1->PT,weight);
          TH1Hists["hist_VR_ElePT"].Fill(el2->PT,weight);
          TH1Hists["hist_VR_EleEta"].Fill(el1->Eta,weight);
          TH1Hists["hist_VR_EleEta"].Fill(el2->Eta,weight);
        }else if(channel=="emu"){
          TH1Hists["hist_VR_ElePT"].Fill(el1->PT,weight);
          TH1Hists["hist_VR_MuPT"].Fill(mu1->PT,weight);
          TH1Hists["hist_VR_EleEta"].Fill(el1->Eta,weight);
          TH1Hists["hist_VR_MuEta"].Fill(mu1->Eta,weight);
        }else if(channel=="mumu"){
          TH1Hists["hist_VR_MuPT"].Fill(mu1->PT,weight);
          TH1Hists["hist_VR_MuPT"].Fill(mu2->PT,weight);
          TH1Hists["hist_VR_MuEta"].Fill(mu1->Eta,weight);
          TH1Hists["hist_VR_MuEta"].Fill(mu2->Eta,weight);
        }
        if(njets!=0){ 
          TH1Hists["hist_VR_JetPT"].Fill(jet->PT,weight);
          TH1Hists["hist_VR_JetEta"].Fill(jet->Eta,weight);
        }
    }
    else {
        if(njets>=4){
          if(channel=="ee"){
            TH1Hists["hist_SR_ElePT_chflip"].Fill(el1->PT,weight);
            TH1Hists["hist_SR_ElePT_chflip"].Fill(el2->PT,weight);
            TH1Hists["hist_SR_JetPT_chflip"].Fill(jet->PT,weight);
          }else if(channel=="emu"){
            TH1Hists["hist_SR_ElePT_chflip"].Fill(el1->PT,weight);
            TH1Hists["hist_SR_MuPT_chflip"].Fill(mu1->PT,weight);
            TH1Hists["hist_SR_JetPT_chflip"].Fill(jet->PT,weight);
          }else if(channel=="mumu"){
            TH1Hists["hist_SR_MuPT_chflip"].Fill(mu1->PT,weight);
            TH1Hists["hist_SR_MuPT_chflip"].Fill(mu2->PT,weight);
            TH1Hists["hist_SR_JetPT_chflip"].Fill(jet->PT,weight);
          }
        } else if(njets==0){
          if(channel=="ee"){
            TH1Hists["hist_CR_ElePT_chflip"].Fill(el1->PT,weight);
            TH1Hists["hist_CR_ElePT_chflip"].Fill(el2->PT,weight);
          }else if(channel=="emu"){ 
            TH1Hists["hist_CR_ElePT_chflip"].Fill(el1->PT,weight);
            TH1Hists["hist_CR_MuPT_chflip"].Fill(mu1->PT,weight);
          }else if(channel=="mumu"){ 
            TH1Hists["hist_CR_MuPT_chflip"].Fill(mu1->PT,weight);
            TH1Hists["hist_CR_MuPT_chflip"].Fill(mu2->PT,weight);
          }
        }
        if(channel=="ee"){
            TH1Hists["hist_VR_ElePT_chflip"].Fill(el1->PT,weight);
            TH1Hists["hist_VR_ElePT_chflip"].Fill(el2->PT,weight);
        }else if(channel=="emu"){
            TH1Hists["hist_VR_ElePT_chflip"].Fill(el1->PT,weight);
            TH1Hists["hist_VR_MuPT_chflip"].Fill(mu1->PT,weight);
        }else if(channel=="mumu"){
            TH1Hists["hist_VR_MuPT_chflip"].Fill(mu1->PT,weight);
            TH1Hists["hist_VR_MuPT_chflip"].Fill(mu2->PT,weight);
        }

        if(njets!=0) TH1Hists["hist_VR_JetPT_chflip"].Fill(jet->PT,weight);
    }

    std::vector<TLorentzVector> leptoQuarks;
    TLorentzVector null4V(0,0,0,0);
    for(unsigned int i=0;i<6;i++) leptoQuarks.push_back(null4V);

    // Analyse missing ET
    if(branchMissingET->GetEntriesFast() > 0)
    {
      MissingET *met = (MissingET*) branchMissingET->At(0);
      //cout << met->MET << endl;
      if(njets>=4){
          TH1Hists["hist_SR_MET"].Fill(met->MET,weight);
      } else if(njets==0){
          TH1Hists["hist_CR_MET"].Fill(met->MET,weight);
      }
      TH1Hists["hist_VR_MET"].Fill(met->MET,weight);
    }

    if(njets<4) continue;

    yields["MoreThan3Jets"] += 1.0;
    yields["MoreThan3Jets"] += weight;

    //Form the jet combinations
    std::vector<TLorentzVector> jets;

    for(UInt_t icomb=0;icomb<3;icomb++)
    {
       std::vector<int> comb;
       comb.clear();
       if(icomb==0){comb.push_back(0);comb.push_back(1);comb.push_back(2);comb.push_back(3);}
       else if (icomb==1){comb.push_back(0);comb.push_back(2);comb.push_back(1);comb.push_back(3);}
       else if (icomb==2){comb.push_back(1);comb.push_back(3);comb.push_back(0);comb.push_back(2);}
       else if (icomb==3){comb.push_back(2);comb.push_back(3);comb.push_back(0);comb.push_back(1);}

       //std::cout << "Combination: " << comb.at(0) << " " << comb.at(1) << " " << comb.at(2) << " " << comb.at(3) << std::endl;

       jets.clear();

       for(UInt_t k = 0; k < jetz.size(); ++k)
       {
          if(k==4) break;
          // Form the different jet combinations
          Jet *jet = jetz.at(k);
          jets.push_back(jet->P4());
       }

       int index=1;
       Int_t ee_charge = 0;
       //Electron *elec1, *elec2;
       //Muon *muon1, *muon2;

       //std::cout << "seg fault is below this line" << electrons.size() << std::endl;

       if(channel=="ee")
       {
            // Take first two electrons
            //elec1 = (Electron *) electrons.at(0);
            //elec2 = (Electron *) electrons.at(1);

            TLorentzVector LQ;
            LQ=null4V;
            LQ += jets.at(0); LQ += jets.at(1); LQ += el1->P4();
            TH1Hists["hist_SR_LQMassE"].Fill(LQ.M(),weight);
            LQ=null4V;
            LQ += jets.at(2); LQ += jets.at(3); LQ += el2->P4();
            TH1Hists["hist_SR_LQMassE"].Fill(LQ.M(),weight);

       }else if(channel=="emu"){
            // Take now 1 electron and 1 muon
            //elec1 = (Electron *) electrons.at(0);
            //muon1 = (Muon *) muons.at(0);

            TLorentzVector LQ;
            LQ=null4V;
            LQ += jets.at(0); LQ += jets.at(1); LQ += el1->P4();
            TH1Hists["hist_SR_LQMassE"].Fill(LQ.M(),weight);
            LQ=null4V;
            LQ += jets.at(2); LQ += jets.at(3); LQ += mu1->P4();
            TH1Hists["hist_SR_LQMassMU"].Fill(LQ.M(),weight);

       }else if(channel=="mumu"){
            // Take now two muons
            //muon1 = (Muon *) muons.at(0);
            //muon2 = (Muon *) muons.at(1);

            TLorentzVector LQ;
            LQ=null4V;
            LQ += jets.at(0); LQ += jets.at(1); LQ += mu1->P4();
            TH1Hists["hist_SR_LQMassMU"].Fill(LQ.M(),weight);
            LQ=null4V;
            LQ += jets.at(2); LQ += jets.at(3); LQ += mu2->P4();
            TH1Hists["hist_SR_LQMassMU"].Fill(LQ.M(),weight);
       }


    }

    Jet *jet1 = jetz.at(0);
    Jet *jet2 = jetz.at(1);
    Jet *jet3 = jetz.at(2);
    Jet *jet4 = jetz.at(3);

    if(channel=="ee"){
      TLorentzVector ee = el1->P4()+el2->P4();
      TH1Hists["hist_SR_Mee"].Fill(ee.M());
      TLorentzVector eejj = el1->P4()+el2->P4()+jet1->P4()+jet2->P4()+jet3->P4()+jet4->P4();
      TH1Hists["hist_SR_eejjjjMass"].Fill(eejj.M(),weight);
      float H = el1->PT+el2->PT+jet1->PT+jet2->PT+jet3->PT+jet4->PT;
      TH1Hists["hist_SR_Heejjjj"].Fill(H, weight);
    }else if(channel=="emu"){
      TLorentzVector emujj = el1->P4()+mu->P4()+jet1->P4()+jet2->P4()+jet3->P4()+jet4->P4();
      TH1Hists["hist_SR_emujjjjMass"].Fill(emujj.M(),weight);
      float H = el1->PT+mu1->PT+jet1->PT+jet2->PT+jet3->PT+jet4->PT;
      TH1Hists["hist_SR_Hemujjjj"].Fill(H, weight);
    }else if(channel=="mumu"){
      TLorentzVector mumujj = mu1->P4()+mu2->P4()+jet1->P4()+jet2->P4()+jet3->P4()+jet4->P4();
      TH1Hists["hist_SR_mumujjjjMass"].Fill(mumujj.M(),weight);
      float H = mu1->PT+mu2->PT+jet1->PT+jet2->PT+jet3->PT+jet4->PT;
      TH1Hists["hist_SR_Hmumujjjj"].Fill(H, weight);
    }


  }

  std::map<std::string, TH1F>::iterator it = TH1Hists.begin();
        for( ; it != TH1Hists.end(); ++it ) it->second.Write();

  f0->Write();
  f0->Close();
  delete f;

  return;
}

void fillHists(std::string process,std::string pathFlag,std::string channel,bool hasMerge, int fileN)
{

  gROOT->SetBatch(1);
  gSystem->Load("libDelphes");

  gROOT->ProcessLine(".L loader.C+");

  std::string path="";

  const char* command=("mkdir "+outDir).c_str();
  gSystem->Exec(command);

  if(pathFlag=="path1") path = "/data/atlas/dbetalhc/";
  else if(pathFlag=="path0") path = "/data/user/e/edson/";

  std::vector<std::string> datasets;
  datasets.push_back(path+process+"/Events/run_01/tag_1_delphes_events.root");

  std::map< std::string, float > yields;
  std::map< std::string, float > yieldsw;

  std::map< std::string, float > initMap;
  initMap["NoSelection"]   = 0.0;
  initMap["2Electrons"]    = 0.0;
  initMap["1Elec1Muon"]    = 0.0;
  initMap["2Muons"]        = 0.0;
  initMap["SSLeptons"]     = 0.0;
  initMap["MoreThan3Jets"] = 0.0;

  yields  = initMap;
  yieldsw = initMap;

  for(UInt_t i=0;i<datasets.size();i++){
    std::cout << "Now running on process " << process << ", dataset:  " << datasets.at(i) << " and channel: " << channel << std::endl;
    if(process.find("zjj") == std::string::npos) processChain(datasets.at(i),process,channel,path,hasMerge,0,fileN,yields,yieldsw);
    else processChain(datasets.at(i),process,channel,path,hasMerge,30,fileN,yields,yieldsw);
  }


  std::map<std::string, float>::iterator it = yields.begin();
        for( ; it != yields.end(); ++it ){
     std::cout << it->first << " == " << it->second << std::endl;
  }

  std::map<std::string, float>::iterator it1 = yieldsw.begin();
        for( ; it1 != yieldsw.end(); ++it1 ){
     std::cout << it1->first << " == " << it1->second << std::endl;
  }

  exit();

}

float readCS(std::string inputFileName,int &nevents)
{
   std::cout << inputFileName << std::endl;
   std::ifstream inputFile(inputFileName.c_str());
   std::string line;
   std::string run_name="";
   std::string tag="";
   float cs = 0.000000;
   float error = 0.0;
   int nb_events = 0;

   float av_cs =0.0;
   int i=0;
   while(getline(inputFile, line)) {
      std::istringstream iss(line);
      iss>>run_name;
      iss>>tag;
      iss>>cs;
      iss>>error;
      iss>>nb_events;
      std::cout << "cs_read: " << cs << std::endl;
      if(cs!=0.0){
        av_cs+=cs;
        i++;
      }
      nevents+=nb_events;
   }
   av_cs=av_cs/i;
   std::cout << "av_cs: " << av_cs << " i: " << i << std::endl;
   return av_cs;
   //std::cout << "cs_read: " << cs << std::endl;
   //return cs;
}

//This is needed if parton shower merging with multiple files have been used
float readAvMergedCS(std::string inputFileName, int &nevents)
{
   std::cout << inputFileName << std::endl;
   std::ifstream inputFile(inputFileName.c_str());
   std::string line;
   std::string run_name="";
   std::string tag="";
   float cs = 0.0;
   float error = 0.0;
   int nb_events = 0;
   float cs_afterm = 0.0;
   int nEventsAM = 0;

   float av_cs =0.0;
   int i=0;
   while(getline(inputFile, line)) {
      std::istringstream iss(line);
      iss>>run_name;
      iss>>tag;
      iss>>cs;
      iss>>error;
      iss>>nb_events;
      iss>>cs_afterm;
      iss>>nEventsAM;
      std::cout << "cs_read: " << cs_afterm << std::endl;
      if(cs_afterm!=0.0){
        av_cs+=cs_afterm;
        i++;
      }
      nevents+=nEventsAM;
   }
   av_cs=av_cs/i;
   std::cout << "av_cs: " << av_cs << " i: " << i << std::endl;
   return av_cs;
}
