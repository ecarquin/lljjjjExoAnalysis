#include "TFile.h"
#include "TH1F.h"

#include <iostream>

void getSF()
{
  TFile *f=TFile::Open("09_03_17/GetScaleFactor.root","READ");
  //TFile *fm=TFile::Open("31_08_16/GetScaleFactor_mumu.root","READ");
  TFile *fatlase = TFile::Open("chargeFlipAtlas_6c.root","READ");
  //TFile *fatlasm = TFile::Open("chargeFlipAtlas_6d.root","READ");

  TFile *fout = TFile::Open("SF.root","RECREATE");
  TH1F* sf_ee=new TH1F("sf_ee","",6,-0.5,5.5);
  TH1F* sf_mumu=new TH1F("sf_mumu","",6,-0.5,5.5);

  TH1F* nJets_ee=(TH1F*)f->Get("hist_VR_nJets_ee");
  TH1F* nJets_mumu=(TH1F*)f->Get("hist_VR_nJets_mumu");

  TH1F* nJets_ee_atlas = (TH1F*)fatlase->Get("lowPt");
  //TH1F* nJets_mumu_atlas = (TH1F*)fatlasm->Get("lowPt");

  float lumi_frac=20.3/100.;
  float sqrts_frac=155.998/505.582;

  float totalFact = lumi_frac*sqrts_frac;

  std::cout << "Lumi factor: " << lumi_frac << std::endl;
  std::cout << "Sqrts factor: " << sqrts_frac << std::endl;
  std::cout << "Total factor: " << totalFact << std::endl;

  int entries_ee[6]={0};
  int entries_atlas_ee[6]={0};

  float entries_mumu[6]={0};
  float entries_atlas_mumu[7]={0.0,215.,145.,100.,50.,26.,25.};

  for(UInt_t ibin=0;ibin< sf_ee->GetNbinsX()+1;ibin++){
    std::cout << "Entries: " << fabs(nJets_ee->GetBinContent(ibin)) << " " << nJets_ee_atlas->GetBinContent(ibin) << std::endl;
    std::cout << "Entries: " << fabs(nJets_mumu->GetBinContent(ibin)) << " " << entries_atlas_mumu[ibin] << std::endl;
    entries_ee[ibin] = nJets_ee->GetBinContent(ibin);
    entries_atlas_ee[ibin] = nJets_ee_atlas->GetBinContent(ibin);
    entries_mumu[ibin] = nJets_mumu->GetBinContent(ibin);
  }

  nJets_ee->Scale(totalFact);
  nJets_mumu->Scale(totalFact);

  float weights_ee[6]={0};
  float weights_mumu[6]={0};

  for(UInt_t ibin=1;ibin< sf_ee->GetNbinsX()+1;ibin++){
    std::cout << "Entries after scaling: " << nJets_ee->GetBinContent(ibin) << " " << nJets_ee_atlas->GetBinContent(ibin) << std::endl;
    weights_ee[ibin] = nJets_ee_atlas->GetBinContent(ibin)/nJets_ee->GetBinContent(ibin);
    //weights_ee[ibin] = (nJets_ee->GetBinContent(ibin)-nJets_ee_atlas->GetBinContent(ibin))/nJets_ee->GetBinContent(ibin);
    std::cout << "Entries after scaling: " << nJets_mumu->GetBinContent(ibin) << " " << entries_atlas_mumu[ibin] << std::endl;
    weights_mumu[ibin] = entries_atlas_mumu[ibin]/nJets_mumu->GetBinContent(ibin);
    //weights_mumu[ibin] = (nJets_mumu->GetBinContent(ibin)-entries_atlas_mumu[ibin])/nJets_mumu->GetBinContent(ibin);

    std::cout << "Weights ee: " << weights_ee[ibin] << " Weights mumu: " << weights_mumu[ibin] << std::endl;
    std::cout << "Bin center: " << sf_ee->GetBinCenter(ibin) << " ibin: " << ibin << std::endl;
    sf_ee->SetBinContent(ibin,1.+fabs(weights_ee[ibin]));
    sf_mumu->SetBinContent(ibin,1.+fabs(weights_mumu[ibin]));
  }
	
  sf_ee->Write();
  sf_mumu->Write();

  

}
