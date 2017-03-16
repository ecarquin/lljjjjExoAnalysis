#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TASImage.h"
#include "TFile.h"
	
void chargeFlipHist()
{

    gStyle->SetOptStat(0);

    TASImage image("chargeFlipPlot.png");
    TFile *f=TFile::Open("chargeFlipAtlas.root","NEW");

    float etaBins[8] = {0.0,0.6,1.1,1.36,1.5,1.8,2.3,2.6};

    TH1F *h1 = new TH1F("h1","Charge flip electrons ATLAS",7,etaBins);
    TH1F *h2 = new TH1F("h2","Charge flip electrons ATLAS",7,etaBins);
    TH1F *h3 = new TH1F("h3","Charge flip electrons ATLAS",7,etaBins);
    TH1F *h4 = new TH1F("h4","Eta Binning",7,etaBins);

    h4->Fill(h4->GetBinCenter(1),1);
    h4->Fill(h4->GetBinCenter(2),2);
    h4->Fill(h4->GetBinCenter(3),3);
    h4->Fill(h4->GetBinCenter(4),4);
    h4->Fill(h4->GetBinCenter(5),5);
    h4->Fill(h4->GetBinCenter(6),6);
    h4->Fill(h4->GetBinCenter(7),7);    

    h1->Fill(h1->GetBinCenter(1),0.00056);
    h1->Fill(h1->GetBinCenter(2),0.00098);
    h1->Fill(h1->GetBinCenter(3),0.002);
    //h1->Fill(h1->GetBinCenter(4),0.0);
    h1->Fill(h1->GetBinCenter(5),0.009);
    h1->Fill(h1->GetBinCenter(6),0.017);
    h1->Fill(h1->GetBinCenter(7),0.0395);

    h2->Fill(h2->GetBinCenter(1),0.00165);
    h2->Fill(h2->GetBinCenter(2),0.0027);
    h2->Fill(h2->GetBinCenter(3),0.004);
    //h2->Fill(h2->GetBinCenter(4),0.0);
    h2->Fill(h2->GetBinCenter(5),0.022);
    h2->Fill(h2->GetBinCenter(6),0.035);
    h2->Fill(h2->GetBinCenter(7),0.072);

    h3->Fill(h3->GetBinCenter(1),0.0037);
    h3->Fill(h3->GetBinCenter(2),0.0027);
    h3->Fill(h3->GetBinCenter(3),0.0157);
    //h2->Fill(h2->GetBinCenter(4),0.0);
    h3->Fill(h3->GetBinCenter(5),0.052);
    h3->Fill(h3->GetBinCenter(6),0.074);
    h3->Fill(h3->GetBinCenter(7),0.12);

    TCanvas *c=new TCanvas("c");

    //TPad *p0 = new TPad("p0","p",0,0,1,1);
    //p0->Draw();
    image.Draw("x");

    //Create a transparent pad filling the full canvas
    TPad *p = new TPad("p","p",0,0,1,1);
    p->SetFillStyle(4000);
    p->SetFrameFillStyle(4000);
    p->SetFrameBorderMode(0);
    p->Draw();
    p->cd();
    p->SetLogy();

    h1->GetYaxis()->SetRangeUser(9e-5,1);
    h1->SetTitle("");

    h1->SetMarkerStyle(4);
    h1->SetMarkerColor(1);

    h2->SetMarkerStyle(25);
    h2->SetMarkerColor(2); 

    h3->SetMarkerStyle(26);
    h3->SetMarkerColor(4);   

    h1->Draw("P");
    h2->Draw("same:P");
    h3->Draw("same:P");

    c->SaveAs("chargeFlipAtlas.png");

    c->Clear();
    c->Modified();

    h4->Draw();
    c->SaveAs("etaBinning.png");

    h1->Write("lowPt");
    h2->Write("mediumPt");
    h3->Write("highPt");
    h4->Write("etaBin");

}

