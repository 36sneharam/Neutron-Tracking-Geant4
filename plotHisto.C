// ROOT macro file for plotting example B4 histograms 
// 
// Can be run from ROOT session:
// root[0] .x plotHisto.C

{
  //gROOT->Reset();
//  gROOT->SetStyle("Plain");
  
  // Draw histos filled by Geant4 simulation 
  //   

  // Open file filled by Geant4 simulation 
  TFile f("HI_500MeV.root");
//  TFile g("HI_100MeV.root");
//  TFile h("HI_50MeV.root");

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  //c1->Divide(2,2);
  
  // Draw Eabs histogram in the pad 1
  //c1->cd(1);
  TH1D* h1_500 = (TH1D*)f.Get("EdepStep");
 // TH1D* h1_100 = (TH1D*)g.Get("EdepStep");
//  TH1D* h1_50 = (TH1D*)h.Get("EdepStep");
    
  h1_500->SetLineColor(kRed);
  h1_500->Draw("HIST");
  /*h1_100->SetLineColor(kBlue);
  h1_100->Draw("HISTsame");
  h1_50->SetLineColor(kBlack);
  h1_50->Draw("HISTsame");*/
    
    
  
 
}  
