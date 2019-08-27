// Builds a graph with errors, displays it and saves it as
// image. First, include some header files
// (not necessary for Cling)

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TFile.h"

void macro_Al27(){
    TFile *MyFile = new TFile("B4-Calor.root","READ");
    //MyFile->ls();
    TH1F *h0 = (TH1F*)MyFile->Get("Interact_Num");
    TH1F *hist = (TH1F*)MyFile->Get("ElasvsInelas");
    //hist->StyleLineWidth(10);
    
    //h0->Draw();
    
    //Getting Bin Content of each column; 
    
    auto E1 = h0->GetBinContent(1); // Interaction_num 10
    auto InE2 = h0->GetBinContent(2); //Interact_Num 11
    auto InE3 = h0->GetBinContent(3); // Interact_Num 12
    auto InE4 = h0->GetBinContent(4); // Interact_Num 13
    //auto InE5 = h0->GetBinContent(5); // Interact_Num 14

    
    
    cout<<"10: "<<E1<<endl;
    cout<<"11: "<<InE2<<endl;
    cout<<"12: "<<InE3<<endl;
    cout<<"13: "<<InE4<<endl;
    //cout<<"14: "<<InE5<<endl;

    
    
   
    
    //Input each bin content into its own histogram 
    //Elastic
    TH1F *h1 = new TH1F("h1", "h1 title", 2, 1, 3);
    h1->AddBinContent(1, E1);
    h1->SetFillColor(kBlue-4);  // Interact_Num 10
    //Inelastic 
    TH1F *h2 = new TH1F("h2", "h2 title", 2, 1, 3);
    h2->AddBinContent(2, InE2);
    h2->SetFillColor(kGreen+1);  // Interact_Num 11
    TH1F *h3 = new TH1F("h3", "h3 title", 2, 1, 3);
    h3->AddBinContent(2, InE3);
    h3->SetFillColor(kOrange+7); // Interact_Num 12
    TH1F *h4 = new TH1F("h4", "h4 title", 2, 1, 3);
    h4->AddBinContent(2, InE4);
    h4->SetFillColor(kViolet+3); // Interact_Num 13
    
    
    
   
    
   
    
    //Create Stack 
    THStack *hs = new THStack("InEvsE_manual", "Inelastic vs Elastic validation");
    
    hs->Add(h1);
    hs->Add(h2);
    hs->Add(h3);
    hs->Add(h4); 
    //hs->Add(h5); 
    
    
    
    
    hist -> SetLineWidth(3);
    
    auto El_total= hist->GetBinContent(1); 
    auto InEl_total= hist->GetBinContent(2); 
    
    cout<<"Inelastic Total: "<<InEl_total<<endl; 
    cout<<"Elastic Total: "<<El_total<<endl; 
    
    THStack *hists = new THStack("Stack", "Inelastic vs Elastic validation");
    
    hists->Add(hist);

    //h0->Draw(); 
    hs->Draw();
    hists->Draw("SAME");
    
    
    double legxtranslate = 0.1;
    double legytranslate = 0; 
    
    auto legend = new TLegend(0.1+legxtranslate,0.7,0.48,0.9);
   //legend->SetHeader("Legend"); // option "C" allows to center the header
   
   legend->AddEntry(h1,"Al27 no gamma");
   legend->AddEntry(h2,"Al27 gamma");
   legend->AddEntry(h3,"Na + alpha + gamma");
   legend->AddEntry(h4,"Al isotope");
   //legend->AddEntry(h5,"Triton + alpha");
   
   
   legend->Draw();
   
    
    
    /*T
    H1F *h1 = new TH1F("h1", "h1 title", 10, 0, 10);
    h1->FillRandom("gaus", 100); 
    h1->SetFillColor(kRed); 
    TH1F *h2 = new TH1F("h2", "h2 title", 10, 0, 10);
    h2->FillRandom("gaus", 500);
    h2->SetFillColor(kBlue);
    THStack *hs = new THStack("hs", "hs title");
    hs->Add(h1); 
    hs->Add(h2); 
    
    TH1F *h3 = new TH1F("h3", "h3 title", 10, 0, 100);
    h3 ->Fill(h1->GetBinContent(1)); 
    
    h3 ->Draw();
    
    cout<<h1->GetBinContent(1)<<endl;
    */
    
    
}


//Get bin content 



