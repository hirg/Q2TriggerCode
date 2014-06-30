#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TMath.h>
#include <TString.h>
#include <stdio.h> 
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>

#include <vector>
#include <map>

using namespace std;
//#include "l1ExtraTree.h"
//#include "l1Tree.h"

long makeKey(long run, long event){
  return (10000000*run + event);
}
//long makeKey(long run, long lumi, long event){
//  return (10000000000*run + 10000000*lumi + event);
//}

int getCentEtSum(int etsum)
{
    if(etsum > 1738){
      return 0;
    } else if(etsum > 1424 && etsum <= 1738){
      return 1;
    } else if(etsum > 599 && etsum <= 1424){
      return 2;
    } else if(etsum > 188 && etsum <= 599){
      return 3;
    } else if(etsum > 2 && etsum <= 188){
      return 4;
    } else if(etsum <= 2){
      return 5;
    } else{
      std::cout << "bad centrality determination" << std::endl;
      return -1;
    }
}

void matching_l12forest_Centq2(bool centPlot=false){
  //const TString l1_input = "/data/ginnocen/minbiasSkim_L1Tree_v1_Saturation/L1Tree_saturation.root";
  //const TString l1_input = "/data/ginnocen/L1Tree_18June_total.root";
  const TString l1_input = "~/Cent_q2_Trigger/CMSSW_7_1_0_pre8/src/L1Upgrade/L1UpgradeAnalyzer/test/L1UpgradeAnalyzer.root";


  TFile *lFile = TFile::Open(l1_input);
  //TTree *lTree = (TTree*)lFile->Get("l1ExtraTreeProducer/L1ExtraTree");
  //TTree *lEvtTree = (TTree*)lFile->Get("l1NtupleProducer/L1Tree");
  Int_t l1Up_evt, l1Up_run, l1Up_et, l1Up_q2;
  TTree *l1UpTree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  l1UpTree->SetBranchStatus("*",0);
  l1UpTree->SetBranchStatus("event",1);
  l1UpTree->SetBranchStatus("run",1);
  l1UpTree->SetBranchStatus("centrality_hwPt",1);
  l1UpTree->SetBranchStatus("v2_hwPt",1);

  l1UpTree->SetBranchAddress("event",&l1Up_evt);
  l1UpTree->SetBranchAddress("run",&l1Up_run);
  l1UpTree->SetBranchAddress("centrality_hwPt",&l1Up_et);
  l1UpTree->SetBranchAddress("v2_hwPt",&l1Up_q2);

  //l1ExtraTree *l1extra = new l1ExtraTree(lTree);
  //l1Tree *l1 = new l1Tree(lEvtTree);

  //const TString forest_input = "/data/richard/0.root";
  //const TString forest_input = "root://xrootd.cmsaf.mit.edu//store/user/luck/L1Emulator/minbiasForest_merged/0.root";
  const TString forest_input = "root://xrootd1.cmsaf.mit.edu//store/user/luck/L1Emulator/minbiasForest_merged/0.root";
  TFile *fFile = TFile::Open(forest_input);
  TTree *fTree = (TTree*)fFile->Get("akPu3CaloJetAnalyzer/t");
  TTree *fTowerTree = (TTree*)fFile->Get("rechitanalyzer/tower");
  fTree->AddFriend(fTowerTree);
  TTree *fEvtTree = (TTree*)fFile->Get("hiEvtAnalyzer/HiTree");
  fTree->AddFriend(fEvtTree);
  TTree *fSkimTree = (TTree*)fFile->Get("skimanalysis/HltTree");
  fTree->AddFriend(fSkimTree);

  Int_t f_evt, f_run, f_lumi, hiBin, f_n;
  Float_t f_eta[5000], f_phi[5000], f_et[5000];
  //std::vector<float> f_eta;
  //std::vector<float> f_phi;
  //std::vector<float> f_et;
  float hiHF;

  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("evt",1);
  fTree->SetBranchStatus("run",1);
  fTree->SetBranchStatus("lumi",1);
  fTree->SetBranchStatus("hiBin",1);
  fTree->SetBranchStatus("hiHF",1);
  fTree->SetBranchStatus("pcollisionEventSelection",1);
  fTree->SetBranchStatus("pHBHENoiseFilter",1);

  //fTowerTree->SetBranchStatus("*",0);
  //fTowerTree->SetBranchStatus("eta",1);
  //fTowerTree->SetBranchStatus("phi",1);
  //fTowerTree->SetBranchStatus("et",1);
  //fTowerTree->SetBranchStatus("n",1);
  fTree->SetBranchStatus("eta",1);
  fTree->SetBranchStatus("phi",1);
  fTree->SetBranchStatus("et",1);
  fTree->SetBranchStatus("n",1);


  fTree->SetBranchAddress("evt",&f_evt);
  fTree->SetBranchAddress("run",&f_run);
  fTree->SetBranchAddress("lumi",&f_lumi);
  fTree->SetBranchAddress("hiBin",&hiBin);
  fTree->SetBranchAddress("hiHF",&hiHF);

  //fTowerTree->SetBranchAddress("eta",&f_eta);
  //fTowerTree->SetBranchAddress("phi",&f_phi);
  //fTowerTree->SetBranchAddress("et",&f_et);
  //fTowerTree->SetBranchAddress("n",&f_n);
  fTree->SetBranchAddress("eta",&f_eta);
  fTree->SetBranchAddress("phi",&f_phi);
  fTree->SetBranchAddress("et",&f_et);
  fTree->SetBranchAddress("n",&f_n);
  
  Int_t pcollisionEventSelection, pHBHENoiseFilter;
  fTree->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  fTree->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  map<long, int> kmap;
  map<long, int> kmapcal;
  
  cout << "Filling the map..." << endl;
  int l1up_entries = l1UpTree->GetEntries();
  //int l_entries = lEvtTree->GetEntries();
  for(long j = 0; j < l1up_entries; ++j){
    if(j % 10000 == 0) printf("%ld / %d\n",j,l1up_entries);
    l1UpTree->GetEntry(j);
    //long key = makeKey(l1->run, l1->lumi, l1->event);
    long key = makeKey(l1Up_run, l1Up_evt);

    pair<long,int> p(key,j);
    kmap.insert(p);
    kmapcal.insert(p);
  }
  cout << "map filled!!!" << endl;

//Centrality histos

  TProfile *profileofflineCentralityVsl1Etsum_Calibration;
  TProfile *profilel1EtsumVsofflineCentrality_Calibration;
  if(centPlot){
    profileofflineCentralityVsl1Etsum_Calibration = new TProfile("profileofflineCentralityVsl1Etsum_Calibration","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",200,0,4000,0,200);
    profilel1EtsumVsofflineCentrality_Calibration = new TProfile("profilel1EtsumVsofflineCentrality_Calibration","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",50,0,200,0,4000);
  }
  
  int countCalib = 0;

//q2 histos
   TH1D* q2L1Hist   = new TH1D("q2L1","Online q2",500,0.0,1.0);
   TH1D* q2OffHist  = new TH1D("q2Off","Offline q2", 500, 0.0, 1.0);
   TH1D* q2OffHist_hiHF  = new TH1D("q2Off_hiHF","Offline q2", 500, 0.0, 1.0);
   TH2D* q2CorrHist = new TH2D("q2CorrHist","Corr between q2 online and q2 off",500,0.0,1.0,500,0.0,1.0);
   TH2D* q2CorrHist_hiHF = new TH2D("q2CorrHist_hiHF","Corr between q2 online and q2 off",500,0.0,1.0,500,0.0,1.0);

   TH2D* q2OffvsCentOffHist = new TH2D("q2OffvsCentOffHist","Corr between q2 off and cent off",1000,0,1000,500,0.0,1.0);
   TH2D* q2OffvsCentOffhiHFHist = new TH2D("q2OffvsCentOffhiHFHist","Corr between q2 off and cent off",1000,0,1000,500,0.0,1.0);

   TH1D* q2L1CentEtSumHist[6]; 
   TH1D* q2OffCentEtSumHist[6]; 
   TH2D* q2CorrCentEtsumHist[6]; 
   TProfile* q2CorrCentEtsumProf[6]; 

   TH1D* q2OffCenthiHFHist[6]; 
   TH2D* q2CorrCenthiHFHist[6]; 
   TProfile* q2CorrCenthiHFProf[6]; 

   //TH1D* q2L1CentOffHist[6]; 
   //TH1D* q2OffCentOffHist[6]; 
   //TH2D* q2CorrCentOffHist[6]; 

   for(int h=0; h<6; h++){
      q2L1CentEtSumHist[h] = new TH1D(Form("q2L1CentEtSumHist_%d",h),
                                      Form("q2 L1 for cent %d from etsum",h),
                                      500, 0., 1.); 
      q2OffCentEtSumHist[h] = new TH1D(Form("q2OffCentEtSumHist_%d",h),
                                       Form("q2 off for cent %d from etsum",h),
                                       500, 0., 1.); 
      q2CorrCentEtsumHist[h] = new TH2D(Form("q2CorrCentEtSumHist_%d",h),
                                        Form("q2 off for cent %d from etsum",h),
                                        500, 0., 1.,500, 0., 1.);
      q2CorrCentEtsumProf[h] = new TProfile(Form("q2CorrCentEtSumProf_%d",h),
                                        Form("q2 off for cent %d from etsum",h),
                                        500, 0., 1.);
                             
      q2OffCenthiHFHist[h] = new TH1D(Form("q2OffCenthiHFHist_%d",h),
                                       Form("q2 off for cent %d from etsum",h),
                                       500, 0., 1.); 
      q2CorrCenthiHFHist[h] = new TH2D(Form("q2CorrCenthiHFHist_%d",h),
                                        Form("q2 off for cent %d from etsum",h),
                                        500, 0., 1.,500, 0., 1.);
      q2CorrCenthiHFProf[h] = new TProfile(Form("q2CorrCenthiHFProf_%d",h),
                                        Form("q2 off for cent %d from etsum",h),
                                        500, 0., 1.);
     // q2L1CentOffHist[h] = new TH1D(Form("q2L1CentOffHist_%d",h),
     //                               Form("q2 L1 for cent %d from off",h),
     //                               500, 0., 1.); 
     // q2OffCentOffHist[h] = new TH1D(Form("q2OffCentOffHist_%d",h),
     //                                Form("q2 for cent %d from off",h),
     //                                500, 0., 1.); 
     // q2CorrCentOffHist[h] = new TH2D(Form("q2CorrCentOffHist_%d",h),
     //                                 Form("q2 off for cent %d from off",h),
     //                                 500, 0., 1.,500, 0., 1.); 
   } 

   Float_t q2Off  = 0.0;
   Float_t q2OffN = 0.0;
   Float_t q2OffN_hiHF = 0.0;
   Float_t q2On   = 0.0; 
   Float_t q2xOff = 0.0;
   Float_t q2yOff = 0.0;
   Float_t etOffSum = 0.0;

  int entries = fTree->GetEntries();
    
//  for(long int j = 1; j < entries; ++j){
  for(long int j = 1; j < 400000; ++j){
    if(j % 10000 == 0) printf("%ld / %d / %d\n",j,entries,countCalib);
    //if(j % 10 == 0) printf("%ld / %d\n",j,entries);
    fTowerTree->GetEntry(j);
    fTree->GetEntry(j);
    long keycal = makeKey(f_run, f_evt);
    //long keycal = makeKey(f_run, f_lumi, f_evt);
    
    map<long,int>::const_iterator gotcal = kmapcal.find(keycal);
    
    if (pcollisionEventSelection==0 || pHBHENoiseFilter ==0) continue;
    
    if(gotcal == kmapcal.end()){
      continue;      
    }
    else {
      l1UpTree->GetEntry(gotcal->second);
      //l1UpTree->GetEntry(kmap[keycal]);
      kmapcal.erase(keycal);
      //printf("id number %ld\n",j);
      //double etsum=l1extra->et[0];

      //Centrality part
      double etsum=l1Up_et;
      if(centPlot){
         profileofflineCentralityVsl1Etsum_Calibration->Fill(etsum,hiBin); 
         profilel1EtsumVsofflineCentrality_Calibration->Fill(hiBin,etsum); 
      }
      countCalib++;

      if(l1Up_et == 0) continue; 
      //q2 part
      for(long int l = 0; l < f_n; l++){
         if(f_eta[l] > -3. && f_eta[l] < 3.) continue;
         if(f_eta[l] < -5. || f_eta[l] > 5.) continue;

         q2xOff += f_et[l]*TMath::Cos(2*f_phi[l]);
         q2yOff += f_et[l]*TMath::Sin(2*f_phi[l]);

         etOffSum += f_et[l];
      }
                
      int centL1etsum = getCentEtSum(etOffSum);
      if(centL1etsum < 0) continue;
 
      q2Off = TMath::Sqrt(q2xOff*q2xOff + q2yOff*q2yOff);
      q2OffN = q2Off/etOffSum;
      q2OffN_hiHF = q2Off/hiHF;

      q2On = TMath::Sqrt((float)l1Up_q2)/(float)l1Up_et;	

      q2CorrHist->Fill(q2On,q2OffN);
      q2CorrHist_hiHF->Fill(q2On,q2OffN_hiHF);
      q2OffHist->Fill(q2OffN);
      q2OffHist_hiHF->Fill(q2OffN_hiHF);
      q2L1Hist->Fill(q2On);
     
      q2L1CentEtSumHist[centL1etsum]->Fill(q2On);
      q2OffCentEtSumHist[centL1etsum]->Fill(q2OffN); 
      q2OffCenthiHFHist[centL1etsum]->Fill(q2OffN_hiHF);
 
      q2CorrCentEtsumHist[centL1etsum]->Fill(q2On,q2OffN);
      q2CorrCentEtsumProf[centL1etsum]->Fill(q2On,q2OffN);
      
      q2CorrCenthiHFHist[centL1etsum]->Fill(q2On,q2OffN_hiHF);
      q2CorrCenthiHFProf[centL1etsum]->Fill(q2On,q2OffN_hiHF);

      //q2L1CentOffHist[]->Fill(q2On);
      //q2OffCentOffHist[]->Fill(q2OffN); 
      //q2CorrCentOffHist[]->Fill(q2On,q2OffN);

      //q2OffvsCentOffHist->Fill(hiBin,q2OffN);
      q2OffvsCentOffHist->Fill(etOffSum,q2OffN);
      q2OffvsCentOffhiHFHist->Fill(etOffSum,q2OffN_hiHF);

      q2Off  = 0.0;
      q2OffN = 0.0;
      q2On   = 0.0; 
      q2xOff = 0.0;
      q2yOff = 0.0;
      etOffSum = 0.0;
    }  
  }  

  printf("Matching entries: %d\n",countCalib);
  printf("step1\n");  
  printf("step2\n");  

  TH2D *hcorrl1EtsumVsofflineCentrality;
  TH2D *hcorrofflineCentralityVsl1Etsum;
  TH2D *hcorrL1CentralityVsfflineCentrality;
  TH1D *hofflineCentrality;
  TH1D *hl1Etsum;
  TH1D *hofflineEtsum;
  TH2D *hcorrOfflineEtsumVsL1Etsum;

  TProfile *profilel1EtsumVsofflineCentrality;
  TProfile *profileofflineCentralityVsl1Etsum;
  TProfile *profilel1CentralityVsofflineCentrality;

  TH1D *hOffline0to5;
  TH1D *hOffline5to10;
  TH1D *hOffline10to30;
  TH1D *hOffline30to50;
  TH1D *hOffline50to90;
  TH1D *hOffline90to100;
  
  TH1D *accepted[12][1];
  TH1D *fCenOffline[1];

  TF1 *fprofileofflineCentralityVsl1Etsum_Calibration;
  TF1 *fprofileofflinel1EtsumVsCentrality_Calibration;

  if(centPlot){
    hcorrl1EtsumVsofflineCentrality = new TH2D("hcorrl1EtsumVsofflineCentrality","L1 Etsum vs offline centrality; offline centrality ; L1 Etsum",100,0,200,100,0,4000);
    hcorrofflineCentralityVsl1Etsum = new TH2D("hcorrofflineCentralityVsl1Etsum","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",100,0,4000,100,0,200);  
    hcorrL1CentralityVsfflineCentrality = new TH2D("hcorrL1CentralityVsfflineCentrality","Online centrality vs offline centrality; L1 Etsum ; offline centrality",200,0,200,200,0,200);  
    hofflineCentrality = new TH1D("hofflineCentrality","Offline cent; Offline centrality; Entries",40,0.,200.);
    hl1Etsum = new TH1D("hl1Etsum","L1 Etsum ; L1 Etsum; Entries ",40,0.,4000.);
    hofflineEtsum = new TH1D("hofflineEtsum","Offline Etsum ; Offline Etsum; Entries ",40,0.,4000.);
    hcorrOfflineEtsumVsL1Etsum = new TH2D("hcorrOfflineEtsumVsL1Etsum","Offline Etsum vs L1 Etsum; Offline Etsum; L1 Etsum",100,0.,4000.,100,0.,4000.);
    
    profilel1EtsumVsofflineCentrality = new TProfile("profilel1EtsumVsofflineCentrality","L1 Etsum vs offline centrality; offline centrality ; L1 Etsum",100,0,200,0,4000);
    profileofflineCentralityVsl1Etsum = new TProfile("profileofflineCentralityVsl1Etsum","Offline centrality vs L1 Etsum; L1 Etsum ; offline centrality",400,0,4000,0,200);
    profilel1CentralityVsofflineCentrality = new TProfile("profilel1CentralityVsofflineCentrality","L1 centrality vs Offline Centrality; Offline Centrality ; L1 centrality",100,0,200,0,200);
  
    hOffline0to5 = new TH1D("hOffline0to5","hOffline0to5 ; Offline Centrality; Entries ",260,-30.,230.);
    hOffline5to10 = new TH1D("hOffline5to10","hOffline5to10 ; Offline Centrality; Entries ",260,-30.,230.);
    hOffline10to30 = new TH1D("hOffline10to30","hOffline10to30 ; Offline Centrality; Entries ",260,-30.,230.);
    hOffline30to50 = new TH1D("hOffline30to50","hOffline30to50 ; Offline Centrality; Entries ",260,-30.,230.);
    hOffline50to90 = new TH1D("hOffline50to90","hOffline50to90 ; Offline Centrality; Entries ",260,-30.,230.);
    hOffline90to100 = new TH1D("hOffline90to100","hOffline90to100 ; Offline Centrality; Entries ",260,-30.,230.);

    const int nBins = 100; 
    const double maxCen = 200.;

    accepted[0][0] = new TH1D("accepted_cen_0","",nBins,0,maxCen);
    accepted[1][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_10");
    accepted[2][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_20");
    accepted[3][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_30");
    accepted[4][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_40");
    accepted[5][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_50");
    accepted[6][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_100");
    accepted[7][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_130");
    accepted[8][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_160");
    accepted[9][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_170");
    accepted[10][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_180");
    accepted[11][0] = (TH1D*)accepted[0][0]->Clone("accepted_cen_190");
    
    fCenOffline[0] = new TH1D("fCenOffline",";Offline Centrality",nBins,0,maxCen);

    fprofileofflineCentralityVsl1Etsum_Calibration = new TF1("fprofileofflineCentralityVsl1Etsum_Calibration","pol9",0,2500);
    profileofflineCentralityVsl1Etsum_Calibration->Fit("fprofileofflineCentralityVsl1Etsum_Calibration");
    fprofileofflinel1EtsumVsCentrality_Calibration = new TF1("fprofileofflinel1EtsumVsCentrality_Calibration","pol9",0,200);
    profilel1EtsumVsofflineCentrality_Calibration->Fit("fprofileofflinel1EtsumVsCentrality_Calibration");
   
    const Double_t L1_THRESHOLD[12] = {0.,10.,20.,30.,40.,50.,100.,130.,160.,170.,180.,190.};
    double L1centrality=0.;
        
    for(long int j = 1; j < entries; ++j){
    
      if(j % 100000 == 0) printf("%ld / %d\n",j,entries);
      fTree->GetEntry(j);
      if (pcollisionEventSelection==0 || pHBHENoiseFilter ==0) continue;
  
    //  long key = makeKey(f_run, f_lumi, f_evt);
      long key = makeKey(f_run, f_evt);
  
      map<long,int>::const_iterator got = kmap.find(key);
      if(got == kmap.end()){
        continue;      
      }
      else {
       l1UpTree->GetEntry(got->second);  
       kmap.erase(key);
       //printf("id number %ld\n",j);
       //double etsum=l1extra->et[0];
       double etsum=l1Up_et;
  
       //printf("et value %f\n",etsum);
       //printf("hiBin %d\n",hiBin);
       //printf("hiHF %f\n",hiHF);
            
       hcorrl1EtsumVsofflineCentrality->Fill(hiBin,etsum);
       hcorrofflineCentralityVsl1Etsum->Fill(etsum,hiBin);
       hofflineCentrality->Fill(hiBin);
       hl1Etsum->Fill(etsum);
       hofflineEtsum->Fill(hiHF);
       hcorrOfflineEtsumVsL1Etsum->Fill(hiHF,etsum); 
       profileofflineCentralityVsl1Etsum->Fill(etsum,hiBin);  
       profilel1EtsumVsofflineCentrality->Fill(hiBin,etsum);  	  
  
       L1centrality=fprofileofflinel1EtsumVsCentrality_Calibration->GetX(etsum);
       printf("L1centrality %f\n",L1centrality);
       hcorrL1CentralityVsfflineCentrality->Fill(hiBin,L1centrality);
       profilel1CentralityVsofflineCentrality->Fill(hiBin,L1centrality);
            
       double Etsum_cen5=fprofileofflinel1EtsumVsCentrality_Calibration->Eval(10.);
       double Etsum_cen10=fprofileofflinel1EtsumVsCentrality_Calibration->Eval(20.);
       double Etsum_cen30=fprofileofflinel1EtsumVsCentrality_Calibration->Eval(60.);
       double Etsum_cen50=fprofileofflinel1EtsumVsCentrality_Calibration->Eval(100.);
       double Etsum_cen90=fprofileofflinel1EtsumVsCentrality_Calibration->Eval(180.);
  
       if(etsum>Etsum_cen5) hOffline0to5->Fill(hiBin);
       if(etsum>Etsum_cen10 && etsum<Etsum_cen5) hOffline5to10->Fill(hiBin);
       if(etsum>Etsum_cen30 && etsum<Etsum_cen10) hOffline10to30->Fill(hiBin);
       if(etsum>Etsum_cen50 && etsum<Etsum_cen30) hOffline30to50->Fill(hiBin);
       if(etsum>Etsum_cen90 && etsum<Etsum_cen50) hOffline50to90->Fill(hiBin);
       if(etsum<Etsum_cen90) hOffline90to100->Fill(hiBin);
  
       fCenOffline[0]->Fill(hiBin);
           
       for(int k = 0; k < 12; ++k){
           if(L1centrality>L1_THRESHOLD[k]) accepted[k][0]->Fill(hiBin);
       }
      }  
    }  
  }

  //TFile*fouput=new TFile("CalibrationEtsum_18June.root","recreate"); 
  TFile*fouput=new TFile("CalibrationEtsum_q2_27June.root","recreate"); 
  fouput->cd();

  if(centPlot){
  
    //Centrality part  
    hcorrl1EtsumVsofflineCentrality->Write();
    hcorrofflineCentralityVsl1Etsum->Write();
    hofflineCentrality->Write();
    hl1Etsum->Write();
    hofflineEtsum->Write();
    hcorrOfflineEtsumVsL1Etsum->Write();
    profilel1EtsumVsofflineCentrality->Write();
    profileofflineCentralityVsl1Etsum->Write();
    profileofflineCentralityVsl1Etsum_Calibration->Write();
    fprofileofflineCentralityVsl1Etsum_Calibration->SetName("fprofileofflineCentralityVsl1Etsum_Calibration");
    fprofileofflineCentralityVsl1Etsum_Calibration->Write();
    fprofileofflinel1EtsumVsCentrality_Calibration->SetName("fprofileofflinel1EtsumVsCentrality_Calibration");
    fprofileofflinel1EtsumVsCentrality_Calibration->Write();
    hcorrL1CentralityVsfflineCentrality->Write();
    profilel1CentralityVsofflineCentrality->Write();
    hOffline0to5->Write();
    hOffline5to10->Write();
    hOffline10to30->Write();
    hOffline30to50->Write();
    hOffline50to90->Write();
    hOffline90to100->Write();
    
  
    TGraphAsymmErrors *a[12][1];
    for(int k = 0; k < 12; ++k)
    {
      for(int l = 0; l < 1; ++l)
      {
        a[k][l] = new TGraphAsymmErrors();
        a[k][l]->BayesDivide(accepted[k][l],fCenOffline[l]);
      }
    }
  
    a[0][0]->SetName("asymm_cen_0_cen");
    a[1][0]->SetName("asymm_cen_10_cen");
    a[2][0]->SetName("asymm_cen_20_cen");
    a[3][0]->SetName("asymm_cen_30_cen");
    a[4][0]->SetName("asymm_cen_40_cen");
    a[5][0]->SetName("asymm_cen_50_cen");
    a[6][0]->SetName("asymm_cen_100_cen");
    a[7][0]->SetName("asymm_cen_130_cen");
    a[8][0]->SetName("asymm_cen_160_cen");
    a[9][0]->SetName("asymm_cen_170_cen");
    a[10][0]->SetName("asymm_cen_180_cen");
    a[11][0]->SetName("asymm_cen_190_cen");
    
    fCenOffline[0]->Write();
    
    for(int k = 0; k < 12; ++k){
      for(int l = 0; l < 1; ++l)
      {
        accepted[k][l]->Write();
        a[k][l]->Write();
      }
    }
  }


//q2 part  
  q2CorrHist->Write();
  q2CorrHist_hiHF->Write();
  q2OffHist->Write();
  q2OffHist_hiHF->Write();
  q2L1Hist->Write();   
  q2OffvsCentOffHist->Write();
  q2OffvsCentOffhiHFHist->Write();

  q2L1CentEtSumHist[0]->Write(); 
  q2L1CentEtSumHist[1]->Write(); 
  q2L1CentEtSumHist[2]->Write(); 
  q2L1CentEtSumHist[3]->Write(); 
  q2L1CentEtSumHist[4]->Write(); 
  q2L1CentEtSumHist[5]->Write(); 

  q2OffCentEtSumHist[0]->Write(); 
  q2OffCentEtSumHist[1]->Write(); 
  q2OffCentEtSumHist[2]->Write(); 
  q2OffCentEtSumHist[3]->Write(); 
  q2OffCentEtSumHist[4]->Write(); 
  q2OffCentEtSumHist[5]->Write(); 

  q2CorrCentEtsumHist[0]->Write(); 
  q2CorrCentEtsumHist[1]->Write(); 
  q2CorrCentEtsumHist[2]->Write(); 
  q2CorrCentEtsumHist[3]->Write(); 
  q2CorrCentEtsumHist[4]->Write(); 
  q2CorrCentEtsumHist[5]->Write(); 

  q2CorrCentEtsumProf[0]->Write(); 
  q2CorrCentEtsumProf[1]->Write(); 
  q2CorrCentEtsumProf[2]->Write(); 
  q2CorrCentEtsumProf[3]->Write(); 
  q2CorrCentEtsumProf[4]->Write(); 
  q2CorrCentEtsumProf[5]->Write(); 

  q2OffCenthiHFHist[0]->Write(); 
  q2OffCenthiHFHist[1]->Write(); 
  q2OffCenthiHFHist[2]->Write(); 
  q2OffCenthiHFHist[3]->Write(); 
  q2OffCenthiHFHist[4]->Write(); 
  q2OffCenthiHFHist[5]->Write(); 

  q2CorrCenthiHFHist[0]->Write(); 
  q2CorrCenthiHFHist[1]->Write(); 
  q2CorrCenthiHFHist[2]->Write(); 
  q2CorrCenthiHFHist[3]->Write(); 
  q2CorrCenthiHFHist[4]->Write(); 
  q2CorrCenthiHFHist[5]->Write(); 

  q2CorrCenthiHFProf[0]->Write(); 
  q2CorrCenthiHFProf[1]->Write(); 
  q2CorrCenthiHFProf[2]->Write(); 
  q2CorrCenthiHFProf[3]->Write(); 


  fouput->Close();
  lFile->Close();
  fFile->Close();
}

int main(){
  matching_l12forest_Centq2(false);
  return 0; 
}
