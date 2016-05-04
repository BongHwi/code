// Analysis Task for the Quality Assurence of Beam Gas Monitoring
//
// This code will check the each event for several parameters,
// and check if it is Background or Signal with below function.
// after that,
//
// Authors
// Alexander Borissov <aborisso@mail.cern.ch>
// Bong-Hwi Lim <bong-hwi.lim@cern.ch>
// Jihye Song <Jihye.Song@cern.ch>
//
// If you have any comment or question of this code,
// Please send a mail to Bong-Hwi or Jihye
//
// Last update: 2016.04.09 (blim)
//
//#include <Riostream.h>
#include <iostream>
#include"AliAnalysisBGMonitorQA.h"
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH1D.h>
#include"TCanvas.h"
#include"TArrayI.h"
#include"AliAnalysisTaskSE.h"
#include"AliAnalysisManager.h"
#include"AliESD.h"
#include"AliESDEvent.h"
#include"AliESDfriend.h"
#include"AliVEvent.h"
#include"AliESDInputHandler.h"
#include"AliLog.h"
#include"AliAnalysisFilter.h"
#include"AliESDtrackCuts.h"
#include"AliESDVertex.h"
#include"AliESDtrack.h"
#include"AliTriggerAnalysis.h"
#include"AliAnalysisCuts.h"
#include"AliMultiplicity.h"
#include"AliESDVZERO.h"
#include"AliESDVZEROfriend.h"
#include"AliESDTZERO.h"
#include"AliAnalysisUtils.h"
#include"AliESDAD.h"
#include"AliESDADfriend.h"
using namespace std;
ClassImp(AliAnalysisBGMonitorQA)

Bool_t IsItBGSPDClusterVsTracklet(AliVEvent *event); // add function info and initial condition (blim)
Bool_t IsItBGSPDClusterVsTracklet2(AliVEvent *event); // add function to check bg in small tracklet 2015.09.14. (blim)

//________________________________________________________________________
AliAnalysisBGMonitorQA::AliAnalysisBGMonitorQA(const char *name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fESDfriend(0x0),
fTreeTrack(0),
fTreeTrack2(0),
fList(0),
fList2(0),
fUseTree(kFALSE),
runNumber(0),
fvertZ(0),
fvertX(0),
fvertY(0),
fvertTPCZ(0),
fvertTPCX(0),
fvertTPCY(0),
fvertZ2(0),
fvertX2(0),
fvertY2(0),
fvertTPCX2(0),
fvertTPCY2(0),
fvertTPCZ2(0),
fv0a(0),
fv0c(0),
fad0a(0),
fad0c(0),
fMulta(0),
fMultc(0),
fTriCha(0),
fTriChc(0),
fV0M(0),
fbx(0),
ftime(0),
fSpdC1(0),
fSpdC2(0),
fSpdT(0),
ntracks(0),
V0A(0),
V0C(0),
V0ABG(0),
V0CBG(0),
VBA(0),
VBC(0),
VGA(0),
VGC(0),
VTX(0),
bgID(0),
bgID2(0),
t0PileUp(0),
spdPileUp(0),
spdPileUpOutOfBunch(0),
triMask(0),
fastORHW(0),
SPD1(0),
SPD2(0),
SPDHw1(0),
SPDHw2(0),
ntr(0),
nbunch(0),
nV0A(0),
nV0C(0),
nV0ABG(0),
nV0CBG(0) // add initiallize 2016.03.31. (blim)n

{
    // Constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class()); //CINT7
    DefineOutput(2, TList::Class()); //V0M, SH2
    DefineOutput(0, TTree::Class()); //RunNumber
    
}

//________________________________________________________________________
void AliAnalysisBGMonitorQA::ConnectInputData(Option_t *)
{
    
    TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
    if (!tree) {
        Printf("ERROR: Could not read chain from input slot 0");
    } else {
        
        AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
        
        if (esdH) {
            fESD = (AliESDEvent*) esdH->GetEvent();
            if(fESD) {
                fESDfriend = (AliESDfriend*)fESD->FindListObject("AliESDfriend");
                if (!fESDfriend){
                    AliError("No friend found");
                }
            }
        } else {
            Printf("ERROR: Could not get ESDInputHandler");
        }
        
    }
}

//________________________________________________________________________
void AliAnalysisBGMonitorQA::CreateOutputObjects()
{

    fTreeTrack2 = new TTree("TreeTrack","Track Properties2");
    fTreeTrack2->Branch("runNumber",&runNumber,"runNumber/I"); //run number
    PostData(0, fTreeTrack2);

    if(fList != NULL){
        delete fList;
        fList = NULL;
    }
    if(fList2 != NULL){
        delete fList2;
        fList2 = NULL;
    }
    if(fList == NULL){
        fList = new TList();
        fList->SetOwner(kTRUE);
    }
    if(fList2 == NULL){
        fList2 = new TList();
        fList2->SetOwner(kTRUE);
    }
    
    //__________CINT7__________
    TH1F *hNumEffPurityBC;
    TH1F *hDenomEffBC;
    TH1F *hDenomPurityBC;
    TH1F *hDenomRejecEffBC;
    TH1F *hNumRejecEffBC;
    TH1F *hSPDNumBC;
    TH1F *hSPDDenomBC;
    TH2F *hNumTrkVsClsSPID;
    TH2F *hDenomTrkVsClsSPID;
    TH2F *hNumV0;
    TH2F *hDenomV0;

    hNumEffPurityBC = new TH1F(Form("hNumEffPurityBC"),"; #V0flags in PF", 35, 0, 35);
    hDenomEffBC = new TH1F(Form("hDenomEffBC"),"; #V0flags in PF", 35, 0, 35);
    hDenomPurityBC = new TH1F(Form("hDenomPurityBC"),"; #V0flags in PF", 35, 0, 35);
    hDenomRejecEffBC = new TH1F(Form("hDenomRejecEffBC"),"; #V0flags in PF", 35, 0, 35);
    hNumRejecEffBC = new TH1F(Form("hNumRejecEffBC"),"; #V0flags in PF", 35, 0, 35);

    hSPDNumBC = new TH1F(Form("hSPDNumBC"),"; Spd tracklet", 200, 0, 200);
    hSPDDenomBC = new TH1F(Form("hSPDDenomBC"),"; Spd tracklet", 200, 0, 200);

    hNumTrkVsClsSPID = new TH2F(Form("hNumTrkVsClsSPID"),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
    hNumTrkVsClsSPID->GetXaxis()->SetTitle("Tracklet");
    hNumTrkVsClsSPID->GetYaxis()->SetTitle("Cluster (fspdC1)");
    hDenomTrkVsClsSPID= new TH2F(Form("hDenomTrkVsClsSPID"),"; Spd : !BGid",140,0,140,500,0,500);
    hDenomTrkVsClsSPID->GetXaxis()->SetTitle("Tracklet");
    hDenomTrkVsClsSPID->GetYaxis()->SetTitle("Cluster (fspdC1)");

    hNumV0 = new TH2F(Form("hNumV0"),"; V0 : !BGid & GoodEvent",600,-300,300,2000,-1000,1000);
    hNumV0->GetXaxis()->SetTitle("V0A-V0C");
    hNumV0->GetYaxis()->SetTitle("V0A+V0C");
    hDenomV0= new TH2F(Form("hDenomV0"),"; V0 : !BGid",600,-300,300,2000,-1000,1000);
    hDenomV0->GetXaxis()->SetTitle("V0A-V0C");
    hDenomV0->GetYaxis()->SetTitle("V0A+V0C");

    fList->Add(hNumEffPurityBC);
    fList->Add(hDenomEffBC);
    fList->Add(hDenomPurityBC);
    fList->Add(hDenomRejecEffBC);
    fList->Add(hNumRejecEffBC);
    fList->Add(hSPDNumBC);
    fList->Add(hSPDDenomBC);
    fList->Add(hNumTrkVsClsSPID);
    fList->Add(hDenomTrkVsClsSPID);
    fList->Add(hNumV0);
    fList->Add(hDenomV0);

    //__________V0M__________
    TH1F *hNumEffPurityBC_V0M;
    TH1F *hDenomEffBC_V0M;
    TH1F *hDenomPurityBC_V0M;
    TH1F *hDenomRejecEffBC_V0M;
    TH1F *hNumRejecEffBC_V0M;
    TH1F *hSPDNumBC_V0M;
    TH1F *hSPDDenomBC_V0M;
    TH2F *hNumTrkVsClsSPID_V0M;
    TH2F *hDenomTrkVsClsSPID_V0M;
    TH2F *hNumV0_V0M;
    TH2F *hDenomV0_V0M;

    hNumEffPurityBC_V0M = new TH1F(Form("hNumEffPurityBC_V0M"),"; #V0flags in PF", 35, 0, 35);
    hDenomEffBC_V0M = new TH1F(Form("hDenomEffBC_V0M"),"; #V0flags in PF", 35, 0, 35);
    hDenomPurityBC_V0M = new TH1F(Form("hDenomPurityBC_V0M"),"; #V0flags in PF", 35, 0, 35);
    hDenomRejecEffBC_V0M = new TH1F(Form("hDenomRejecEffBC_V0M"),"; #V0flags in PF", 35, 0, 35);
    hNumRejecEffBC_V0M = new TH1F(Form("hNumRejecEffBC_V0M"),"; #V0flags in PF", 35, 0, 35);

    hSPDNumBC_V0M = new TH1F(Form("hSPDNumBC_V0M"),"; Spd tracklet", 200, 0, 200);
    hSPDDenomBC_V0M = new TH1F(Form("hSPDDenomBC_V0M"),"; Spd tracklet", 200, 0, 200);

    hNumTrkVsClsSPID_V0M = new TH2F(Form("hNumTrkVsClsSPID_V0M"),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
    hNumTrkVsClsSPID_V0M->GetXaxis()->SetTitle("Tracklet");
    hNumTrkVsClsSPID_V0M->GetYaxis()->SetTitle("Cluster (fspdC1)");
    hDenomTrkVsClsSPID_V0M= new TH2F(Form("hDenomTrkVsClsSPID_V0M"),"; Spd : !BGid",140,0,140,500,0,500);
    hDenomTrkVsClsSPID_V0M->GetXaxis()->SetTitle("Tracklet");
    hDenomTrkVsClsSPID_V0M->GetYaxis()->SetTitle("Cluster (fspdC1)");

    hNumV0_V0M = new TH2F(Form("hNumV0_V0M"),"; V0 : !BGid & GoodEvent",600,-300,300,2000,-1000,1000);
    hNumV0_V0M->GetXaxis()->SetTitle("V0A-V0C");
    hNumV0_V0M->GetYaxis()->SetTitle("V0A+V0C");
    hDenomV0_V0M= new TH2F(Form("hDenomV0_V0M"),"; V0 : !BGid",600,-300,300,2000,-1000,1000);
    hDenomV0_V0M->GetXaxis()->SetTitle("V0A-V0C");
    hDenomV0_V0M->GetYaxis()->SetTitle("V0A+V0C");

    fList2->Add(hNumEffPurityBC_V0M);
    fList2->Add(hDenomEffBC_V0M);
    fList2->Add(hDenomPurityBC_V0M);
    fList2->Add(hDenomRejecEffBC_V0M);
    fList2->Add(hNumRejecEffBC_V0M);
    fList2->Add(hSPDNumBC_V0M);
    fList2->Add(hSPDDenomBC_V0M);
    fList2->Add(hNumTrkVsClsSPID_V0M);
    fList2->Add(hDenomTrkVsClsSPID_V0M);
    fList2->Add(hNumV0_V0M);
    fList2->Add(hDenomV0_V0M);

    //__________SH2__________
    TH1F *hNumEffPurityBC_SH2;
    TH1F *hDenomEffBC_SH2;
    TH1F *hDenomPurityBC_SH2;
    TH1F *hDenomRejecEffBC_SH2;
    TH1F *hNumRejecEffBC_SH2;
    TH1F *hSPDNumBC_SH2;
    TH1F *hSPDDenomBC_SH2;
    TH2F *hNumTrkVsClsSPID_SH2;
    TH2F *hDenomTrkVsClsSPID_SH2;
    TH2F *hNumV0_SH2;
    TH2F *hDenomV0_SH2;


    hNumEffPurityBC_SH2 = new TH1F(Form("hNumEffPurityBC_SH2"),"; #V0flags in PF", 35, 0, 35);
    hDenomEffBC_SH2 = new TH1F(Form("hDenomEffBC_SH2"),"; #V0flags in PF", 35, 0, 35);
    hDenomPurityBC_SH2 = new TH1F(Form("hDenomPurityBC_SH2"),"; #V0flags in PF", 35, 0, 35);
    hDenomRejecEffBC_SH2 = new TH1F(Form("hDenomRejecEffBC_SH2"),"; #V0flags in PF", 35, 0, 35);
    hNumRejecEffBC_SH2 = new TH1F(Form("hNumRejecEffBC_SH2"),"; #V0flags in PF", 35, 0, 35);

    hSPDNumBC_SH2 = new TH1F(Form("hSPDNumBC_SH2"),"; Spd tracklet", 200, 0, 200);
    hSPDDenomBC_SH2 = new TH1F(Form("hSPDDenomBC_SH2"),"; Spd tracklet", 200, 0, 200);

    hNumTrkVsClsSPID_SH2 = new TH2F(Form("hNumTrkVsClsSPID_SH2"),"; Spd : !BGid & GoodEvent",140,0,140,500,0,500);
    hNumTrkVsClsSPID_SH2->GetXaxis()->SetTitle("Tracklet");
    hNumTrkVsClsSPID_SH2->GetYaxis()->SetTitle("Cluster (fspdC1)");
    hDenomTrkVsClsSPID_SH2= new TH2F(Form("hDenomTrkVsClsSPID_SH2"),"; Spd : !BGid",140,0,140,500,0,500);
    hDenomTrkVsClsSPID_SH2->GetXaxis()->SetTitle("Tracklet");
    hDenomTrkVsClsSPID_SH2->GetYaxis()->SetTitle("Cluster (fspdC1)");

    hNumV0_SH2 = new TH2F(Form("hNumV0_SH2"),"; V0 : !BGid & GoodEvent",600,-300,300,2000,-1000,1000);
    hNumV0_SH2->GetXaxis()->SetTitle("V0A-V0C");
    hNumV0_SH2->GetYaxis()->SetTitle("V0A+V0C");
    hDenomV0_SH2= new TH2F(Form("hDenomV0_SH2"),"; V0 : !BGid",600,-300,300,2000,-1000,1000);
    hDenomV0_SH2->GetXaxis()->SetTitle("V0A-V0C");
    hDenomV0_SH2->GetYaxis()->SetTitle("V0A+V0C");

    fList2->Add(hNumEffPurityBC_SH2);
    fList2->Add(hDenomEffBC_SH2);
    fList2->Add(hDenomPurityBC_SH2);
    fList2->Add(hDenomRejecEffBC_SH2);
    fList2->Add(hNumRejecEffBC_SH2);
    fList2->Add(hSPDNumBC_SH2);
    fList2->Add(hSPDDenomBC_SH2);
    fList2->Add(hNumTrkVsClsSPID_SH2);
    fList2->Add(hDenomTrkVsClsSPID_SH2);
    fList2->Add(hNumV0_SH2);
    fList2->Add(hDenomV0_SH2);

    //_______________________


    //_______________________________________
    
    TH1F *runNumber_hist;
    runNumber_hist = new TH1F("runNumber_hist","runNum", 1, 0, 1);
    fList->Add(runNumber_hist);
    TH1F *runNumber_hist_V0M;
    runNumber_hist_V0M = new TH1F("runNumber_hist_V0M","runNum", 1, 0, 1);
    fList2->Add(runNumber_hist_V0M);
    TH1F *runNumber_hist_SH2;
    runNumber_hist_SH2 = new TH1F("runNumber_hist_SH2","runNum", 1, 0, 1);
    fList2->Add(runNumber_hist_SH2);
    //______________________________
    
    TH2F *hTotalTrkVsClsSPID = new TH2F("hTotalTrkVsClsSPID","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTotalTrkVsClsSPID);

    TH2F *hTotalTrkVsClsSPID_PF2 = new TH2F("hTotalTrkVsClsSPID_PF2","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_PF2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_PF2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTotalTrkVsClsSPID_PF2);

    TH2F *hTotalTrkVsClsSPID_PF10 = new TH2F("hTotalTrkVsClsSPID_PF10","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_PF10->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_PF10->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTotalTrkVsClsSPID_PF10);

    //______________________________
    TH2F *hTotalTrkVsClsSPID_V0M = new TH2F("hTotalTrkVsClsSPID_V0M","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_V0M->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_V0M->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_V0M); 

    TH2F *hTotalTrkVsClsSPID_V0M_PF2 = new TH2F("hTotalTrkVsClsSPID_V0M_PF2","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_V0M_PF2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_V0M_PF2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_V0M_PF2); 

        TH2F *hTotalTrkVsClsSPID_V0M_PF10 = new TH2F("hTotalTrkVsClsSPID_V0M_PF10","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_V0M_PF10->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_V0M_PF10->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_V0M_PF10); 
    //______________________________
    TH2F *hTotalTrkVsClsSPID_SH2 = new TH2F("hTotalTrkVsClsSPID_SH2","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_SH2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_SH2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_SH2); 

        TH2F *hTotalTrkVsClsSPID_SH2_PF2 = new TH2F("hTotalTrkVsClsSPID_SH2_PF2","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_SH2_PF2->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_SH2_PF2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_SH2_PF2); 

        TH2F *hTotalTrkVsClsSPID_SH2_PF10 = new TH2F("hTotalTrkVsClsSPID_SH2_PF10","; Spd : total",140,0,140,500,0,500);
    hTotalTrkVsClsSPID_SH2_PF10->GetXaxis()->SetTitle("Tracklet");
    hTotalTrkVsClsSPID_SH2_PF10->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTotalTrkVsClsSPID_SH2_PF10); 
    //______________________________

    TH2F *hTotalV0 = new TH2F("hTotalV0","; V0 : total",600,-300,300,2000,-1000,1000);
    hTotalV0->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0->GetYaxis()->SetTitle("V0A+V0C");
    fList->Add(hTotalV0);
    
    TH2F *hTotalV0_V0M = new TH2F("hTotalV0_V0M","; V0 : total",600,-300,300,2000,-1000,1000);
    hTotalV0_V0M->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0_V0M->GetYaxis()->SetTitle("V0A+V0C");
    fList2->Add(hTotalV0_V0M); 
    
    TH2F *hTotalV0_SH2 = new TH2F("hTotalV0_SH2","; V0 : total",600,-300,300,2000,-1000,1000);
    hTotalV0_SH2->GetXaxis()->SetTitle("V0A-V0C");
    hTotalV0_SH2->GetYaxis()->SetTitle("V0A+V0C");
    fList2->Add(hTotalV0_SH2); 
    //______________________________
    TH2F *hTotalAD = new TH2F("hTotalAD","; AD : total",400,-200,200,2000,-1000,1000);
    hTotalAD->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD->GetYaxis()->SetTitle("ADA+ADC");
    fList->Add(hTotalAD);
    
    TH2F *hTotalAD_V0M = new TH2F("hTotalAD_V0M","; AD : total",400,-200,200,2000,-1000,1000);
    hTotalAD_V0M->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD_V0M->GetYaxis()->SetTitle("ADA+ADC");
    fList2->Add(hTotalAD_V0M); 
    
    TH2F *hTotalAD_SH2 = new TH2F("hTotalAD_SH2","; AD : total",400,-200,200,2000,-1000,1000);
    hTotalAD_SH2->GetXaxis()->SetTitle("ADA-ADC");
    hTotalAD_SH2->GetYaxis()->SetTitle("ADA+ADC");
    fList2->Add(hTotalAD_SH2); 

    //histogram for event list(blim)
    TH1F *hNumEvents  = new TH1F("hNumEvents","total event",10,0,10);
    fList->Add(hNumEvents);
    
    //_________histogram for modified cut____________
    
    TH2F *hTrkVsClsSPIDSlopeM = new TH2F("hTrkVsClsSPIDSlopeM","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlopeM);
    
    TH2F *hTrkVsClsSPIDSlopeM_V0M = new TH2F("hTrkVsClsSPIDSlopeM_V0M","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM_V0M->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_V0M->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTrkVsClsSPIDSlopeM_V0M); 

    TH2F *hTrkVsClsSPIDSlopeM_SH2 = new TH2F("hTrkVsClsSPIDSlopeM_SH2","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM_SH2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_SH2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTrkVsClsSPIDSlopeM_SH2); 

    
    TH2F *hTrkVsClsSPIDSlopeM2 = new TH2F("hTrkVsClsSPIDSlopeM2","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList->Add(hTrkVsClsSPIDSlopeM2);
    
    TH2F *hTrkVsClsSPIDSlopeM_V0M2 = new TH2F("hTrkVsClsSPIDSlopeM_V0M2","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM_V0M2->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_V0M2->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTrkVsClsSPIDSlopeM_V0M2); 
    
    TH2F *hTrkVsClsSPIDSlopeM_SH22 = new TH2F("hTrkVsClsSPIDSlopeM_SH22","; Spd : total",140,0,140,500,0,500);
    hTrkVsClsSPIDSlopeM_SH22->GetXaxis()->SetTitle("Tracklet");
    hTrkVsClsSPIDSlopeM_SH22->GetYaxis()->SetTitle("Cluster (fspdC1+fspdC2)");
    fList2->Add(hTrkVsClsSPIDSlopeM_SH22); 

    PostData(1, fList);
    PostData(2, fList2);
    
}

//________________________________________________________________________
void AliAnalysisBGMonitorQA::Exec(Option_t *)
{
    // Called for each event
    if (!fESD) {
        Printf("ERROR: fESD not available");
        return;
    }
    
    Int_t iEv= 0;
    iEv = fESD->GetEventNumberInFile();
    runNumber = fESD->GetRunNumber();

    ((TH1F*)fList->FindObject("runNumber_hist"))->SetBinContent(1,runNumber);
    ((TH1F*)fList2->FindObject("runNumber_hist_V0M"))->SetBinContent(1,runNumber);
    
    UInt_t timeGDC=fESD->GetTimeStamp();
    ftime=timeGDC;
    Int_t timeStampBX = fESD->GetBunchCrossNumber();
    fbx=timeStampBX;
    ntr = 10;
    nbunch = 21;
   // ofstream ftxt;
    
    static AliTriggerAnalysis * triggerAnalysis = new AliTriggerAnalysis();
    
    V0A = 0;
    V0C = 0;
    V0ABG = 0;
    V0CBG = 0;
    bgID = 0;
    
    // additional value initialize (blim)
    bgID2=0;
    
    VBA = 0;
    VBC = 0;
    VGA = 0;
    VTX = 0;
    fastORHW = 0;
    SPD1 = 0;
    SPD2 = 0;
    SPDHw1 = 0;
    SPDHw2 = 0;
    
    fastORHW = triggerAnalysis->EvaluateTrigger(fESD, AliTriggerAnalysis::kSPDGFO); // SPD number of chips from trigger bits (!)
    SPD1 = triggerAnalysis->SPDFiredChips(fESD,0,kFALSE,1);  //SPD Fired Chips in layer 1 (from cluster)
    SPD2 = triggerAnalysis->SPDFiredChips(fESD,0,kFALSE,2);  //SPD Fired Chips in layer 2 (from cluster)
    SPDHw1 = triggerAnalysis->SPDFiredChips(fESD,1,kFALSE,1);  //SPD Fired Chips in layer 1 (from hardware bit)
    SPDHw2 = triggerAnalysis->SPDFiredChips(fESD,1,kFALSE,2);  //SPD Fired Chips in layer 2 (from hardware bit)
    t0PileUp = triggerAnalysis->EvaluateTrigger(fESD, (AliTriggerAnalysis::Trigger) (AliTriggerAnalysis::kOfflineFlag | AliTriggerAnalysis::kT0Pileup)); //T0 pile-up
    
    AliVVZERO *vzero = fESD->GetVZEROData();
    V0A   = (vzero->GetV0ADecision()==AliVVZERO::kV0BB);
    V0ABG = (vzero->GetV0ADecision()==AliVVZERO::kV0BG);
    V0C   = (vzero->GetV0CDecision()==AliVVZERO::kV0BB);
    V0CBG = (vzero->GetV0CDecision()==AliVVZERO::kV0BG);
    
    AliAnalysisUtils *utils = new AliAnalysisUtils();
    //    bgID = utils->IsSPDClusterVsTrackletBG(fESD);
    
    // modified slope cut. the function is in below of this source(blim)
    bgID = IsItBGSPDClusterVsTracklet(fESD); // original modified function
    bgID2 = IsItBGSPDClusterVsTracklet2(fESD); // modified modified function
    
    spdPileUp = utils->IsPileUpSPD(fESD);
    spdPileUpOutOfBunch = utils->IsOutOfBunchPileUp(fESD);
    
    
    //CTP inputs
    VTX = fESD->GetHeader()->IsTriggerInputFired("0TVX");
    VGA = fESD->GetHeader()->IsTriggerInputFired("0VGA");
    VGC = fESD->GetHeader()->IsTriggerInputFired("0VGC");
    VBA = fESD->GetHeader()->IsTriggerInputFired("0VBA");
    VBC = fESD->GetHeader()->IsTriggerInputFired("0VBC");
    triMask = fESD->GetHeader()->GetTriggerMask();
    
    //--- vertex
    const AliESDVertex *vertSPD=fESD->GetPrimaryVertexSPD();
    if(vertSPD->GetNContributors()>0){
        fvertZ=vertSPD->GetZ();
        fvertX=vertSPD->GetX();
        fvertY=vertSPD->GetY();
    }
    else{
        fvertZ=-99999;
        fvertX=-99999;
        fvertY=-99999;
    }
    
    const AliESDVertex *vertTPC=fESD->GetPrimaryVertexTracks();
    if(vertTPC->GetNContributors()>0){
        fvertTPCZ=vertTPC->GetZ();
        fvertTPCX=vertTPC->GetX();
        fvertTPCY=vertTPC->GetY();
    }
    else{
        fvertTPCZ=-99999;
        fvertTPCX=-99999;
        fvertTPCY=-99999;
    }
    
    //--- SPD cluster and tracklets
    const AliMultiplicity* mult = fESD->GetMultiplicity();
    
    fSpdC1 = 0;
    fSpdC2 = 0;
    //for(Int_t ilayer = 0; ilayer < 2; ilayer++){
    //  fSpdC += mult->GetNumberOfITSClusters(ilayer);
    //}
    fSpdC1 = mult->GetNumberOfITSClusters(0);
    fSpdC2 = mult->GetNumberOfITSClusters(1);
    
    fSpdT = mult->GetNumberOfTracklets();
    
    //--- V0 data
    //AliESDVZERO* vzero = fESD->GetVZEROData();
    fv0a = vzero->GetV0ATime();  //V0A time
    fv0c = vzero->GetV0CTime();  //V0C time
    fMulta = vzero->GetMTotV0A();  //V0A multiplicity
    fMultc = vzero->GetMTotV0C();  //V0C multiplicity
    fTriCha = vzero->GetTriggerChargeA();  //Sum of the trigger (clock=10) charge on A side (Ring 0 excluded)
    fTriChc = vzero->GetTriggerChargeC();  //Sum of the trigger (clock=10) charge on A side
    
    
    
    //"online" V0 flags
    nV0A = 0;
    nV0ABG = 0;
    for (Int_t i = 32; i < 64; ++i) {
        if (vzero->GetBBFlag(i)) nV0A++;
        if (vzero->GetBGFlag(i)) nV0ABG++;
    }
    nV0C = 0;
    nV0CBG = 0;
    for (Int_t i = 0; i < 32; ++i) {
        if (vzero->GetBBFlag(i)) nV0C++;
        if (vzero->GetBGFlag(i)) nV0CBG++;
    }
    
     memset(BGFlagA, 0, sizeof(Float_t)*nbunch);
     memset(BBFlagA, 0, sizeof(Float_t)*nbunch);
     memset(BGFlagC, 0, sizeof(Float_t)*nbunch);
     memset(BBFlagC, 0, sizeof(Float_t)*nbunch);
    
    AliESDVZEROfriend *esdV0friend = fESDfriend->GetVZEROfriend();
    if(esdV0friend) {
        for(Int_t j = 0; j < 20; j++){
            //V0 --- infor
            for (Int_t i = 32; i < 64; ++i) {
                //BBFlagA[j] |= esdV0friend->GetBBFlag(i,j);
                //BGFlagA[j] |= esdV0friend->GetBGFlag(i,j);
                if(esdV0friend->GetBBFlag(i,j)) BBFlagA[j]++;
                if(esdV0friend->GetBGFlag(i,j)) BGFlagA[j]++;
            }
            for (Int_t i = 0; i < 32; ++i) {
                //BBFlagC[j] |= esdV0friend->GetBBFlag(i,j);
                //BGFlagC[j] |= esdV0friend->GetBGFlag(i,j);
                if(esdV0friend->GetBBFlag(i,j)) BBFlagC[j]++;
                if(esdV0friend->GetBGFlag(i,j)) BGFlagC[j]++;
            }
        }
    } else {
        Printf("No esdV0friend available");
        return;
    }
    
    ntracks = fESD->GetNumberOfTracks(); // number of tracks (no quality cuts)
    
    //--- Trigger classes --//
    memset(ftrigger, 0, sizeof(Float_t)*ntr);
    
    //Minimum Bias
    if(fESD->IsTriggerClassFired("CINT7-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT7-S-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-B-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT1-S-NOPF-ALLNOTRD") || fESD->IsTriggerClassFired("CINT7-A-NOPF-CENT") || fESD->IsTriggerClassFired("CINT7-B-NOPF-CENT") || fESD->IsTriggerClassFired("CINT7-C-NOPF-CENT") || fESD->IsTriggerClassFired("CINT7-E-NOPF-CENT")) ftrigger[0] = 1; // CINT7 trigger

    if(fESD->IsTriggerClassFired("CVHMV0M-A-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-B-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-C-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-E-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMV0M-B-SPD1-CENT") || fESD->IsTriggerClassFired("CVHMV0M-A-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-B-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-C-NOPF-CENTNOTRD") || fESD->IsTriggerClassFired("CVHMV0M-E-NOPF-CENTNOTRD")) ftrigger[2] = 1; // VOM trigger

    if(fESD->IsTriggerClassFired("CVHMSH2-A-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-B-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-C-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-E-NOPF-CENT") || fESD->IsTriggerClassFired("CVHMSH2-B-NOPF-ALL") || fESD->IsTriggerClassFired("CVHMSH2-B-NOPF-CENTNOTRD")) ftrigger[3] = 1; // SH2 trigger
    
    // count total event number (blim)
    if(ftrigger[0]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(1);
    if(ftrigger[2]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(2);
    if(ftrigger[3]==1)((TH1F*)fList->FindObject("hNumEvents"))->Fill(3);
    
    
    int nAfterBunch = 3;
    int nV0 = 3;
    int nFlag = 3;
    
    // Bool_t goodEvent;
    static Bool_t SelGoodEvent;
    
    
    if(ftrigger[0]) {  // trigger class for MB
        Printf("CINT7 triggred");        
        ((TH1F*)fList->FindObject("hTotalTrkVsClsSPID"))->Fill(fSpdT, fSpdC1+fSpdC2);
        ((TH1F*)fList->FindObject("hTotalV0"))->Fill(fv0a-fv0c, fv0a+fv0c);
        
        
        // Modified Cut result, added by blim
        if (!bgID) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlopeM"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID
        }
        if (!bgID2) {
            ((TH1F*)fList->FindObject("hTrkVsClsSPIDSlopeM2"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID2
        }
        
        for(Int_t ii=1; ii<33; ii++){
            
            //___________
            SelGoodEvent = BBFlagA[11]<ii  &  BBFlagA[12]<ii  &  BBFlagA[13]<ii  &  BBFlagA[14]<ii  &  BBFlagA[15]<ii  &  BBFlagA[16]<ii  & BBFlagA[17]<ii; //BB-A 11-17 
            SelGoodEvent &= BBFlagC[11]<ii  &  BBFlagC[12]<ii  &  BBFlagC[13]<ii  &  BBFlagC[14]<ii  &  BBFlagC[15]<ii  &  BBFlagC[16]<ii  &  BBFlagC[17]<ii; //BB-C 11-17
            SelGoodEvent &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  &  BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii; //BB-A 3-9
            SelGoodEvent &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  &  BBFlagC[7]<ii  &  BBFlagC[6]<ii  &  BBFlagC[5]<ii  &  BBFlagC[4]<ii  & BBFlagC[3]<ii; //BB-C 3-9
            //___________
            //printf(SelGoodEvent);
            if(SelGoodEvent) {
                ((TH1F*)fList->FindObject(Form("hDenomPurityBC")))->Fill(ii-1);
                if(ii == 2){
                            ((TH1F*)fList->FindObject("hTotalTrkVsClsSPID_PF2"))->Fill(fSpdT, fSpdC1+fSpdC2);
                }
                if(ii == 10){
                            ((TH1F*)fList->FindObject("hTotalTrkVsClsSPID_PF10"))->Fill(fSpdT, fSpdC1+fSpdC2);
                }
            }
            
            if(!bgID) {
                ((TH1F*)fList->FindObject(Form("hDenomEffBC")))->Fill(ii-1);
                ((TH1F*)fList->FindObject(Form("hSPDDenomBC")))->Fill(fSpdT);
                ((TH1F*)fList->FindObject(Form("hDenomTrkVsClsSPID")))->Fill(fSpdT, fSpdC1+fSpdC2);
                ((TH1F*)fList->FindObject(Form("hDenomV0")))->Fill(fv0a-fv0c, fv0a+fv0c);
                
                if(SelGoodEvent){
                    ((TH1F*)fList->FindObject(Form("hNumEffPurityBC")))->Fill(ii-1);
                    ((TH1F*)fList->FindObject(Form("hSPDNumBC")))->Fill(fSpdT);
                    ((TH1F*)fList->FindObject(Form("hNumTrkVsClsSPID")))->Fill(fSpdT, fSpdC1+fSpdC2);
                    ((TH1F*)fList->FindObject(Form("hNumV0")))->Fill(fv0a-fv0c, fv0a+fv0c);
                }
            }
            if(bgID){
                ((TH1F*)fList->FindObject(Form("hDenomRejecEffBC")))->Fill(ii-1);
                if(!SelGoodEvent){
                    ((TH1F*)fList->FindObject(Form("hNumRejecEffBC")))->Fill(ii-1);
                }
            }
        } // end of V0 flag loop
    } // end of events in trigger loop
    //-------------------------------------------------------V0M-------------------------------------------------------
    if(ftrigger[2]) {  // trigger class for HM // add new List for both result 2015.08.20. (blim)
        Printf("V0M triggred");
        ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_V0M"))->Fill(fSpdT, fSpdC1+fSpdC2);
        ((TH1F*)fList2->FindObject("hTotalV0_V0M"))->Fill(fv0a-fv0c, fv0a+fv0c);

        // Modified Cut result, added by blim
        if (!bgID) {
            ((TH1F*)fList2->FindObject("hTrkVsClsSPIDSlopeM_V0M"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID3->slope3
        }
        if (!bgID2) {
            ((TH1F*)fList2->FindObject("hTrkVsClsSPIDSlopeM_V0M2"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID2
        }
        
        for(Int_t ii=1; ii<33; ii++){
            
            
                       //___________
            SelGoodEvent = BBFlagA[11]<ii  &  BBFlagA[12]<ii  &  BBFlagA[13]<ii  &  BBFlagA[14]<ii  &  BBFlagA[15]<ii  &  BBFlagA[16]<ii  & BBFlagA[17]<ii; //BB-A 11-17 
            SelGoodEvent &= BBFlagC[11]<ii  &  BBFlagC[12]<ii  &  BBFlagC[13]<ii  &  BBFlagC[14]<ii  &  BBFlagC[15]<ii  &  BBFlagC[16]<ii  &  BBFlagC[17]<ii; //BB-C 11-17
            SelGoodEvent &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  &  BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii; //BB-A 3-9
            SelGoodEvent &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  &  BBFlagC[7]<ii  &  BBFlagC[6]<ii  &  BBFlagC[5]<ii  &  BBFlagC[4]<ii  & BBFlagC[3]<ii; //BB-C 3-9
            //___________

            //___________
            if(SelGoodEvent) {
                ((TH1F*)fList2->FindObject(Form("hDenomPurityBC_V0M")))->Fill(ii-1);
                if(ii == 2){
                            ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_V0M_PF2"))->Fill(fSpdT, fSpdC1+fSpdC2);
                }
                if(ii == 10){
                            ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_V0M_PF10"))->Fill(fSpdT, fSpdC1+fSpdC2);
                }
            }
            
            if(!bgID) {
                ((TH1F*)fList2->FindObject(Form("hDenomEffBC_V0M")))->Fill(ii-1);
                ((TH1F*)fList2->FindObject(Form("hSPDDenomBC_V0M")))->Fill(fSpdT);
                ((TH1F*)fList2->FindObject(Form("hDenomTrkVsClsSPID_V0M")))->Fill(fSpdT, fSpdC1+fSpdC2);
                ((TH1F*)fList2->FindObject(Form("hDenomV0_V0M")))->Fill(fv0a-fv0c, fv0a+fv0c);
                
                if(SelGoodEvent){
                    ((TH1F*)fList2->FindObject(Form("hNumEffPurityBC_V0M")))->Fill(ii-1);
                    ((TH1F*)fList2->FindObject(Form("hSPDNumBC_V0M")))->Fill(fSpdT);
                    ((TH1F*)fList2->FindObject(Form("hNumTrkVsClsSPID_V0M")))->Fill(fSpdT, fSpdC1+fSpdC2);
                    ((TH1F*)fList2->FindObject(Form("hNumV0_V0M")))->Fill(fv0a-fv0c, fv0a+fv0c);

                }
            }
            if(bgID){
                ((TH1F*)fList2->FindObject(Form("hDenomRejecEffBC_V0M")))->Fill(ii-1);
                if(!SelGoodEvent){
                    ((TH1F*)fList2->FindObject(Form("hNumRejecEffBC_V0M")))->Fill(ii-1);
                }

            }
        } // end of V0 flag loop  
    } // end of events in trigger loop
    //-------------------------------------------------------SH2-------------------------------------------------------
        if(ftrigger[3]) {  // trigger class for HM // add new List for both result 2015.08.20. (blim)
        Printf("SH2 triggred");
        ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_SH2"))->Fill(fSpdT, fSpdC1+fSpdC2);
        ((TH1F*)fList2->FindObject("hTotalV0_SH2"))->Fill(fv0a-fv0c, fv0a+fv0c);

        // Modified Cut result, added by blim
        if (!bgID) {
            ((TH1F*)fList2->FindObject("hTrkVsClsSPIDSlopeM_SH2"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID3->slope3
        }
        if (!bgID2) {
            ((TH1F*)fList2->FindObject("hTrkVsClsSPIDSlopeM_SH22"))->Fill(fSpdT, fSpdC1+fSpdC2); // bgID2
        }
        
        for(Int_t ii=1; ii<33; ii++){
            
            
                       //___________
            SelGoodEvent = BBFlagA[11]<ii  &  BBFlagA[12]<ii  &  BBFlagA[13]<ii  &  BBFlagA[14]<ii  &  BBFlagA[15]<ii  &  BBFlagA[16]<ii  & BBFlagA[17]<ii; //BB-A 11-17 
            SelGoodEvent &= BBFlagC[11]<ii  &  BBFlagC[12]<ii  &  BBFlagC[13]<ii  &  BBFlagC[14]<ii  &  BBFlagC[15]<ii  &  BBFlagC[16]<ii  &  BBFlagC[17]<ii; //BB-C 11-17
            SelGoodEvent &= BBFlagA[9]<ii  &  BBFlagA[8]<ii  &  BBFlagA[7]<ii  & BBFlagA[6]<ii  & BBFlagA[5]<ii  & BBFlagA[4]<ii  & BBFlagA[3]<ii; //BB-A 3-9
            SelGoodEvent &= BBFlagC[9]<ii  &  BBFlagC[8]<ii  &  BBFlagC[7]<ii  &  BBFlagC[6]<ii  &  BBFlagC[5]<ii  &  BBFlagC[4]<ii  & BBFlagC[3]<ii; //BB-C 3-9
            //___________

            //___________
            if(SelGoodEvent) {
                ((TH1F*)fList2->FindObject(Form("hDenomPurityBC_SH2")))->Fill(ii-1);
                if(ii == 2){
                            ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_SH2_PF2"))->Fill(fSpdT, fSpdC1+fSpdC2);
                }
                if(ii == 10){
                            ((TH1F*)fList2->FindObject("hTotalTrkVsClsSPID_SH2_PF10"))->Fill(fSpdT, fSpdC1+fSpdC2);
                }
            }
            
            if(!bgID) {
                ((TH1F*)fList2->FindObject(Form("hDenomEffBC_SH2")))->Fill(ii-1);
                ((TH1F*)fList2->FindObject(Form("hSPDDenomBC_SH2")))->Fill(fSpdT);
                ((TH1F*)fList2->FindObject(Form("hDenomTrkVsClsSPID_SH2")))->Fill(fSpdT, fSpdC1+fSpdC2);
                ((TH1F*)fList2->FindObject(Form("hDenomV0_SH2")))->Fill(fv0a-fv0c, fv0a+fv0c);
                
                if(SelGoodEvent){
                    ((TH1F*)fList2->FindObject(Form("hNumEffPurityBC_SH2")))->Fill(ii-1);
                    ((TH1F*)fList2->FindObject(Form("hSPDNumBC_SH2")))->Fill(fSpdT);
                    ((TH1F*)fList2->FindObject(Form("hNumTrkVsClsSPID_SH2")))->Fill(fSpdT, fSpdC1+fSpdC2);
                    ((TH1F*)fList2->FindObject(Form("hNumV0_SH2")))->Fill(fv0a-fv0c, fv0a+fv0c);

                }
            }
            if(bgID){
                ((TH1F*)fList2->FindObject(Form("hDenomRejecEffBC_SH2")))->Fill(ii-1);
                if(!SelGoodEvent){
                    ((TH1F*)fList2->FindObject(Form("hNumRejecEffBC_SH2")))->Fill(ii-1);
                }

            }
        } // end of V0 flag loop  
    } // end of events in trigger loop

    PostData(1, fList);
    PostData(2, fList2);
    fTreeTrack2->Fill();
    PostData(0, fTreeTrack2);
}


//________________________________________________________________________
void AliAnalysisBGMonitorQA::Terminate(Option_t *)
{
    //   fList = dynamic_cast<TList*> (GetOutputData(1));
    //   if(!fList)    Printf("ERROR: fList is not available");
}



//______________________________________________________________________ Modified cut function(blim)
Bool_t IsItBGSPDClusterVsTracklet(AliVEvent *event)
{
    /*
     Int_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
     Int_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
     Int_t nTracklets      = event->GetMultiplicity()->GetNumberOfTracklets();
     if (nClustersLayer0 + nClustersLayer1 > 65. + nTracklets*slope) return kTRUE;
     return kFALSE;
     */
    Double_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
    Double_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
    Double_t trk = event->GetMultiplicity()->GetNumberOfTracklets();
    Double_t cls = nClustersLayer0 + nClustersLayer1;
    Bool_t spdBg = kFALSE;
    
    if (trk < 1.5) {
        spdBg = kTRUE;
    }
    else if (trk >= 1.5 && trk < 26 && cls > 20.0 + (378-20)/(26-1.5)*(trk-1.5)) {
        spdBg = kTRUE;
    }
    else if (trk >= 26 && trk < 40 && cls > 378.0 + (505-378)/(40-26.)*(trk-26.)) {
        spdBg = kTRUE;
    }
    else if (trk >= 40 && cls > 505.0 + (1770-505)/(300-40.)*(trk-40.)) {
        spdBg = kTRUE;
    }
    
    return spdBg;
}
//______________________________________________________________________ Modified cut function(blim)
Bool_t IsItBGSPDClusterVsTracklet2(AliVEvent *event)
{
    Double_t nClustersLayer0 = event->GetNumberOfITSClusters(0);
    Double_t nClustersLayer1 = event->GetNumberOfITSClusters(1);
    Double_t trk = event->GetMultiplicity()->GetNumberOfTracklets();
    Double_t cls = nClustersLayer0 + nClustersLayer1;
    Bool_t spdBg = kFALSE;
    
    if (trk < 5.944 && cls > 65 + 4.*trk) {
        spdBg = kTRUE;
    }
    else if (trk >= 5.944 && trk < 26 && cls > 20.0 + (378-20)/(26-1.5)*(trk-1.5)) {
        spdBg = kTRUE;
    }
    else if (trk >= 26 && trk < 40 && cls > 378.0 + (505-378)/(40-26.)*(trk-26.)) {
        spdBg = kTRUE;
    }
    else if (trk >= 40 && cls > 505.0 + (1770-505)/(300-40.)*(trk-40.)) {
        spdBg = kTRUE;
    }
    
    return spdBg;
}