#ifndef AliAnalysisBGMonitorQA_h
#define AliAnalysisBGMonitorQA_h

#include "AliAnalysisTaskSE.h"


class AliESDEvent;
class AliESDfriend;
class AliAnalysisCuts;
class TH1D;
class TH1F;
class TH2F;
class TH2D;

class AliAnalysisBGMonitorQA : public AliAnalysisTaskSE {
 public:
  AliAnalysisBGMonitorQA(const char *name = "AliAnalysisBGMonitorQA");
  virtual ~AliAnalysisBGMonitorQA() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

    
  virtual void SelectGoodEventWithV0Variation(Int_t bunchrange, Int_t v0variation ,Int_t flagvariation, Int_t ii); // add function to select good event with 3 types of variation  2015.08.12. (blim)
  virtual void SelectADGoodEventWithV0Variation(Int_t bunchrange, Int_t v0variation ,Int_t flagvariation, Int_t ii); // add function to select good event with 3 types of variation in AD 2015.08.12. (blim)

 // virtual void   Terminate(Option_t *);
    
 private: 
  AliESDEvent *fESD;        //! ESD event
  AliESDfriend* fESDfriend; //! ESDfriend   
  TTree *fTreeTrack;        //! tree
  TTree *fTreeTrack2;        //! tree
  TList *fList;             //! list
  TList *fList2;             //! list for additional data 2015.08.20 (blim)
  Int_t fUseTree;


  Int_t runNumber = 0;
  Int_t ftrigger[kMaxUShort]; 
  Double_t fvertZ = 0;
  Double_t fvertX = 0;
  Double_t fvertY = 0;
  Double_t fvertTPCZ = 0;
  Double_t fvertTPCX = 0;
  Double_t fvertTPCY = 0;
  Double_t fvertZ2 = 0;
  Double_t fvertX2 = 0;
  Double_t fvertY2 = 0;
  Double_t fvertTPCZ2 = 0;
  Double_t fvertTPCX2 = 0;
  Double_t fvertTPCY2 = 0;
  Float_t fv0a = 0;
  Float_t fv0c = 0;
  Float_t fad0a = 0;
  Float_t fad0c = 0;
  Float_t fMulta = 0;
  Float_t fMultc = 0;
  Float_t fTriCha = 0;
  Float_t fTriChc = 0;
  Float_t fV0M = 0;
  
  Int_t fbx = 0;
  Int_t ftime = 0;
  Int_t fSpdC1 = 0;
  Int_t fSpdC2 = 0;
  Int_t fSpdT = 0;
  Int_t ntracks = 0;
  Int_t V0A = 0;
  Int_t V0C = 0;
  Int_t V0ABG = 0;
  Int_t V0CBG = 0;
  Int_t nV0A = 0;
  Int_t nV0C = 0;
  Int_t nV0ABG = 0;
  Int_t nV0CBG = 0;
  Int_t VBA = 0;
  Int_t VBC = 0;
  Int_t VGA = 0;
  Int_t VGC = 0;
  Int_t VTX = 0;
  Int_t bgID = 0;
    
  Int_t bgID2 = 0;
    
  Int_t t0PileUp = 0;
  Int_t spdPileUp = 0;
  Int_t spdPileUpOutOfBunch = 0;
  Int_t triMask = 0; 
  Int_t fastORHW = 0;  
  Int_t SPD1 = 0;
  Int_t SPD2 = 0;
  Int_t SPDHw1 = 0;
  Int_t SPDHw2 = 0; 
  Int_t BGFlagA[kMaxUShort];
  Int_t BGFlagC[kMaxUShort];
  Int_t BBFlagA[kMaxUShort];
  Int_t BBFlagC[kMaxUShort];
  Int_t ADBGFlagA[kMaxUShort];
  Int_t ADBGFlagC[kMaxUShort];
  Int_t ADBBFlagA[kMaxUShort];
  Int_t ADBBFlagC[kMaxUShort];
  
  UShort_t ntr = 0;
  UShort_t nbunch = 0;

  AliAnalysisBGMonitorQA(const AliAnalysisBGMonitorQA&); // not implemented
  AliAnalysisBGMonitorQA& operator=(const AliAnalysisBGMonitorQA&); // not implemented
  ClassDef(AliAnalysisBGMonitorQA, 2);// example of analysis
};

#endif
