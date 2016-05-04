AliAnalysisBGMonitorQA *AddTaskBGMonitorQA(Bool_t UseTree = kFALSE)
{
  //
  //This macro configures the task for the Beam Gas Monitoring QA
  //==============================================================================

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskMBVeto", "No analysis manager to connect to.");
    return NULL;
  }   


    // Input handlers
    AliESDInputHandler* esdH = new AliESDInputHandler();
    
    esdH->SetReadFriends(kTRUE);
    mgr->SetInputEventHandler(esdH);

    
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskMBVeto", "This task requires an input event handler");
    return NULL;
  }
    
  // Create and configure the task
  AliAnalysisBGMonitorQA *taskMBVeto = new AliAnalysisBGMonitorQA("taskMBVeto");
  mgr->AddTask(taskMBVeto);
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("cchain1",TChain::Class(), AliAnalysisManager::kInputContainer);
  mgr->ConnectInput(taskMBVeto, 0, mgr->GetCommonInputContainer());
  
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("cOutputH_CINT7", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName()));
  mgr->ConnectOutput(taskMBVeto, 1, coutput1);

    AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("cOutputH_V0MandSH2", TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName()));
    mgr->ConnectOutput(taskMBVeto, 2, coutput2);

    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer("cOutputT", TTree::Class(), AliAnalysisManager::kOutputContainer, Form("%s:BeamGasMon", mgr->GetCommonFileName()));
    mgr->ConnectOutput(taskMBVeto, 0, coutput3);
    
  return taskMBVeto;
  
}

