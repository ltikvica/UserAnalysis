import FWCore.ParameterSet.Config as cms

wzTreeMaker = cms.EDAnalyzer('WZTreeMaker',
                              
    outputFileName     = cms.untracked.string("outputTree.root"),
    hltProcessName     = cms.untracked.string("*"              ),
    debug              = cms.untracked.bool  (False            ),

    #primaryVertexTag = cms.untracked.InputTag("offlinePrimaryVertices" ),
    primaryVertexTag = cms.untracked.InputTag("goodVertices" ),
     genParticlesTag = cms.untracked.InputTag("prunedGenParticles"           ),
# genParticlesTag = cms.untracked.InputTag("genParticles"           ),
          genInfoTag = cms.untracked.InputTag("generator"              ),
     patTrigEventTag = cms.untracked.InputTag("patTriggerEvent"        ),
       hltResultsTag = cms.untracked.InputTag("TriggerResults"         ),
       hltSummaryTag = cms.untracked.InputTag("hltTriggerSummaryAOD"   ),
           pileupTag = cms.untracked.InputTag("addPileupInfo"          ),
              zllTag = cms.untracked.InputTag("zToLL"                  ),
         electronTag = cms.untracked.InputTag("selectedPatElectrons"       ),
#electronTag = cms.untracked.InputTag("userPatElectrons"       ),
             muonTag = cms.untracked.InputTag("selectedPatMuons"           ),
#muonTag = cms.untracked.InputTag("userPatMuons"           ),
              jetTag = cms.untracked.InputTag("selectedPatJets"        ),
            trackTag = cms.untracked.InputTag("generalTracks"          ),

    metPSets  = cms.untracked.VPSet(
        cms.PSet(name = cms.string("met"),
                 tag  = cms.InputTag("patMETsPF")),
#             tag  = cms.InputTag("patMETs")),
        cms.PSet(name = cms.string("tcMet"),
                 tag  = cms.InputTag("patMETsPF")),
#                         tag  = cms.InputTag("patMETsTC")),
        cms.PSet(name = cms.string("pfMet"),
#                 tag  = cms.InputTag("patMETs")),
                 tag  = cms.InputTag("patMETsPF")),
        cms.PSet(name = cms.string("t1Met"),
                 tag  = cms.InputTag("patMETsPF")),
# tag  = cms.InputTag("patMETs")),
#tag  = cms.InputTag("patMETsPFType1")),
    ),
   
    metT1Label = cms.untracked.InputTag("pfType1CorrectedMet"),

    numEventsNames = cms.untracked.vstring(),
    weightProducerName = cms.untracked.string("userWeightLB"),
    contentToDrop    = cms.vstring(),
    contentToKeep    = cms.vstring(),

    genIdsToStore    = cms.vint32(11, 12, 13, 14, 15, 16, 22, 23, 24),

    triggersToStore  = cms.vstring(
         "HLT_PhysicsDeclared",
         "HLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v2"
    ),
    hltMuPaths = cms.vstring(),
    hltElePaths = cms.vstring(),
    hltUpdatedForRun = cms.int32(0)

#    MCPUFile = cms.string('${CMSSW_BASE}/src/UserAnalysis/WZAnalysis/test/MCSummer11S4PUDist.root'),
#    DataPUFile = cms.string('${CMSSW_BASE}/src/UserAnalysis/WZAnalysis/test/DataPUDist.root'), #3BX method
#    DataPUFileTruth = cms.string('${CMSSW_BASE}/src/UserAnalysis/WZAnalysis/test/DataPUDistTruth.root'),#Fall11 Truth method
#    DataPUFile3D = cms.string('${CMSSW_BASE}/src/UserAnalysis/WZAnalysis/test/DataPUDistTruth_v2_finebin.root'),#3D method
#    MCPUHist = cms.string('PoissonOneXDist'),  #Summer11 3BX method
#    MCPUHist3D = cms.string('probdistFlat10'), #Summer11 3D method
#    MCPUHist3D = cms.string('Fall2011'),       #Fall11 3D method
#    MCPUHistTruth = cms.string('Fall2011'),    #Fall11 Truth method
#    DataPUHist = cms.string('pileup'),
)

