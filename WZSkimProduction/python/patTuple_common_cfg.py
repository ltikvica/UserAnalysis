from PhysicsTools.PatAlgos.patTemplate_cfg import *
from PhysicsTools.PatAlgos.tools.coreTools import *
from PhysicsTools.PatAlgos.tools.trigTools import switchOnTrigger
from PhysicsTools.PatAlgos.tools.metTools import addPfMET
from PhysicsTools.PatAlgos.tools.pfTools import addPFCandidates
from PhysicsTools.PatAlgos.tools.pfTools import *

def common_config(process, reportEveryNum=100, maxEvents=-1, runOnData=False) :
    usePFIso( process )
    if runOnData:
        removeMCMatching(process, ['All'])

    process.patElectrons.pfElectronSource = 'particleFlow'
    addPfMET(process, 'PF')

    switchOnTrigger(process)
    process.patTrigger.addL1Algos = cms.bool(True)
                
    # this is needed so we can correct the pfMET by adjusting the e/mu-pt
    # when switching to one of the dedicated Heep/TeV muon reconstructors
    addPFCandidates(process, 'particleFlow')
    process.selectedPatPFParticles.cut = "abs(pdgId())==11 || abs(pdgId())==13"

    process.load("FWCore.MessageLogger.MessageLogger_cfi")
    process.MessageLogger.cerr.FwkReport.reportEvery = reportEveryNum

    process.maxEvents.input = maxEvents    ##  (e.g. -1 to run on all events)
#    process.GlobalTag.globaltag = cms.string('GR_R_52_V8::All')
    process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')
    #                                         ##
    process.out.outputCommands = [
    # GEN
        'keep *_prunedGenParticles_*_*',
        'keep GenEventInfoProduct_*_*_*',
        'keep GenRunInfoProduct_*_*_*',
    # TRIGGER
        'keep edmTriggerResults_TriggerResults*_*_*',
        'keep *_hltTriggerSummaryAOD_*_*',
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*',
        'keep *_userPat*_*_*',
        'keep *_patTrigger_*_*',
        'keep *_patTriggerEvent_*_*',
    # PILEUP
        'keep *_addPileupInfo_*_*',     
    # PF CANDS
        'keep *_selectedPatPFParticles*_*_*'
        ]
    
##  (to suppress the long output at the end of the job)    
    process.options.wantSummary = True        


def common_filters(process) :
    ##Common Filters ############
    #CSC
    process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')
    process.p *= process.CSCTightHaloFilter
    
    #HBHE
    process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')
    process.p *= process.HBHENoiseFilter

    #HCAL Laser
    process.load("RecoMET.METFilters.hcalLaserEventFilter_cfi")
    process.p *= process.hcalLaserEventFilter

    #ECAL deal cell TP
    process.load('RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi')
    ## For AOD and RECO recommendation to use recovered rechits
    process.EcalDeadCellTriggerPrimitiveFilter.tpDigiCollection = cms.InputTag("ecalTPSkimNA")

    # The section below is for the filter on Boundary Energy. Available in AOD in CMSSW>44x
    process.load('RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi')
    process.EcalDeadCellBoundaryEnergyFilter.taggingMode = cms.bool(False)
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEB=cms.untracked.double(10)
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyDeadCellsEE=cms.untracked.double(10)
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEB=cms.untracked.double(100)
    process.EcalDeadCellBoundaryEnergyFilter.cutBoundEnergyGapEE=cms.untracked.double(100)
    process.EcalDeadCellBoundaryEnergyFilter.enableGap=cms.untracked.bool(False)
    process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEB = cms.vint32(12,14)
    process.EcalDeadCellBoundaryEnergyFilter.limitDeadCellToChannelStatusEE = cms.vint32(12,14)
    # End of Boundary Energy filter configuration

    # The line below is the default recommendation
    process.p *= process.EcalDeadCellTriggerPrimitiveFilter

    #Tracking Failure
    process.goodVertices = cms.EDFilter(
          "VertexSelector",
            filter = cms.bool(False),
            src = cms.InputTag("offlinePrimaryVertices"),
            cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.rho < 2")
          )

    process.load('RecoMET.METFilters.trackingFailureFilter_cfi')
    process.out.outputCommands.append('keep *_goodVertices_*_*')

    process.p *= process.goodVertices*process.trackingFailureFilter

    #Bad EE SC
    process.load('RecoMET.METFilters.eeBadScFilter_cfi')

    process.p *= process.eeBadScFilter

    #ECAL Laser

    #Tracking POG
