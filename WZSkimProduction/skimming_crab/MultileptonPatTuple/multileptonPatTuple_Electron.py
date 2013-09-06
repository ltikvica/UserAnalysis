## Configuration based on PhysicsTools.PatAlgos.patTemplate_cfg

import FWCore.ParameterSet.Config as cms

## This cfg file sets up the process object, loads the geometry & detector
## conditions, and calls the standard PAT sequences.
from PhysicsTools.PatAlgos.patTemplate_cfg import process

## Load functions to tailor PAT for our needs
from UserAnalysis.WZAnalysis.wzPatTools import *

## The 'analysis mode' selected below will determine selections about modules
## and filters to include.  For generated samples, mode should be "MC"; for
## data, mode should be "EG" or "Mu" depending on the dataset.  The hltMuFilter
## will be used to ensure no overlap between the EG and Mu datasets.
mode = analysisMode("EG")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.MessageLogger.suppressWarning += ["patTrigger"]
process.MessageLogger.cerr.limit = cms.uint32(1)

## Options and Output Report
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

## GlobalTag
#process.GlobalTag.globaltag = 'GR_P_V14::All'
process.GlobalTag.globaltag = 'GR_R_311_V2::All'

## Source
process.source.fileNames = cms.untracked.vstring(
#    "/store/mc/Fall10/WZTo3LNu_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0012/2C675164-F7C9-DF11-9E22-001E0B62A9E0.root"
#'/store/mc/Spring11/QCDPt9_MuElPt15EtaFilter_7TeV-pythia6/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0008/34B4E7BB-B54F-E011-9CAB-003048678D9A.root'

'/store/data/Run2011A/SingleElectron/AOD/PromptReco-v1/000/160/835/080DD672-CB53-E011-A2BB-000423D9997E.root'

# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/1A62CAA7-8153-E011-B37E-0030487CAF5E.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/2058A542-8353-E011-9612-001D09F28D54.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/38E2DC8A-7F53-E011-99EF-003048F11C58.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/5A89378A-7F53-E011-9E7D-003048F11942.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/627BCB80-7853-E011-8EA3-003048F118C4.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/92506D2B-7C53-E011-9BF3-003048F024C2.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/92CA84F4-8053-E011-A631-003048F118C4.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/9E2AE9A5-8153-E011-8A56-0030487C6A66.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/AE9D6D42-8353-E011-A261-001D09F24934.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/DAB13D30-8553-E011-9DF8-0019DB2F3F9A.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/EACA29D4-CB53-E011-BC56-003048F11CF0.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/F493F6F2-8053-E011-98E3-001D09F295A1.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/F4C79623-7C53-E011-B9E7-001D09F34488.root',
# '/store/data/Run2011A/SingleMu/AOD/PromptReco-v1/000/160/835/F6EB44FD-7D53-E011-AF67-001D09F2AD4D.root',
)

## Output
process.out.fileName = cms.untracked.string('TrileptonPatTuple.root')
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p3Lep'))
process.out.outputCommands = cms.untracked.vstring('drop *')
process.out.dropMetaData = cms.untracked.string('DROPPED')

process.out2Lep = process.out.clone()
process.out2Lep.fileName = cms.untracked.string('DileptonPatTuple.root')
process.out2Lep.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p2Lep'))
process.out2Lep.outputCommands = cms.untracked.vstring('drop *')
process.out2Lep.dropMetaData = cms.untracked.string('DROPPED')

## Load PAT
optimizePat(process, cleaning=False, trigger=True, mode=mode)

## Counters
process.nEventsTotal    = cms.EDProducer("EventCountProducer")
process.nEventsHLT      = cms.EDProducer("EventCountProducer")
process.nEventsFiltered = cms.EDProducer("EventCountProducer")
process.nEventsPat      = cms.EDProducer("EventCountProducer")

## Cross-Section
process.crossSection = cms.EDProducer("DoubleProducer", value=cms.double(0.))

## Reduce PAT file size
process.countPatLeptons.minNumber = 3
process.selectedPatJets.cut = "et > 20. & abs(eta) < 2.5"
process.selectedPatMuons.cut = "pt > 10 & abs(eta) < 2.4 & isGlobalMuon"
process.selectedPatElectrons.cut = "et > 10 & abs(eta) < 2.5"
process.out.outputCommands += wzEventContent
process.out2Lep.outputCommands += wzEventContent

## Prepare output
process.p = cms.Sequence(
    process.nEventsTotal *
    process.hltMuFilter *
    process.nEventsHLT *
    process.dataFilters *
    process.nEventsFiltered *
    process.userPatSequence *
    process.nEventsPat *
    process.wzPreselectionProducer *
    process.crossSection
)
configureFilters(process, mode)

process.outpath = cms.EndPath(process.out)
process.p3Lep = cms.Path(process.p)

process.outpath2Lep = cms.EndPath(process.out2Lep)
process.p2Lep = cms.Path(process.p)

process.countPatLeptons2Lep = process.countPatLeptons.clone()
process.countPatLeptons2Lep.minNumber = 2
process.countPatLeptons2Lep.maxNumber = 2
process.p2Lep.replace(process.countPatLeptons, process.countPatLeptons2Lep)

