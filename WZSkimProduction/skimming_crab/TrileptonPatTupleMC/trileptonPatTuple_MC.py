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
mode = analysisMode("MC")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 50
process.MessageLogger.suppressWarning += ["patTrigger"]
process.MessageLogger.cerr.limit = cms.uint32(1)

## Options and Output Report
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

## Maximal Number of Events
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

## Source
process.source.fileNames = cms.untracked.vstring(
#    "/store/mc/Fall10/WZTo3LNu_TuneZ2_7TeV-pythia6/GEN-SIM-RECO/START38_V12-v1/0012/2C675164-F7C9-DF11-9E22-001E0B62A9E0.root"
'/store/mc/Spring11/QCDPt9_MuElPt15EtaFilter_7TeV-pythia6/GEN-SIM-RECO/PU_S1_START311_V1G1-v1/0008/34B4E7BB-B54F-E011-9CAB-003048678D9A.root'
)

## Output
process.out.fileName = cms.untracked.string('patTuple.root')
process.out.SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('p'))
process.out.outputCommands = cms.untracked.vstring('drop *')
process.out.dropMetaData = cms.untracked.string('DROPPED')

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

## Prepare output
process.outpath = cms.EndPath(process.out)
process.p = cms.Path(
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
