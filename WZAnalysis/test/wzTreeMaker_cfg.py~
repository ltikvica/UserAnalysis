import FWCore.ParameterSet.Config as cms

import os
import sys

process = cms.Process("ANALYSIS")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load("UserAnalysis.WZAnalysis.wzTreeMaker_cfi")


# include electron energy corrections here
#process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")

# dataset to correct
#process.calibratedPatElectrons.isMC = cms.bool(False)
#process.calibratedPatElectrons.inputDataset = cms.string("Moriond2013")
#process.calibratedPatElectrons.updateEnergyError = cms.bool(True)
#process.calibratedPatElectrons.correctionsType = cms.int32(1)
#process.calibratedPatElectrons.combinationType = cms.int32(3)
#process.calibratedPatElectrons.lumiRatio = cms.double(1.0)
#process.calibratedPatElectrons.verbose = cms.bool(True)
#process.calibratedPatElectrons.synchronization = cms.bool(True)
# end of  electron energy corrections here



from UserAnalysis.WZAnalysis.wzCommonParameters import *

## Parse arguments from the command-line
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.setDefault('inputFiles', 'file:patTuple.root') 
options.setDefault('outputFile', 'outputTree.root')
#options.setDefault('outputFile', 'outputTreeMC5.root')
#options.setDefault('outputFile', 'outputTreeData2.root')
#options.setDefault('outputFile', 'outputTreeMC.root')

options.register('debug', False, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "Enable debug output")
options.parseArguments()

process.wzTreeMaker.numEventsNames = ['nEventsTotal',
                                      'nEventsHLT',
                                      'nEventsFiltered',
                                      'nEventsPat',
                                      'nEventsZ']

# HLT paths to use for trigger matching, etc.
# Underscores are removed to match the labels of the matching modules
process.wzTreeMaker.hltMuPaths = singleMuonPaths + doubleMuonPaths
process.wzTreeMaker.hltElePaths = singleElectronPaths + doubleElectronPaths
process.wzTreeMaker.triggersToStore += muegPaths

# Manage text output
process.MessageLogger.cerr.FwkReport.reportEvery = 50
# Accessing prescales sometimes generates copious errors on HLTConfigData
process.MessageLogger.categories += ['HLTConfigData']
process.MessageLogger.cerr.HLTConfigData = cms.untracked.PSet(
    limit = cms.untracked.int32(0))

# Load GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.PyReleaseValidation.autoCond import autoCond
#process.GlobalTag.globaltag = cms.string( autoCond[ 'startup' ] )
#global tag for data:
#process.GlobalTag.globaltag = 'FT_R_53_V18::All'
#global tag for MC:
process.GlobalTag.globaltag = 'START53_V7G::All'

# Any branch whose name matches one of these perl-style regular expressions
# will not be generated in the ouput tree (use '.*' for wildcard)
process.wzTreeMaker.contentToDrop = []
process.wzTreeMaker.contentToKeep = []

#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

#using the 'fillIn' function to process all files in a specified directory
#from UserAnalysis.WZAnalysis.wzPatTools import *
#process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring(fillIn('/scratch/fantasia/forVuko/data/patTuple_1_1_DBX.root'))
#    )
#)

## Prepare output
process.wzTreeMaker.outputFileName = options.outputFile

process.p = cms.Path(process.wzTreeMaker)

