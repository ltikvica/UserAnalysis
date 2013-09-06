import FWCore.ParameterSet.Config as cms

process = cms.Process("RUN")

readFiles = cms.untracked.vstring()
process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            #Run: 1 Event: 80333
                            duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
                            fileNames = readFiles
                            )
readFiles.extend([
        #'file:cmgTuple.root'
        #'/store/cmst3/user/sbologne/Higgs/cmgTuples/ZJetsMAD.root'
        #'file:/cmsrm/pc24_2/emanuele/data/DYeeSummer11.root',
#         '/store/user/fladias/V350-Summer11_RSZZmmjj_m1750/patTuple_5_1_jLM.root'
#    '/store/user/fladias/V350-Summer11_DYJetsToLL_PtZ100/patTuple_27_1_e4Z.root'
    '/store/user/fladias/V390C-Fall11-DYJetsToLL_PtZ-100/patTuple_451_1_LYt.root'
         ])
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1))

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V24::All'

process.readAK7PF = cms.EDAnalyzer('JetCorrectorDBReader',
                                   payloadName    = cms.untracked.string('AK7PF'), # this is the communication to the databas
                                   printScreen    = cms.untracked.bool(False),
                                   createTextFile = cms.untracked.bool(True),
                                   globalTag      = cms.untracked.string('DY-GR_R_42_V24') # this is used ONLY for the name of the printed txt files. You can use any name that you like, but it is recommended  .
                                   )
process.p = cms.Path(process.readAK7PF)
