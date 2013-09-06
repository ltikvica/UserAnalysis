from UserAnalysis.WZSkimProduction.patTuple_common_cfg import *
from UserAnalysis.WZSkimProduction.patTuple_el_cfg import *
from UserAnalysis.WZSkimProduction.patTuple_mu_cfg import *
from UserAnalysis.WZSkimProduction.patTuple_met_cfg import *
from UserAnalysis.WZSkimProduction.patTuple_jetlep_cfg import addFastJet

def trilep_config(process, reportEveryNum=100, maxEvents=-1, runOnData=False) :
    common_config(process, reportEveryNum, maxEvents, runOnData)
    el_config(process)
    mu_config(process)
    met_config(process)
    addFastJet(process)

    # redefine selectedPatMuons (isGlobalMuon not included in std definition)
    process.selectedPatElectrons.cut = "pt > 10. & abs(eta) < 2.5"
    process.selectedPatMuons.cut = "pt > 10. & abs(eta) < 2.4 & (isGlobalMuon | isTrackerMuon)"

    process.selectedPatJets.cut = "pt > 30. & abs(eta) < 3.0"
    process.out.outputCommands.append('keep *_selectedPatJets_*_*')

    # keep all events with 3 leptons above 10 GeV
    process.countPatLeptons.electronSource = "selectedPatElectrons"
    process.countPatLeptons.muonSource     = "selectedPatMuons"
    process.countPatLeptons.minNumber = 3

    process.selectedPatPFParticles.cut = "abs(pdgId())==13 || abs(pdgId())==11"
    process.out.outputCommands.append('keep *_selectedPatMuons_*_*')
    process.out.outputCommands.append('keep *_selectedPatElectrons_*_*')

    ###### Z Producers
    process.dimu = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("selectedPatMuons@+ selectedPatMuons@-"),
                                  cut   = cms.string("60 < mass < 120"),
                                  )
    process.diel = cms.EDProducer("CandViewShallowCloneCombiner",
                                  decay = cms.string("selectedPatElectrons@+ selectedPatElectrons@-"),
                                  cut   = cms.string("60 < mass < 120"),
                                  )

    # merge muons and electrons into leptons
    process.dilep = cms.EDProducer("CandViewMerger",
                               src = cms.VInputTag("dimu", "diel")
                               )

    ## Di-muons
    process.zmass = cms.EDFilter("CandViewCountFilter",
                                 src = cms.InputTag('dilep'),
                                 minNumber = cms.uint32(1)
                                 )
    
    process.diLepSelSeq = cms.Sequence(process.dimu * process.diel * process.dilep * process.zmass)

    #this is for the MVA

    process.load('EgammaAnalysis.ElectronTools.electronIdMVAProducer_cfi')
    process.mvaID = cms.Sequence(  process.mvaTrigV0 + process.mvaTrigNoIPV0 + process.mvaNonTrigV0 )

    #Electron ID
    process.patElectrons.electronIDSources = cms.PSet(
        #MVA
        mvaTrigV0 = cms.InputTag("mvaTrigV0"),
        mvaNonTrigV0 = cms.InputTag("mvaNonTrigV0"),
        mvaTrigNoIPV0 = cms.InputTag("mvaTrigNoIPV0"),
        )
    
    process.patElectrons.addElectronID = cms.bool(True)
    #add pat conversions
    process.patConversions = cms.EDProducer("PATConversionProducer",
                                            # input collection
                                            #electronSource = cms.InputTag("gsfElectrons"),
                                            electronSource = cms.InputTag("cleanPatElectrons")  
                                            #electronSource = cms.InputTag("selectedPatElectrons")  
                                            # this should be your last selected electron collection name since currently index is used to match with electron later. We can fix this using reference pointer. ,

                                            )
    ## let it run
    process.p = cms.Path(
        process.mvaID *
        process.patDefaultSequence *
        process.patConversions *
        process.diLepSelSeq *
        (process.kt6PFJets*process.patJetCorrFactors) 
        )
    common_filters(process)
    #print process.p
    
    process.out.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('p')
        )

    #Keep NVtxs
    process.out.outputCommands.append('keep *_offlinePrimaryVertices_*_*')
    
    
    process.patJetCorrFactors.rho = cms.InputTag("kt6PFJetsPFlow", "rho")
    process.out.outputCommands.append('keep *_kt6PFJets_rho_PAT')
    process.out.outputCommands.append('keep patElectrons_selected*PFlowLoose*_*_*')
    process.out.outputCommands.append('keep *_cleanPatElectrons_*_*')
    process.out.outputCommands.append('keep *_diel_*_*')
    process.out.outputCommands.append('keep *_dimu_*_*')
    process.out.outputCommands.append('keep *_mvaNonTrigV0_*_PAT')
    process.out.outputCommands.append('keep *_mvaTrigNoIPV0_*_PAT')
    process.out.outputCommands.append('keep *_mvaTrigV0_*_PAT')
    
#    process.out.outputCommands.append('keep *_*_*_*')
    
