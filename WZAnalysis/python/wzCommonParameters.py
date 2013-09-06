import FWCore.ParameterSet.Config as cms
singleMuonPaths = [
    #'HLT_Mu11',            # from 2010 data
    #'HLT_Mu15_v*',
    #'HLT_Mu17_v*',
    #'HLT_Mu20_v*',
    #'HLT_Mu24_v*',
    #'HLT_Mu30_v*',
    #'HLT_Mu40_v*',         # 2.5e33 unprescaled
    #'HLT_Mu40_eta2p1_v*',
    #'HLT_IsoMu15_v*',
    #'HLT_IsoMu17_v*',
    #'HLT_IsoMu24_v*',
    #'HLT_IsoMu30_v*',
    #'HLT_IsoMu40_v*',
    #'HLT_IsoMu24_eta2p1_v*',
]

doubleMuonPaths = [
    #'HLT_DoubleMu5_v*',
    #'HLT_DoubleMu7_v*',
    #'HLT_TripleMu5_v*',
    #'HLT_Mu13_Mu8_v*',     # 1e33 unprescaled
    #'HLT_Mu17_Mu8_v*',     # 4e33 unprescaled
    'HLT_DoubleMu7_v*',
    'HLT_Mu13_Mu8_v*',   #1e33 unprescaled
    'HLT_Mu17_Mu8_v*',    #3e33 unprescaled
    'HLT_Mu17_TkMu8_v*', #2011B
    'HLT_Mu22_TkMu8_v*', #2012 single l1 seeded
    'HLT_Mu22_TkMu22_v*',#2012
    ]

singleElectronPaths = [
    #'HLT_Ele17_SW_L1R',  # from 2010 data
    #'HLT_Ele17_SW_Isol_L1R_v*',
    #'HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    #'HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*',
    #'HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v*', #1e33 unprescaled (removed after end of June TS)
    #'HLT_Ele52_CaloIdVT_TrkIdT_v*',                  #1e33 unprescaled
    #'HLT_Ele65_CaloIdVT_TrkIdT_v*',                  #3e33 unprescaled
    #'HLT_Ele27_WP80_v*',
    #'HLT_Ele27_WP80_PFMET_MT50_v*',
    #'HLT_Ele80_CaloIdVT_TrkIdT_v*',
    ]

doubleElectronPaths = [
    #'HLT_DoubleEle17_SW_L1R_v*',
    #'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    #'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*', 	# variant used in Summer11 MC and earlier data!
    #'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',  # data
    #'HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v*', # for T&P
    'HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v*',
    'HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v*',#early version, MC
    'HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v*',
    'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*',
    'HLT_DoubleEle33_CaloIdT_v*',
  

    ]

# For special processing of MuEG datasets
muegPaths = [
    #'HLT_Mu5_Ele13_v2', # available in Spring11 MC
    #'HLT_Mu5_Ele8_CaloIdT_TrkIdVL_Ele8_CaloIdL_TrkIdVL_v*',
    #'HLT_Mu8_Ele17_CaloIdL_v*',
    #'HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v*'
    ]

hltUpdatedForRun = cms.int32(173692)
