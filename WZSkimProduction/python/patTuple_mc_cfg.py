def mc_config(process, cms):
### Prune the GEN particle collection ###
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.prunedGenParticles = cms.EDProducer("GenParticlePruner",
            src = cms.InputTag("genParticles"), select = cms.vstring(
           "drop  *",
           #keeps all particles from the hard matrix element
           "keep status = 3",
           #keeps all stable muons (13) and electrons (11) + neutrinos (14, 12)
           # + Z (23) + W (24) + W' (34)  and their (direct) mothers.
           "+keep (abs(pdgId) = 11 | abs(pdgId) = 13 | abs(pdgId) = 12 | abs(pdgId) = 14 | abs(pdgId) = 23 | abs(pdgId) = 24 | abs(pdgId) = 34) & status = 1"
    )
)
    process.GlobalTag.globaltag = 'START53_V7A::All'
    process.p *= process.prunedGenParticles
