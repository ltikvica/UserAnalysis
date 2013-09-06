// -*- C++ -*-
//
// Package:    UserAnalysis/WZAnalysis
// Class:      WZTreeMaker
// 
/**\class WZTreeMaker WZTreeMaker.cc UserAnalysis/WZAnalysis/src/WZTreeMaker.cc

Description: Analyzer for associated W and Z boson production

Implementation:
Documentation available at https://twiki.cern.ch/twiki/bin/view/CMS/WZSoftwareFramework
*/
//
// Original Author:  Jeff Klukas, Srecko Morovic
//         Created:  Thu Jul 16 10:07:35 CDT 2009
// $Id: WZTreeMaker.cc,v 1.3 2013/04/26 08:54:37 ltikvica Exp $
//
//
#include "UserAnalysis/WZProducts/interface/WZProducts.h"

#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h" 

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/TriggerEvent.h" 
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"

//non-pat MET T1 corrections
#include "DataFormats/METReco/interface/PFMET.h" 

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TLorentzVector.h"


class WZTreeMaker : public edm::EDAnalyzer {

public:
  explicit WZTreeMaker(const edm::ParameterSet&);
  ~WZTreeMaker();

private:
  virtual void beginJob();
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endLuminosityBlock(const edm::LuminosityBlock &, 
                                  const edm::EventSetup &);
  virtual void beginRun(const edm::Run &, const edm::EventSetup &);
  virtual void endRun(const edm::Run &, const edm::EventSetup &);
  virtual void endJob();

  void clearMaps();
  void initializeTree();
  void fillTree();

  bool shouldMakeBranch(TString);

  void storeI(std::string, int  );
  void storeF(std::string, float);
  void storeToMaps(const reco::Candidate *, std::string);
  void storeToMaps(const pat::Electron &, wz::Point, std::string = "electron");
  void storeToMaps(const pat::Muon &, wz::Point, std::string = "muon");
  void storeToMaps(const pat::Jet &, int i);
  void storeToMaps(const reco::Track &, wz::Point);
  void storeToMaps(const pat::MET &, std::string);
  void storeToMaps(const reco::Vertex &);

  void storeGenSummary(const reco::GenParticleCollection *);
  void storeGenWZInfo(const reco::GenParticleCollection *);
  int  findTauDecayMode(const reco::Candidate *,const reco::Candidate *);

  void storeTriggerSummary(const edm::Event *, const edm::EventSetup *, const edm::TriggerResults *, const trigger::TriggerEvent *);
  void storeTriggerMatchingInfo(const edm::Event * iEvent, const edm::Handle<pat::TriggerEvent > *patTriggerEvent, edm::Handle<wz::ElectronV> *electrons, edm::Handle<wz::MuonV> *muons);


  void storePileup(const std::vector<PileupSummaryInfo> *, const edm::Event * iEvent);
  void storeZCand(const reco::Candidate *, wz::Point, const wz::ElectronV &, const wz::MuonV &);
  void storeWCand(const wz::WCandidate *, wz::Point, const wz::ElectronV &, const wz::MuonV &, std::string);
  void storeWZCands(const wz::ZCandidate *, const wz::WCandidate *, std::string);

  // ----------member data ---------------------------
  
  TFile * outputFile_;
  TTree * outputTree_;
  TH1F  * numEventsHist_;
  TH1D  * userWeightHist_;

  std::map< std::string, float               >       floatMap_;
  std::map< std::string, int                 >         intMap_;
  std::map< std::string, long                >         longMap_;
  std::map< std::string, std::vector<float> *> floatVectorMap_;
  std::map< std::string, std::vector<int  > *>   intVectorMap_;

  double rho_;
  double rhoCentral_;
  double rhoCentralNoPU_;

  std::string outputFileName_;
  std::string hltProcessName_;
  bool        debug_;
  std::vector<std::string> numEventsNames_;
  std::string weightProducerName_;

  HLTConfigProvider hltConfigProvider_;
  //edm::LumiReWeighting LumiWeights_;
  //edm::LumiReWeighting LumiWeightsTruth_;
  //edm::Lumi3DReWeighting LumiWeights3D_;

  edm::InputTag       vertexLabel_;
  edm::InputTag genParticlesLabel_;
  edm::InputTag      genInfoLabel_;
  edm::InputTag patTrigEventLabel_;
  edm::InputTag   hltResultsLabel_;
  edm::InputTag     hltEventLabel_;
  edm::InputTag       pileupLabel_;
  edm::InputTag          zllLabel_;
  edm::InputTag     electronLabel_;
  edm::InputTag         muonLabel_;
  edm::InputTag          jetLabel_;
  edm::InputTag        trackLabel_;

  edm::InputTag  metT1Label_;
  std::vector<edm::ParameterSet> metPSets_;

  std::vector<int        >       genIdsToStore_;
  std::vector<std::string>       contentToDrop_;
  std::vector<std::string>       contentToKeep_;
  std::vector<std::string>     triggersToStore_;
  std::vector<std::string>          hltMuPaths_;
  std::vector<std::string>         hltElePaths_;
  int                         hltUpdatedForRun_;

  std::vector<unsigned int> jetMatchingInfo_;

  std::vector<std::string> triggersToStoreAux_;
  std::vector<std::string> hltElePathsAux_;
  std::vector<std::string> hltMuPathsAux_;
  std::vector<std::string> triggersToStoreFull_;
  std::vector<std::string> hltElePathsFull_;
  std::vector<std::string> hltMuPathsFull_;
  std::vector<std::string> trNamesShort_;
  std::vector<std::string> trStoreNames_;
  std::vector<int> hltEleIsTracker_;
  std::vector<int> trsEleIsTracker_;
  std::vector<int> trMatchExactName_;
  std::vector<int> trigMatchInfoEle_pass_;
  std::vector<int> trigMatchInfoEle_prescale_;
  std::vector<int> trigMatchInfoMu_pass_;
  std::vector<int> trigMatchInfoMu_prescale_;

//  std::string MCPUFile_;
//  std::string DataPUFile_;
//  std::string DataPUFile3D_;
//  std::string DataPUFileTruth_;
//  std::string MCPUHist_;
//  std::string MCPUHist3D_;
//  std::string MCPUHistTruth_;
//  std::string DataPUHist_;
};



using namespace wz;
using namespace edm;
using namespace std;
using namespace reco;
using namespace trigger;



WZTreeMaker::WZTreeMaker(const ParameterSet& pSet) :

  outputFileName_(pSet.getUntrackedParameter<string>("outputFileName")),
  hltProcessName_(pSet.getUntrackedParameter<string>("hltProcessName")),
  debug_         (pSet.getUntrackedParameter<bool  >("debug"         )),
  numEventsNames_(pSet.getUntrackedParameter< vector<string> >("numEventsNames")),
  weightProducerName_(pSet.getUntrackedParameter<string>("weightProducerName")),

        vertexLabel_(pSet.getUntrackedParameter<InputTag>("primaryVertexTag")),
  genParticlesLabel_(pSet.getUntrackedParameter<InputTag>( "genParticlesTag")),
       genInfoLabel_(pSet.getUntrackedParameter<InputTag>(      "genInfoTag")),
  patTrigEventLabel_(pSet.getUntrackedParameter<InputTag>( "patTrigEventTag")),
    hltResultsLabel_(pSet.getUntrackedParameter<InputTag>(   "hltResultsTag")),
      hltEventLabel_(pSet.getUntrackedParameter<InputTag>(   "hltSummaryTag")),
        pileupLabel_(pSet.getUntrackedParameter<InputTag>(       "pileupTag")),
           zllLabel_(pSet.getUntrackedParameter<InputTag>(          "zllTag")),
      electronLabel_(pSet.getUntrackedParameter<InputTag>(     "electronTag")),
          muonLabel_(pSet.getUntrackedParameter<InputTag>(         "muonTag")),
           jetLabel_(pSet.getUntrackedParameter<InputTag>(          "jetTag")),
         trackLabel_(pSet.getUntrackedParameter<InputTag>(        "trackTag")),

     metT1Label_(pSet.getUntrackedParameter<edm::InputTag>(      "metT1Label")),
   metPSets_(pSet.getUntrackedParameter< vector<ParameterSet> >( "metPSets" )),

        genIdsToStore_(pSet.getParameter<vector<int   > >(   "genIdsToStore")),
        contentToDrop_(pSet.getParameter<vector<string> >(   "contentToDrop")),
        contentToKeep_(pSet.getParameter<vector<string> >(   "contentToKeep")),
      triggersToStore_(pSet.getParameter<vector<string> >( "triggersToStore")),
           hltMuPaths_(pSet.getParameter<vector<string> >(      "hltMuPaths")),
          hltElePaths_(pSet.getParameter<vector<string> >(     "hltElePaths")),
          hltUpdatedForRun_(pSet.getParameter<int >    ("hltUpdatedForRun")),
  trNamesShort_()

//  MCPUFile_(pSet.getParameter<string>("MCPUFile")),
//  DataPUFile_(pSet.getParameter<string>("DataPUFile")),
//  DataPUFile3D_(pSet.getParameter<string>("DataPUFile3D")),
//  DataPUFileTruth_(pSet.getParameter<string>("DataPUFileTruth")),
//  MCPUHist_(pSet.getParameter<string>("MCPUHist")),
//  MCPUHist3D_(pSet.getParameter<string>("MCPUHist3D")),
//  MCPUHistTruth_(pSet.getParameter<string>("MCPUHistTruth")),
//  DataPUHist_(pSet.getParameter<string>("DataPUHist")),

{

  if (debug_) cout << ">>>>>>>>>>>>>>>>>>> In constructor" << endl;
//  LumiWeights_ = edm::LumiReWeighting(MCPUFile_,DataPUFile_,
//                                      MCPUHist_,DataPUHist_);
//  LumiWeightsTruth_ = edm::LumiReWeighting(MCPUFile_,DataPUFileTruth_,
//                                      MCPUHistTruth_,DataPUHist_);
//  LumiWeights3D_ = edm::Lumi3DReWeighting(MCPUFile_,DataPUFile3D_,
//                                      MCPUHist3D_,DataPUHist_);
//  LumiWeights3D_.weight3D_init(1.0);
}



WZTreeMaker::~WZTreeMaker()
{
}



void
WZTreeMaker::analyze(const Event& iEvent, const EventSetup& iSetup)
{

  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // Begin Event
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  
  if (debug_) cout << ">>>>>>>>>>>>>>>>>>> In analyze" << endl;

  clearMaps();


  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // Get objects
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
 
  intMap_["runNumber"] = iEvent.id().run();
  intMap_["eventID"  ] = iEvent.id().event();
  longMap_["eventIDL" ] = iEvent.id().event();
  intMap_["lumiBlock"] = iEvent.id().luminosityBlock();

  Handle<VertexCollection      > primaryVerticesHandle;
  Handle<GenEventInfoProduct   > genInfoHandle        ;
  Handle<pat::TriggerEvent     > patTriggerEvent      ;
  Handle<TriggerEvent          > triggerEvent         ;
  Handle<TriggerResults        > triggerResults       ;
  Handle<GenParticleCollection > genParticles         ;
  Handle<std::vector< PileupSummaryInfo > >  PupInfo  ;
  Handle<ElectronV             > electrons            ;
  Handle<MuonV                 > muons                ;
  Handle<JetV                  > jets                 ;
  Handle<TrackCollection       > tracks               ;

  string hltName = hltConfigProvider_.processName();
  InputTag hltEventTag(hltEventLabel_.label(), "", hltName);
  InputTag hltResultsTag(hltResultsLabel_.label(), "", hltName);

  iEvent.getByLabel(      vertexLabel_, primaryVerticesHandle);
  iEvent.getByLabel(     genInfoLabel_, genInfoHandle        );
  iEvent.getByLabel(patTrigEventLabel_, patTriggerEvent      );
  iEvent.getByLabel(       hltEventTag, triggerEvent         );
  iEvent.getByLabel(     hltResultsTag, triggerResults       );

  iEvent.getByLabel("prunedGenParticles", genParticles         );
  // iEvent.getByLabel(genParticlesLabel_, genParticles         );
  iEvent.getByLabel(      pileupLabel_, PupInfo              );
  iEvent.getByLabel(    electronLabel_, electrons            );
  iEvent.getByLabel(        muonLabel_, muons                );
  iEvent.getByLabel(         jetLabel_, jets                 );
  iEvent.getByLabel(       trackLabel_, tracks               );

  Handle<double> rhoHandle;
 
  iEvent.getByLabel(InputTag("kt6PFJets", "rho"), rhoHandle);
  rho_ = * rhoHandle;
  floatMap_["PU_rhoCorr"] = rho_;
  try {
    //Handle<double> rhoCentralH;
  //iEvent.getByLabel(edm::InputTag("kt6PFJetsCentral", "rho", "PAT"), rhoCentralH);
  rhoCentral_ = -999;//* rhoCentralH;
  floatMap_["PU_rhoCorr_Central"] = rhoCentral_;

  //Handle<double> rhoCentralNoPUH;
  //iEvent.getByLabel(edm::InputTag("kt6PFJetsCentralNoPU", "rho", "PAT"), rhoCentralNoPUH);
  rhoCentralNoPU_ = -999; //* rhoCentralNoPUH;
  floatMap_["PU_rhoCorr_CentralNoPU"] = rhoCentralNoPU_;
  }
  catch (...) {}
  //met t1
  Handle<std::vector<reco::PFMET> > t1metHandle;
  iEvent.getByLabel(metT1Label_,t1metHandle);
  floatMap_["t1MET_et"]=t1metHandle->at(0).et();
  floatMap_["t1MET_phi"]=t1metHandle->at(0).phi();

  METV mets;
  vector<string> metNames;
  vector<string> wNames;
  for (size_t i = 0; i < metPSets_.size(); i++) {
    InputTag label = metPSets_[i].getParameter<InputTag>("tag");
    Handle<METV> handle;
    iEvent.getByLabel(label, handle);
    mets.push_back(handle->at(0));
    metNames.push_back(metPSets_[i].getParameter<string>("name"));
    if (metNames[i].substr(0, 3) == "met")
      wNames.push_back("W");
    else
      wNames.push_back(metNames[i].substr(0, 2) + "W");
  }

  if (debug_) cout << muons->size()      << " muons, " 
                   << electrons->size()  << " electrons, "
                   << jets->size()       << " jets, "
                   << endl;

  Point vertex = Point();
  if (primaryVerticesHandle.isValid()) {
    vertex = primaryVerticesHandle->begin()->position();
    for (size_t i = 0; i < primaryVerticesHandle->size(); i++)
      storeToMaps(primaryVerticesHandle->at(i));
  }
  else if (debug_) 
    cout << "No primary vertex found!" << endl;

  if (genInfoHandle.isValid()) {
    const vector<double> * binningVals = & genInfoHandle->binningValues();
    if (binningVals->size() > 0) {
      floatMap_["genEventScale"] = binningVals->at(0);
      if (debug_) cout << "genEventScale = " << binningVals->at(0) << endl;
    }
  }
  else if (debug_) cout << "No genEventScale" << endl;


  if (triggerResults.isValid() && triggerEvent.isValid())
    storeTriggerSummary(& iEvent, & iSetup, & * triggerResults, & * triggerEvent);

  if (debug_)
    std::cout << " ------------- begin Event: " 
      <<"runNumber=" << iEvent.id().run()
      << " eventID=" << iEvent.id().event() << std::endl;

  for (size_t i = 0; i < electrons->size(); i++) {
    storeToMaps(electrons->at(i), vertex);
  }
  
  for (size_t i = 0; i < muons->size(); i++) {
    storeToMaps(muons->at(i), vertex);
  }
  if (patTriggerEvent.isValid())
    storeTriggerMatchingInfo( &iEvent, &patTriggerEvent, &electrons, &muons);

  if (debug_)
    std::cout << " --------------end Event " << iEvent.id().event() << std::endl;

  for (size_t i = 0; i < jets->size(); i++) 
    storeToMaps(jets->at(i), i);
  // for (size_t i = 0; i < tracks->size(); i++) 
  //   if (tracks->at(i).pt() > 0.5)
  //     storeToMaps(tracks->at(i), vertex);
  for (size_t i = 0; i < mets.size(); i++)
    storeToMaps(mets[i], metNames[i]);

  if (genParticles.isValid()) { 
    storeGenSummary(& * genParticles);
    storeGenWZInfo(& * genParticles);
  }

  if (PupInfo.isValid()) {
    storePileup(& * PupInfo, & iEvent);
  }

  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // Reconstruct the W, Z, and resonance
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

  // Get only leptons passing ID
  MuonV looseMuons, tightMuons;
  for (size_t i = 0; i < muons->size(); i++) {
    if (muons->at(i).userInt("wp95") == 7)
      looseMuons.push_back(muons->at(i));
    if (muons->at(i).userInt("wp80") == 5 || muons->at(i).userInt("wp80") == 7)
      tightMuons.push_back(muons->at(i));
  }  
  ElectronV looseElectrons, tightElectrons;
  for (size_t i = 0; i < electrons->size(); i++) {
    bool skipElectron = false;
    for (size_t j = 0; j < looseMuons.size(); j++) {
      if (deltaR(electrons->at(i).eta(), electrons->at(i).phi(),
                looseMuons.at(j).eta(), looseMuons.at(j).phi()) < 0.01)
        skipElectron = true;
    }
    if (skipElectron) continue;
    if (electrons->at(i).userInt("wp95") == 7)
      looseElectrons.push_back(electrons->at(i));
    if (electrons->at(i).userInt("wp80") == 5 || electrons->at(i).userInt("wp80") == 7)
      tightElectrons.push_back(electrons->at(i));
  }

  ZCandV zCands = getZCands(looseElectrons, looseMuons, 60., 120.);

  intMap_["numberOfZs"] = zCands.size();
  if (zCands.size() > 0) {
    storeZCand(& zCands[0], vertex, * electrons, * muons);
    for (size_t i = 0; i < wNames.size(); i++) {
      WCandidate wCand = getWCand(tightElectrons, tightMuons, mets[i], zCands[0]);
      if (wCand) {
        storeWCand(& wCand, vertex, * electrons, * muons, wNames[i]);
        storeWZCands(& zCands[0], & wCand, wNames[i]);
      }
    }
  }

  fillTree();
  if (debug_) cout << ">>>>>>>>>>>>>>>>>>> Out analyze" << endl;

}


void 
WZTreeMaker::clearMaps()
{

  map< string, int            >::iterator  imi;
  map< string, long           >::iterator  lmi;
  map< string, float          >::iterator  fmi;
  map< string, vector<int  > *>::iterator ivmi;
  map< string, vector<float> *>::iterator fvmi;

  for (imi =   intMap_.begin(); imi !=   intMap_.end(); imi++) 
    imi->second = 0 ;

  for (lmi =   longMap_.begin(); lmi != longMap_.end(); lmi++) 
    lmi->second = 0 ;
 
  for (fmi = floatMap_.begin(); fmi != floatMap_.end(); fmi++) 
    fmi->second = 0.;

  for (ivmi =   intVectorMap_.begin(); ivmi !=   intVectorMap_.end(); ivmi++) 
    ivmi->second->clear();
  for (fvmi = floatVectorMap_.begin(); fvmi != floatVectorMap_.end(); fvmi++) 
    fvmi->second->clear();

}



void 
WZTreeMaker::initializeTree()
{

  intMap_["runNumber" ] = 0;
  intMap_["eventID"   ] = 0;
  longMap_["eventIDL" ] = 0;
  intMap_["lumiBlock" ] = 0;
  intMap_["numberOfZs"] = 0;

  floatMap_["PU_rhoCorr"   ] = 0.;
  floatMap_["PU_rhoCorr_Central"   ] = 0.;
  floatMap_["PU_rhoCorr_CentralNoPU"] = 0.;
  floatMap_["genEventScale"] = 0.;
  floatMap_["weight"       ] = 1.;

  floatMap_["t1MET_et"]=0;
  floatMap_["t1MET_phi"]=0;


  storeToMaps(pat::Muon(), Point());
  storeToMaps(pat::Electron(), Point());
  storeToMaps(pat::Jet(),-1);
  //storeToMaps(reco::Track(), Point());

  for (size_t i = 0; i < metPSets_.size(); i++)
    storeToMaps(pat::MET(), metPSets_[i].getParameter<string>("name"));
  storeToMaps(reco::Vertex());

  storeGenSummary(0);
  storeGenWZInfo(0);
  storeTriggerSummary(0,0,0,0);
  storeTriggerMatchingInfo(0,0,0,0);
  storePileup(0,0);

  storeZCand  (0, Point(), ElectronV(), MuonV());
  storeWCand  (0, Point(), ElectronV(), MuonV(), "W");
  storeWCand  (0, Point(), ElectronV(), MuonV(), "pfW");
  storeWCand  (0, Point(), ElectronV(), MuonV(), "tcW");
  storeWZCands(0, 0, "W");
  storeWZCands(0, 0, "pfW");
  storeWZCands(0, 0, "tcW");

  clearMaps();

  map< string, int            >::iterator  imi;
  map< string, long           >::iterator  lmi;
  map< string, float          >::iterator  fmi;
  map< string, vector<int  > *>::iterator ivmi;
  map< string, vector<float> *>::iterator fvmi;

  // Speeds up searching if all branches are to be kept
  if (contentToDrop_.size() == 0) contentToKeep_.resize(1, ".*");

  for (imi =   intMap_.begin(); imi !=   intMap_.end(); imi++) {
    string name     = imi->first.c_str();
    string leafList = name + "/I";
    if (shouldMakeBranch(name))
      outputTree_->Branch(name.c_str(), &(imi->second), leafList.c_str());
  }


  for (lmi =   longMap_.begin(); lmi !=   longMap_.end(); lmi++) {
    string name     = lmi->first.c_str();
    string leafList = name + "/L";
    if (shouldMakeBranch(name))
      outputTree_->Branch(name.c_str(), &(lmi->second), leafList.c_str());
  }


  for (fmi = floatMap_.begin(); fmi != floatMap_.end(); fmi++) {
    string name     = fmi->first.c_str();
    string leafList = name + "/F";
    if (shouldMakeBranch(name))
      outputTree_->Branch(name.c_str(), &(fmi->second), leafList.c_str());
  }

  for (ivmi =   intVectorMap_.begin(); ivmi !=   intVectorMap_.end(); ivmi++)
    if (shouldMakeBranch(ivmi->first.c_str()))
      outputTree_->Branch(ivmi->first.c_str(), &(ivmi->second));

  for (fvmi = floatVectorMap_.begin(); fvmi != floatVectorMap_.end(); fvmi++)
    if (shouldMakeBranch(fvmi->first.c_str()))
      outputTree_->Branch(fvmi->first.c_str(), &(fvmi->second));

}



void 
WZTreeMaker::fillTree()
{

  outputTree_->Fill();

}



bool
WZTreeMaker::shouldMakeBranch(TString branchName)
{

  for (size_t i = 0; i < contentToKeep_.size(); i++) {
    TPRegexp regexp = TPRegexp(contentToKeep_[i]);
    if (branchName.Contains(regexp)) return true;
  }

  for (size_t i = 0; i < contentToDrop_.size(); i++) {
    TPRegexp regexp = TPRegexp(contentToDrop_[i]);
    if (branchName.Contains(regexp)) return false;
  }

  return true;

}



void
WZTreeMaker::storeI(string key, int value)
{ 

  if (intVectorMap_[key] == 0) intVectorMap_[key] = new vector<int>;
  intVectorMap_[key]->push_back(value); 

}



void
WZTreeMaker::storeF(string key, float value)
{ 

  if (floatVectorMap_[key] == 0) floatVectorMap_[key] = new vector<float>;
  floatVectorMap_[key]->push_back(value); 

}



void
WZTreeMaker::storeToMaps(const Candidate * p, std::string prefix)
{

  string pre = prefix + "_";

  bool store;
  const bool isNull = (p == 0);

  store = !isNull;
  storeI(pre + "pdgId" , store ? p->pdgId()  : 0 );
  storeF(pre + "pt"    , store ? p->pt()     : 0.);
  storeF(pre + "et"    , store ? p->et()     : 0.);
  storeF(pre + "eta"   , store ? p->eta()    : 0.);
  storeF(pre + "phi"   , store ? p->phi()    : 0.);
  storeF(pre + "px"    , store ? p->px()     : 0.);
  storeF(pre + "py"    , store ? p->py()     : 0.);
  storeF(pre + "pz"    , store ? p->pz()     : 0.);
  storeF(pre + "mass"  , store ? p->mass()   : 0.);
  storeF(pre + "energy", store ? p->energy() : 0.);
  storeF(pre + "charge", store ? p->charge() : 0.);
}



void
WZTreeMaker::storeToMaps(const pat::Electron & p, Point vertex, std::string prefix)
{

  bool store;
  std::string pre = prefix + "_";

  storeToMaps(           &p             , prefix              );
  //storeToMaps(            p.genLepton() , prefix + "Gen"      );
  //storeToMaps( findMother(p.genLepton()), prefix + "GenMother");

  float pIn  = p.trackMomentumAtVtx().R();
  float pOut = p.trackMomentumOut  ().R();

  storeF(pre + "chargedHadronIso", p.chargedHadronIso()    );
  storeF(pre + "neutralHadronIso", p.neutralHadronIso()    );
  storeF(pre + "photonIso"       , p.photonIso()           );
  storeF(pre + "caloIso"      , p.caloIso()                       );
  storeF(pre + "ecalIso"      , p.ecalIso()                       );
  storeF(pre + "hcalIso"      , p.hcalIso()                       );
  storeF(pre + "pfChargedHadrons04",  p.userFloat("pfChargedHadrons04")   );
  storeF(pre + "pfChargedHadrons03",p.userFloat("pfChargedHadrons03") );
  storeF(pre + "pfNeutralHadrons04",  p.userFloat("pfNeutralHadrons04")   );
  storeF(pre + "pfNeutralHadrons03",p.userFloat("pfNeutralHadrons03") );
  storeF(pre + "pfPhotons04",  p.userFloat("pfPhotons04")                 );
  storeF(pre + "pfPhotons03",p.userFloat("pfPhotons03")               );
  storeF(pre + "pfAEff03",p.userFloat("pfAEff04"));
  storeF(pre + "pfAEff04",p.userFloat("pfAEff03"));
  storeF(pre + "trackIso"     , p.trackIso()                      );
  storeF(pre + "eOverP"       , p.eSuperClusterOverP()            );
  storeF(pre + "hOverE"       , p.hadronicOverEm()                );
  storeF(pre + "pIn"          , pIn                               ); 
  storeF(pre + "pOut"         , pOut                              ); 
  storeF(pre + "fBrem"        , p.fbrem()                         );
  storeF(pre + "deltaEtaIn"   , p.deltaEtaSuperClusterTrackAtVtx());
  storeF(pre + "deltaPhiIn"   , p.deltaPhiSuperClusterTrackAtVtx());
  storeF(pre + "eSeedOverPOut", p.eSeedClusterOverPout()          );
  storeF(pre + "sigmaEtaEta"  , p.scSigmaEtaEta()                 );
  storeF(pre + "sigmaIEtaIEta", p.sigmaIetaIeta()                 );
  storeF(pre + "convDist"     , p.convDist()                      );
  storeF(pre + "convDcot"     , p.convDcot()                      );
  


  //additional MVA variables
  storeF(pre + "kfChi2", p.userFloat("fMVAVar_kfchi2"));
  storeF(pre + "kfHits", p.userFloat("fMVAVar_kfhits"));
  storeF(pre + "kfHitsAll", p.userFloat("fMVAVar_kfhitsall"));
  storeF(pre + "deltaEtaCalo"   , p.deltaEtaSeedClusterTrackAtCalo());
  storeF(pre + "deltaPhiCalo"   , p.deltaPhiSeedClusterTrackAtCalo());
  storeF(pre + "sigmaIPhiIPhi", p.userFloat("fMVAVar_sigmaIPhiIPhi"));
  storeF(pre + "sigmaIEtaIPhi", p.userFloat("fMVAVar_sigmaIEtaIPhi"));
  storeF(pre + "etaWidth", p.core().isNonnull() ? p.superCluster()->etaWidth():-999); 
  storeF(pre + "phiWidth", p.core().isNonnull() ? p.superCluster()->phiWidth():-999); 
  storeF(pre + "e1x5e5x5", p.userFloat("fMVAVar_e1x5e5x5"));
  storeF(pre + "R9", p.userFloat("R9"));
  storeF(pre + "nBrems", p.core().isNonnull() ? p.numberOfBrems():-999);
  storeF(pre + "IoEmIoP", p.userFloat("fMVAVar_IoEmIoP"));
  storeF(pre + "eleEoPout",p.userFloat("fMVAVar_eleEoPout"));
  storeF(pre + "PreShowerOverRaw",p.userFloat("fMVAVar_PreShowerOverRaw"));
  storeF(pre + "EoPin",p.eSeedClusterOverP());
  storeF(pre + "mvaD0",p.userFloat("fMVAVAR_d0"));
  storeF(pre + "mvaIp3d",p.userFloat("fMVAVAR_ip3d"));
  storeF(pre + "mvaIp3dErr",p.userFloat("fMVAVAR_ip3dErr"));

  //2012 & HZZ working point MVA result
  storeF(pre + "BDT_MVA_NonTrig",p.userFloat("BDT_MVA_NonTrig"));
  storeF(pre + "BDT_MVA_Trig",p.userFloat("BDT_MVA_Trig"));
  storeI(pre + "MVA_pass_HZZ",p.userInt("MVA_pass_HZZ")); 
  storeI(pre + "MVA_pass_HZZ_relaxed",p.userInt("MVA_pass_HZZ_relaxed"));
  storeI(pre + "passPreselMVA2012",p.userInt("passPreselMVA2012"));

  //2011 & HWW working point MVA result
  storeF(pre+"BDT_MVA_HWW2011_Trig",p.userFloat("BDT_MVA_HWW2011_Trig"));
  storeI(pre+"BDT_MVA_HWW2011_Trig_pass",p.userInt("BDT_MVA_HWW2011_Trig_pass"));
  storeI(pre+"passPreselMVA2011",p.userInt("passPreselMVA2011"));

  //uncorrected MVA and variables (selection and preselection)
  //presel uncorrected specific
  storeF(pre+"gsf_dr03TkSumPt",p.userFloat("gsf_dr03TkSumPt"));
  storeF(pre+"gsf_dr03EcalRechitSumEt",p.userFloat("gsf_dr03EcalRechitSumEt"));
  storeF(pre+"gsf_dr03HcalTowerSumEt",p.userFloat("gsf_dr03HcalTowerSumEt"));

  //new MVA  
  storeF(pre+"mvaTrig",p.electronID("mvaTrigV0"));

  //MVA vars uncorrected gsf
  storeF(pre+"gsf_pt",p.userFloat("gsf_pt"));
  storeF(pre+"gsf_p",p.userFloat("gsf_p"));
  storeF(pre+"gsf_et",p.userFloat("gsf_et"));
  storeF(pre+"gsf_e",p.userFloat("gsf_e"));
  storeF(pre+"gsf_eSC",p.userFloat("gsf_eSC"));
  storeF(pre+"gsf_etaSC",p.userFloat("gsf_etaSC"));
  storeF(pre+"gsf_sigmaIEtaIEta",p.userFloat("gsf_sigmaIEtaIEta"));
  storeF(pre+"gsf_dEtaIn",p.userFloat("gsf_dEtaIn"));
  storeF(pre+"gsf_dPhiIn",p.userFloat("gsf_dPhiIn"));
  storeF(pre+"gsf_hoe",p.userFloat("gsf_hoe"));
  storeF(pre+"gsf_d0vtx",p.userFloat("gsf_d0vtx"));
  storeF(pre+"gsf_fbrem",p.userFloat("gsf_fbrem"));
  storeF(pre+"gsf_eSCOverP",p.userFloat("gsf_eSCOverP"));
  storeF(pre+"gsf_eSeedOverPout",p.userFloat("gsf_eSeedOverPout"));
  storeF(pre+"gsf_sigmaIPhiIPhi",p.userFloat("gsf_sigmaIPhiIPhi"));
  storeI(pre+"gsf_nBrems",p.userInt(   "gsf_nBrems"));
  storeF(pre+"gsf_IoEmIoP",p.userFloat("gsf_IoEmIoP"));
  storeF(pre+"gsf_eSeedOverP",p.userFloat("gsf_eSeedOverP"));
  storeF(pre+"gsf_ip3d",p.userFloat("gsf_ip3d"));
  storeF(pre+"gsf_ip3dErr",p.userFloat("gsf_ip3dErr"));

  storeF(pre+"gsf_BDT_MVA_HWW2011_Trig",p.userFloat("gsf_BDT_MVA_HWW2011_Trig"));
  storeI(pre+"gsf_BDT_MVA_HWW2011_Trig_pass",p.userInt("gsf_BDT_MVA_HWW2011_Trig_pass"));
  storeI(pre+"gsf_passPreselMVA2011",p.userInt("gsf_passPreselMVA2011"));

  // Store trigger matching info
  store = (p.triggerObjectMatches().size() > 0);
  storeF(pre + "hltMatchEt", store ? p.triggerObjectMatches()[0].et() : 0.);

  float hltMatchEtCluster = 0;
  float hltMatchEtElectron = 0;

  for (size_t i = 0; i < hltElePaths_.size(); i++) {
    string shortpath = hltElePaths_[i];
    if (shortpath.find("_v*") != string::npos)
        shortpath = shortpath.substr(0, shortpath.length() - 3);
    string longpath = i<hltElePathsFull_.size() ? hltElePathsFull_[i] : hltElePaths_[i];
    size_t nmatches = p.triggerObjectMatchesByPath(longpath, true, false).size();
    storeI(pre + "match" + shortpath, nmatches);
    if (nmatches) {
      if (hltEleIsTracker_.at(i)) hltMatchEtElectron = p.triggerObjectMatchesByPath(longpath, true, false)[0].et();
      else hltMatchEtCluster = p.triggerObjectMatchesByPath(longpath, true, false)[0].et();
    }
  }
  for (size_t i = 0; i < triggersToStore_.size(); i++) {
    string shortpath = triggersToStore_[i];
    if (shortpath.find("_v*") != string::npos)
        shortpath = shortpath.substr(0, shortpath.length() - 3);
    string longpath = i<triggersToStoreFull_.size() ? triggersToStoreFull_[i] : triggersToStore_[i];
    size_t nmatches = p.triggerObjectMatchesByPath(longpath, true, false).size();
    storeI(pre + "match" + shortpath, nmatches);
    if (nmatches) {
      if (trsEleIsTracker_.at(i)) hltMatchEtElectron = p.triggerObjectMatchesByPath(longpath, true, false)[0].et();
      else hltMatchEtCluster = p.triggerObjectMatchesByPath(longpath, true, false)[0].et();
    }
  }
  storeF(pre + "hltMatchEtCluster", hltMatchEtCluster);
  storeF(pre + "hltMatchEtElectron", hltMatchEtElectron);

  // Store Electron ID
  store = p.core().isNonnull(); // not a default-constructed electron
  storeI(pre+"ecalDrivenSeed",store? p.ecalDrivenSeed() : -1);
  string s;
  s = "wp60"; storeI(pre + s, store ? p.userInt(s) : 0);
  s = "wp70"; storeI(pre + s, store ? p.userInt(s) : 0);
  s = "wp80"; storeI(pre + s, store ? p.userInt(s) : 0);
  s = "wp85"; storeI(pre + s, store ? p.userInt(s) : 0);
  s = "wp90"; storeI(pre + s, store ? p.userInt(s) : 0);
  s = "wp95"; storeI(pre + s, store ? p.userInt(s) : 0);

  storeI(pre + "isEB", store ? p.isEB() :0); 
  storeI(pre + "isEE", store ? p.isEE() :0); 
  // Number of expected inner hits (for conversion rejection)
  int nHits = -10;
  if (p.core().isNonnull() && p.gsfTrack().isNonnull())
    nHits = p.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
  storeI(pre + "numHits" , nHits);

  int nLostHits = -999;
  if (p.core().isNonnull() && p.gsfTrack().isNonnull())
    nLostHits = p.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
  storeI(pre + "numLostHits" , nLostHits);
  storeI(pre + "gsf_numLostHits",p.userInt("gsf_nLostHits"));

  //  storeI(pre+"passConvVeto2012", store ? p.userInt("passConvVeto2012"):-1);
  storeI(pre+"passConvVeto2012", store ? p.passConversionVeto():-1);
  storeI(pre+"gsf_passConvVeto2012", store ? p.userInt("gsf_passConvVeto2012"):-1);
  storeF(pre+"hcalDepth1",store ? p.userFloat("hcalDepth1"):0);
  storeF(pre+"hcalDepth2",store ? p.userFloat("hcalDepth2"):0);
  storeF(pre+"hOverE2012",store ? p.userFloat("hOverE2012"):0);
  storeF(pre+"hcalIso2012",store ? p.userFloat("hcalIso2012"):0);

  double gsfTrackD0 = 0;
  double gsfTrackDz = 0;
  double gsfTrackD0Error = 0;
  double gsfTrackDzError = 0;

  //alternate calculation
  double gsfD0vtx = 0;
  double gsfDzvtx = 0;

  if (p.core().isNonnull() && p.gsfTrack().isNonnull()) {
    gsfTrackD0 = p.gsfTrack()->dxy(vertex) * (-1);
    gsfTrackD0Error = p.gsfTrack()->d0Error();
    gsfTrackDz = p.gsfTrack()->dz(vertex) * (-1);
    gsfTrackDzError = p.gsfTrack()->dzError();

    gsfD0vtx = p.gsfTrack()->d0() - vertex.x() * sin(p.gsfTrack()->phi()) + vertex.y() * cos(p.gsfTrack()->phi());
    gsfDzvtx = (p.vz() - vertex.z()) - ((p.vx() - vertex.x()) * p.px() + (p.vy() - vertex.y()) * p.py()) / p.pt() 
                     * p.pz()/p.pt();
  }
  storeF(pre + "gsfTrackD0" , gsfTrackD0);
  storeF(pre + "gsfTrackD0Error" , gsfTrackD0Error);
  storeF(pre + "gsfTrackDz" , gsfTrackDz);
  storeF(pre + "gsfTrackDzError" , gsfTrackDzError);
  storeF(pre + "gsfTrackD0vtx", gsfD0vtx);
  storeF(pre + "gsfTrackDzvtx", gsfDzvtx);

  double scEt,scEta;
  scEt = scEta = -999;
  if (p.core().isNonnull() ) {
      const reco::SuperClusterRef sc = p.superCluster();
      scEt = sc.get() ? sc.get()->energy() / cosh(sc.get()->eta()) : -1.;
      scEta = sc.get() ? sc.get()->eta() : -999;
  }
  storeF(pre+"ScEt", scEt);
  storeF(pre+"ScEta", scEta);

  storeF(pre + "dr03TkSumPt",store ? p.dr03TkSumPt():0.);
  storeF(pre + "dr03HcalTowerSumEt", store ? p.dr03HcalTowerSumEt() : 0.);
  storeF(pre + "dr03EcalRecHitSumEt",store ? p.dr03EcalRecHitSumEt() : 0. );
  storeF(pre + "combRelIso", store ? p.userFloat("combRelIso") : 0);
  storeF(pre + "modCombRelIso", store ? p.userFloat("modCombRelIso") : 0);
  storeF(pre + "corrEt", store ? p.userFloat("corrEt") : 0);

  storeI(pre + "wp12Veto03"  ,store? p.userInt("wp12Veto03"):0);
  storeI(pre + "wp12Loose03" ,store? p.userInt("wp12Loose03"):0);
  storeI(pre + "wp12Medium03",store? p.userInt("wp12Medium03"):0);
  storeI(pre + "wp12Tight03" ,store? p.userInt("wp12Tight03"):0);
  storeI(pre + "wp12Veto04"  ,store? p.userInt("wp12Veto04"):0);
  storeI(pre + "wp12Loose04" ,store? p.userInt("wp12Loose04"):0);
  storeI(pre + "wp12Medium04",store? p.userInt("wp12Medium04"):0);
  storeI(pre + "wp12Tight04" ,store? p.userInt("wp12Tight04"):0);

  storeI(pre + "wp12VetoNoIso"  ,store? p.userInt("wp12VetoNoIso"):0);
  storeI(pre + "wp12LooseNoIso" ,store? p.userInt("wp12LooseNoIso"):0);
  storeI(pre + "wp12MediumNoIso",store? p.userInt("wp12MediumNoIso"):0);
  storeI(pre + "wp12TightNoIso" ,store? p.userInt("wp12TightNoIso"):0);


  if (store) {
    //extract iso deposit
    /*
    const pat::IsoDeposit * pfChargedDep = p.isoDeposit(pat::PfChargedHadronIso);
    const pat::IsoDeposit * pfNeutralDep = p.isoDeposit(pat::PfNeutralHadronIso);
    const pat::IsoDeposit * pfPhotonDep  = p.isoDeposit(pat::PfGammaIso);
    reco::isodeposit::Direction dir = reco::isodeposit::Direction(p.eta(),p.phi());

    reco::IsoDeposit::AbsVetos vetosCharged;
    reco::IsoDeposit::AbsVetos vetosNeutral;
    reco::IsoDeposit::AbsVetos vetosPhoton;
    reco::isodeposit::ConeVeto v1( dir, 0.01 );
    reco::isodeposit::ConeVeto v2( dir, 0.07 );
    reco::isodeposit::RectangularEtaPhiVeto v3(dir,-0.025,0.025,-0.5,0.5);

    vetosCharged.push_back(&v1);
    vetosNeutral.push_back(&v2);
    vetosPhoton.push_back(&v3);

    float  isovalue_charged = pfChargedDep->sumWithin(0.4,vetosCharged);
    float  isovalue_neutral = pfNeutralDep->sumWithin(0.4,vetosNeutral);
    float  isovalue_photon = pfPhotonDep->sumWithin(0.4,vetosPhoton);

    */
    storeF(pre+"pfIsoEleSmurf",p.userFloat("eleSmurfPF"));
    storeF(pre+"gsf_pfIsoEleSmurf",p.userFloat("eleSmurfPFgsf"));
  }
  else {
    storeF(pre+"pfIsoEleSmurf",0);
    storeF(pre+"gsf_pfIsoEleSmurf",0);
  } 


  
  if (debug_)
  if (p.pt()>10)
  std::cout << "electron pt=" << p.pt()
            << " passMVAPresel=" << p.userInt("passPreselMVA2011")
            << " MVA BDT var=" << p.userFloat("BDT_MVA_HWW2011_Trig")
            << " pfIsoEle=" << p.userFloat("eleSmurfPF")
            << std::endl;
  
}



void
WZTreeMaker::storeToMaps(const pat::Muon & p, Point vertex, std::string prefix)
{

  bool store;
  std::string pre = prefix + "_";

  storeToMaps(           &p             , prefix              );
  //storeToMaps(            p.genLepton() , prefix + "Gen"      );
  //storeToMaps( findMother(p.genLepton()), prefix + "GenMother");

  // storeF(pre+"t0"                  , p.t0()               );
  storeF(pre+"caloIso"             , p.caloIso()          );
  storeF(pre+"trackIso"            , p.trackIso()         );
  const reco::MuonPFIsolation & iso = p.pfIsolationR04();
  storeF(pre + "pfChargedHadrons", iso.sumChargedHadronPt);
  storeF(pre + "pfNeutralHadrons", iso.sumNeutralHadronEt);
  storeF(pre + "pfPhotons"       , iso.sumPhotonEt       );
  storeF(pre + "deltaBeta"       , iso.sumPUPt           );
  
  //storeF(pre + "pfChargedHadrons",p.userFloat("pfChargedHadrons") );
  //storeF(pre + "pfNeutralHadrons",p.userFloat("pfNeutralHadrons") );
  //storeF(pre + "pfPhotons",p.userFloat("pfPhotons") );
  //storeF(pre + "deltaBeta",p.userFloat("deltaBeta") );

  /*
  if (debug_)
  if (p.pt()>10)
  std::cout << "muon pt=" << p.pt() << " pfChargedHadrons=" << p.userFloat("pfChargedHadrons")
            << " pfNeutralHadrons=" << p.userFloat("pfNeutralHadrons") << " pfPhotons="<<p.userFloat("pfPhotons") << " deltabeta="
            << p.userFloat("deltaBeta") << std::endl;
  */

  storeF(pre+"caloCompatibility"   , p.caloCompatibility());
  storeF(pre+"segmentCompatibility", muon::segmentCompatibility(p));

  float isoEtEcalDr30 = p.isolationR03().emEt;
  float isoEtHcalDr30 = p.isolationR03().hadEt;
  float isoEtTrackerDr30 = p.isolationR03().sumPt;

  storeF(pre+"isoEtEcalDr30"   , isoEtEcalDr30);
  storeF(pre+"isoEtHcalDr30"   , isoEtHcalDr30);
  storeF(pre+"isoEtTrackerDr30", isoEtTrackerDr30);

  // Combined relative isolation with pileup correction.
  storeF(pre+"combRelIso",  p.userFloat("combRelIso"));
  storeF(pre+"modCombRelIso",  p.userFloat("modCombRelIso"));
  storeF(pre+"combRelPFIso",  p.userFloat("combRelPFIso"));
  storeF(pre+"combRelPFIsoCorr",  p.userFloat("combRelPFIsoCorr"));

  if      (p.isGlobalMuon()    ) storeI(pre+"fitType", 1);
  else if (p.isTrackerMuon()   ) storeI(pre+"fitType", 2);
  else if (p.isStandAloneMuon()) storeI(pre+"fitType", 3);
  else                           storeI(pre+"fitType", 0);

  storeI(pre+"isGlobal", p.isGlobalMuon());
  storeI(pre+"isTracker", p.isTrackerMuon());
  storeI(pre+"numGlobalMatches", p.numberOfMatches()); 
  storeI(pre+"numChambers", p.numberOfChambers());

  // Store trigger matching info
  store = p.triggerObjectMatches().size() > 0;
  storeF(pre + "hltMatchPt", store ? p.triggerObjectMatches()[0].pt() : 0.);
  for (size_t i = 0; i < hltMuPaths_.size(); i++) {
    string shortpath = hltMuPaths_[i];
    if (shortpath.find("_v*") != string::npos)
        shortpath = shortpath.substr(0, shortpath.length() - 3);
    string longpath = i<hltMuPathsFull_.size() ? hltMuPathsFull_[i] : hltMuPaths_[i];
    storeI(pre + "match" + shortpath,
           p.triggerObjectMatchesByPath(longpath, true, false).size());
  }
  for (size_t i = 0; i < triggersToStore_.size(); i++) {
    string shortpath = triggersToStore_[i];
    if (shortpath.find("_v*") != string::npos)
        shortpath = shortpath.substr(0, shortpath.length() - 3);
    string longpath = i<triggersToStoreFull_.size() ? triggersToStoreFull_[i] : triggersToStore_[i];
    size_t nmatches = p.triggerObjectMatchesByPath(longpath, true, false).size();
    storeI(pre + "match" + shortpath, nmatches);
  }

  // Store Muon ID
  store = (p.isGlobalMuon() || p.isStandAloneMuon());
  string s;
  s = "wp95"; storeI(pre + s, store ? p.userInt(s) : 0);
  s = "wp80"; storeI(pre + s, store ? p.userInt(s) : 0);
 
  store = p.isGlobalMuon();
  TrackRef global = store ? p.globalTrack()  : TrackRef();
  storeF(pre+"globalChi2"   , store ? global->chi2()              : 0.);
  storeF(pre+"globalNdof"   , store ? global->ndof()              : 0.);
  storeF(pre+"globalD0"     , store ? global->dxy(vertex) * (-1)  : 0.);
  storeF(pre+"globalD0Error", store ? global->d0Error()           : 0.);
  storeF(pre+"globalDz"     , store ? global->dz(vertex) * (-1)   : 0.);
  storeF(pre+"globalDzError", store ? global->dzError()           : 0.);
  storeI(pre+"numValidMuonHits", store ? global->hitPattern().numberOfValidMuonHits() : 0);
  storeI(pre+"globalNumTrackerHits",store ?  global->hitPattern().numberOfValidTrackerHits() : 0);
  storeI(pre+"globalNumTrackerLayers" , store ? global->hitPattern().trackerLayersWithMeasurement() : 0);
  storeI(pre+"globalNumPixelHits" , store ?  global->hitPattern().numberOfValidPixelHits()   : 0);
  storeI(pre+"globalNumStripHits" , store ?  global->hitPattern().numberOfValidStripHits()   : 0);

  //muon trajectory kink
  storeF(pre+"trkKink", store ? p.combinedQuality().trkKink:-1);

  store = (p.isGlobalMuon() || p.isTrackerMuon());
  TrackRef inner  = store ? p.innerTrack() : TrackRef();
  storeF(pre+"innerChi2"       , store ?  inner->chi2()              : 0.);
  storeF(pre+"innerNdof"       , store ?  inner->ndof()              : 0.);
  storeF(pre+"innerD0"         , store ?  inner->dxy(vertex) * (-1)  : 0.);
  storeF(pre+"innerD0Error"    , store ?  inner->d0Error()           : 0.);
  storeF(pre+"innerDz"         , store ?  inner->dz(vertex) * (-1)   : 0.);
  storeF(pre+"innerDzError"    , store ?  inner->dzError()           : 0.);
  storeI(pre+"numInnerHits"    , store ?  inner->numberOfValidHits() : 0 );
  storeI(pre+"numTrackerHits"  ,store  ?  inner->hitPattern().numberOfValidTrackerHits() : 0);
  storeI(pre+"numTrackerLayers",store  ?  inner->hitPattern().trackerLayersWithMeasurement() : 0);
  storeI(pre+"numPixelHits"    , store ?  inner->hitPattern().numberOfValidPixelHits()   : 0);
  storeI(pre+"numStripHits"    , store ?  inner->hitPattern().numberOfValidStripHits()   : 0);

  store = (p.isGlobalMuon() || p.isStandAloneMuon());
  TrackRef outer  = store ? p.outerTrack() : TrackRef();
  storeF(pre+"outerChi2"    , store ?  outer->chi2()              : 0.);
  storeF(pre+"outerNdof"    , store ?  outer->ndof()              : 0.);
  storeF(pre+"outerD0"      , store ?  outer->dxy(vertex) * (-1)  : 0.);
  storeF(pre+"outerD0Error" , store ?  outer->d0Error()           : 0.);
  storeF(pre+"outerDz"      , store ?  outer->dz(vertex) * (-1)    : 0.);
  storeF(pre+"outerDzError" , store ?  outer->dzError()           : 0.);
  storeI(pre+"numOuterHits" , store ?  outer->numberOfValidHits() : 0 );

  storeF(pre+"sip"       , p.userFloat("sip")       );
  storeI(pre+"nSegments" , p.userFloat("nSegments") );

}



void
WZTreeMaker::storeToMaps(const pat::Jet & p, int i)
{

  string pre = "jet_";

  storeToMaps(&p         , "jet"   );
  storeToMaps( p.genJet(), "jetGen");

  string bDiscTag = "simpleSecondaryVertexBJetTags";

  storeI(pre + "partonFlavour" , p.partonFlavour()         );
  storeF(pre + "bDiscriminator", p.bDiscriminator(bDiscTag));

  bool store = jetMatchingInfo_.size();
  storeI(pre + "trigMatchMask", store ? jetMatchingInfo_[i] : 0);

}



void
WZTreeMaker::storeToMaps(const reco::Track & p, Point vertex)
{

  string pre = "track_";

  storeF(pre + "pt", p.pt()              );
  storeF(pre + "d0", p.dxy(vertex) * (-1));
  
  const reco::HitPattern & pattern = p.hitPattern();
  storeI(pre + "nPixelHits", pattern.numberOfValidPixelHits()); 
  storeI(pre + "nStripHits", pattern.numberOfValidStripHits());

}



void
WZTreeMaker::storeToMaps(const pat::MET & p, string prefix)
{

  string pre = prefix + "_";

  floatMap_[pre + "et"    ] = p.et() ;
  floatMap_[pre + "phi"   ] = p.phi();
  
  if (prefix == "met") {
    floatMap_["metGen_et" ] = p.genMET() ? p.genMET()->et()  : 0.;
    floatMap_["metGen_phi"] = p.genMET() ? p.genMET()->phi() : 0.;
  }

}



void
WZTreeMaker::storeToMaps(const Vertex & p)
{

  string pre = "vertex_";

  bool store = p.isValid();

  storeI(pre + "isValid", store);
  storeF(pre + "x", store ? p.x() : 0.);
  storeF(pre + "y", store ? p.y() : 0.);
  storeF(pre + "z", store ? p.z() : 0.);
  storeF(pre + "xError", store ? p.xError() : 0.);
  storeF(pre + "yError", store ? p.yError() : 0.);
  storeF(pre + "zError", store ? p.zError() : 0.);
  storeF(pre + "nTracks", store ? p.nTracks() : 0.);
  storeF(pre + "ndof", store ? p.ndof() : 0.);
  storeF(pre + "chi2", store ? p.chi2() : 0.);

}



void
WZTreeMaker::storeZCand(const Candidate * p, Point vertex, ElectronV const & ev, MuonV const & mv)
{

  string pre = "Z_";

  const bool isNull = (p == 0);
  bool       store  = !isNull;

  int flavor = store ? abs(p->daughter(0)->pdgId()) : 0;
  vector<int> leptonIndices;
  
  if (!isNull && flavor == 11)
    for (size_t i = 0; i < ev.size(); i++)
      if (areOverlapping(* p, ev[i]))
        leptonIndices.push_back(i);

  if (!isNull && flavor == 13)
    for (size_t i = 0; i < mv.size(); i++)
      if (areOverlapping(* p, mv[i]))
        leptonIndices.push_back(i);

  store = (leptonIndices.size() == 2);
  intMap_  [pre + "pdgId"       ] = store ? p->pdgId()       : 0 ;
  intMap_  [pre + "flavor"      ] = store ? flavor           : 0 ;
  intMap_  [pre + "leptonIndex1"] = store ? leptonIndices[0] : 0 ;
  intMap_  [pre + "leptonIndex2"] = store ? leptonIndices[1] : 0 ;
  floatMap_[pre + "pt"          ] = store ? p->pt()          : 0.;
  floatMap_[pre + "eta"         ] = store ? p->eta()         : 0.;
  floatMap_[pre + "phi"         ] = store ? p->phi()         : 0.;
  floatMap_[pre + "mass"        ] = store ? p->mass()        : 0.;

}



void
WZTreeMaker::storeWCand(const WCandidate * p, Point vertex,
                        ElectronV const & ev, MuonV const & mv, string prefix)
{

  string pre = prefix + "_";

  const bool isNull = (p == 0);
  bool       store  = !isNull;

  const Candidate * wLep  = store ? p->daughter(0) : 0;
  int flavor = store ? abs(wLep->pdgId()) : 0;
  int charge = store ? wLep->pdgId() / abs(wLep->pdgId()) : 0;
  vector<int> leptonIndices;

  if (!isNull && flavor == 11)
    for (size_t i = 0; i < ev.size(); i++)
      if (areOverlapping(* p, ev[i]))
        leptonIndices.push_back(i);

  if (!isNull && flavor == 13)
    for (size_t i = 0; i < mv.size(); i++)
      if (areOverlapping(* p, mv[i]))
        leptonIndices.push_back(i);

  const int iZLep1  = intMap_["Z_leptonIndex1"];
  const int iZLep2  = intMap_["Z_leptonIndex2"];
  const int zFlavor = intMap_["Z_flavor"      ];
  double dR1 = 999.;
  double dR2 = 999.;

  if (zFlavor && wLep) {
    const Candidate * zLep1 = 0;
    const Candidate * zLep2 = 0;
    if (zFlavor == 11) zLep1 = & ev[iZLep1];
    if (zFlavor == 11) zLep2 = & ev[iZLep2];
    if (zFlavor == 13) zLep1 = & mv[iZLep1];
    if (zFlavor == 13) zLep2 = & mv[iZLep2];
    dR1 = deltaR(zLep1->eta(), zLep1->phi(), wLep->eta(), wLep->phi());
    dR2 = deltaR(zLep2->eta(), zLep2->phi(), wLep->eta(), wLep->phi());
  }

  store        = (leptonIndices.size() == 1);
  intMap_  [pre + "pdgId"       ] = store ? 24 * charge      : 0 ;
  intMap_  [pre + "flavor"      ] = store ? flavor           : 0 ;
  intMap_  [pre + "leptonIndex" ] = store ? leptonIndices[0] : 0 ;
  floatMap_[pre + "pt"          ] = store ? p->pt()          : 0.;
  floatMap_[pre + "eta"         ] = store ? p->eta()         : 0.;
  floatMap_[pre + "phi"         ] = store ? p->phi()         : 0.;
  floatMap_[pre + "transMass"   ] = store ? p->mt()          : 0.;
  floatMap_["Wlep_Zlep1_dR"     ] = store ? dR1              : 0.;
  floatMap_["Wlep_Zlep2_dR"     ] = store ? dR2              : 0.;

}



void
WZTreeMaker::storeWZCands(const ZCandidate * Z, const WCandidate * W, 
                          string prefix)
{

  const bool isNull = (Z == 0 || W == 0);

  WZCandidate wzCand;
  if (!isNull)
    wzCand = WZCandidate(* Z, * W);

  string pre;
  pre = prefix + "_neutrino_";
  floatMap_[pre + "pzMinAngle"] = wzCand.neutrinoPz("minAngle");
  floatMap_[pre + "pzMaxAngle"] = wzCand.neutrinoPz("maxAngle");
  floatMap_[pre + "pzMaxPz"   ] = wzCand.neutrinoPz("maxPz"   );
  floatMap_[pre + "pzMinPz"   ] = wzCand.neutrinoPz("minPz"   );
  pre = prefix + "Z_";
  floatMap_[pre + "transMass"      ] = wzCand.transMass();
  floatMap_[pre + "invMassMinAngle"] = wzCand.mass("minAngle");
  floatMap_[pre + "invMassMaxAngle"] = wzCand.mass("maxAngle");
  floatMap_[pre + "invMassMaxPz"   ] = wzCand.mass("maxPz"   );
  floatMap_[pre + "invMassMinPz"   ] = wzCand.mass("minPz"   );

}



void
WZTreeMaker::storeGenSummary(const GenParticleCollection * genParticles)
{

  const bool isNull = (genParticles == 0);
  if (isNull) {
    storeToMaps((Candidate *) 0, "genParticle");
    storeToMaps((Candidate *) 0, "genMother"  );
    //storeToMaps((Candidate *) 0, "genAncestor"  );
    return;
  }

  for (unsigned int i = 0; i < genParticles->size(); i++) {
    const Candidate *     p      = &genParticles->at(i);
    int                   pdgId  = abs(p->pdgId() );
    int                   status = abs(p->status());
    vector<int> &         ids    = genIdsToStore_;
    vector<int>::iterator iter   = std::find(ids.begin(), ids.end(), pdgId);
    if (iter != ids.end())
    if (status == 1 || ((abs(pdgId)==24 || pdgId==23) && status!=2)) {
      storeToMaps(           p , "genParticle");
      storeToMaps(findMother(p), "genMother"  );
      //storeToMaps(findAncestor(p), "genAncestor"  );
    }
  }
}

void
WZTreeMaker::storeGenWZInfo(const GenParticleCollection * genParticles)
{
  intMap_["genWZ_zLeptId1"]=0;
  intMap_["genWZ_zLept1TauDecay"]=0;
  intMap_["genWZ_zLeptId2"]=0;
  intMap_["genWZ_zLept2TauDecay"]=0;
  intMap_["genWZ_wLeptId"]=0;
  intMap_["genWZ_wLeptTauDecay"]=0;
  intMap_["genWZ_wNeutrinoId"]=0;

  floatMap_["genWZ_zLept1Pt"]=0;
  floatMap_["genWZ_zLept1Eta"]=0;
  floatMap_["genWZ_zLept1Phi"]=0;
  floatMap_["genWZ_zLept1Energy"]=0;
  intMap_["genWZ_zLept1Charge"]=0;

  floatMap_["genWZ_zLept2Pt"]=0;
  floatMap_["genWZ_zLept2Eta"]=0;
  floatMap_["genWZ_zLept2Phi"]=0;
  floatMap_["genWZ_zLept2Energy"]=0;
  intMap_["genWZ_zLept2Charge"]=0;

  floatMap_["genWZ_wLeptPt"]=0;
  floatMap_["genWZ_wLeptEta"]=0;
  floatMap_["genWZ_wLeptPhi"]=0;
  floatMap_["genWZ_wLeptEnergy"]=0;
  intMap_["genWZ_wLeptCharge"]=0;

  floatMap_["genWZ_wNeutrinoPt"]=0;
  floatMap_["genWZ_wNeutrinoEta"]=0;
  floatMap_["genWZ_wNeutrinoPhi"]=0;
  floatMap_["genWZ_wNeutrinoEnergy"]=0;

  floatMap_["genWZ_zLept1TauDecayPt"]=0;
  floatMap_["genWZ_zLept1TauDecayEta"]=0;
  floatMap_["genWZ_zLept1TauDecayPhi"]=0;
  floatMap_["genWZ_zLept2TauDecayPt"]=0;
  floatMap_["genWZ_zLept2TauDecayEta"]=0;
  floatMap_["genWZ_zLept2TauDecayPhi"]=0;
  floatMap_["genWZ_wLeptTauDecayPt"]=0;
  floatMap_["genWZ_wLeptTauDecayEta"]=0;
  floatMap_["genWZ_wLeptTauDecayPhi"]=0;

  floatMap_["genWZ_FSRZMass"]=0;
  floatMap_["genWZ_FSRZMass_wFake"]=0;

  if (!genParticles) return;

  bool zfound=false;
  bool zFoundConsistent=false;
  bool zfoundInMassRange=false;
  bool wfound=false;
  double zMass=0;
  double zMass_wFake=0;
  int zLeptType=0;
  int zLeptType1=0;
  int zLeptType2=0;
  int wLeptType=0;

  TLorentzVector zI,zJ;

  const Candidate * zLept1=0;
  const Candidate * zLept2=0;
  const Candidate * wLept=0;
  const Candidate * wNeutrino=0;

  vector<int> zLeptonIds_;
  vector<int> wLeptonIds_;

  zLeptonIds_.push_back(11);
  zLeptonIds_.push_back(13);
  //zLeptonIds_.push_back(15);
  wLeptonIds_.push_back(11);
  wLeptonIds_.push_back(13);
  //wLeptonIds_.push_back(15);

  //find signal (e/mu Z decay)
  for (size_t i=0; i<genParticles->size();i++) {
   const Candidate *p = &genParticles->at(i);
   const Candidate *mother = 0;
   int pdg=0;
   for (size_t j=0;j< zLeptonIds_.size();j++) {
     if (abs(zLeptonIds_[j])==abs(p->pdgId()) && p->status()==1) pdg=p->pdgId();
   }
   if (!pdg) continue;

   mother = p->mother();
   while (mother) {
     if (mother->pdgId()==p->pdgId()) {mother=mother->mother();continue;}
     if (mother->pdgId()==23) {
       zfound=true;
       zLeptType1=p->pdgId();
       zLept1=p;
       break;
     }
     break;
   }
   if (zfound) {
     for (size_t j=i+1;j<genParticles->size();j++) {
       const Candidate *r = &genParticles->at(j);
       //skip if not final state and same flavour
       if (r->status()!=1 ||r->pdgId()!=-p->pdgId()) continue;
       mother = r->mother();
       while (mother) {
	 if (mother->pdgId()==r->pdgId()) {mother=mother->mother();continue;}
	 if (mother->pdgId()==23) {
	   //found both leptons, reconstruct mass
	   zFoundConsistent=true;
	   TLorentzVector pi(p->p4().px(),p->p4().py(),p->p4().pz(),p->energy()) ;
	   TLorentzVector pj(r->p4().px(),r->p4().py(),r->p4().pz(),r->energy()) ;
           zI=pi;zJ=pj;
	   zMass = (pi+pj).M();
	   zLeptType=abs(pdg);
           zLeptType2=r->pdgId();
           zLept2=r;
	   break;
	 }
	 break;
       }
     }
     break;
   }
  } //Z lepton loop
  //consistency check and Z mass tange check
  if (zfound) {
    if (!zFoundConsistent) std::cout << " WARNING: only a single Z lepton found in this event\n";
    else {
      //check Z mass range
      if (zMass>60 && zMass<120) zfoundInMassRange=true;
      int zcount=0;
      bool consistent=false;
      for (size_t i=0;i<genParticles->size();i++) {
        const Candidate *p = &genParticles->at(i);
        if (p->pdgId()==23) {
            if (p->status()!=2
                && abs(p->daughter(0)->pdgId())==zLeptType && abs(p->daughter(1)->pdgId())==zLeptType) consistent=true;
            zcount++;
        }
      }
      //if (zcount>1) cout << "WARNING: multiple Z cands found\n";
      //if (!consistent) cout << "WARNING: error (bug) in MC Z search\n";
    }
  }
  else { //look for taus
    for (size_t i=0;i<genParticles->size();i++) {
      const Candidate *p = &genParticles->at(i);
      if (p->pdgId()==23) {
	if (p->status()!=2
	    && abs(p->daughter(0)->pdgId())==15 && abs(p->daughter(1)->pdgId())==15) {
              zLeptType1=p->daughter(0)->pdgId();
              zLeptType2=p->daughter(1)->pdgId();
              const reco::Candidate * fsrFromTau1 = 0;
              const reco::Candidate * fsrFromTau2 = 0;
              intMap_["genWZ_zLept1TauDecay"]=findTauDecayMode(p->daughter(0),fsrFromTau1);
              intMap_["genWZ_zLept2TauDecay"]=findTauDecayMode(p->daughter(1),fsrFromTau2);
              if (abs(intMap_["genWZ_zLept1TauDecay"])==13 || abs(intMap_["genWZ_zLept1TauDecay"])==11) {
                if (fsrFromTau1) {
                 floatMap_["genWZ_zLept1TauDecayPt"]=fsrFromTau1->pt();
                 floatMap_["genWZ_zLept1TauDecayEta"]=fsrFromTau1->eta();
                 floatMap_["genWZ_zLept1TauDecayPhi"]=fsrFromTau1->phi();
                }
              }
              if (abs(intMap_["genWZ_zLept2TauDecay"])==13 || abs(intMap_["genWZ_zLept2TauDecay"])==11) {
                if (fsrFromTau2) {
                 floatMap_["genWZ_zLept2TauDecayPt"]=fsrFromTau2->pt();
                 floatMap_["genWZ_zLept2TauDecayEta"]=fsrFromTau2->eta();
                 floatMap_["genWZ_zLept2TauDecayPhi"]=fsrFromTau2->phi();
                }
              }
              break;
        }
      }
    }
  }

  for (size_t i=0; i<genParticles->size();i++) {
   const Candidate *p = &genParticles->at(i);
   const Candidate *mother = 0;
   int pdg=0;
   for (size_t j=0;j< wLeptonIds_.size();j++) {
     if (abs(wLeptonIds_[j])==abs(p->pdgId()) && p->status()==1) pdg=p->pdgId();
     mother = p->mother();
   }
   if (!pdg) continue;
   //find W
   while (mother) {
     if (mother->pdgId()==pdg) {mother=mother->mother();continue;}
     if (mother->pdgId()==23) break;//do not let reach virtual W
     if (abs(mother->pdgId())==24) {
       wfound=true;wLeptType=p->pdgId();wLept=p;
       TLorentzVector wp(p->p4().px(),p->p4().py(),p->p4().pz(),p->energy());
       if (zLept1) { //build fake Z from W and Z lept
         if (zLept1->pdgId()==-p->pdgId())
           zMass_wFake=(wp+zI).M();
       }
       if (zLept2) {
         if (zLept2->pdgId()==-p->pdgId())
           zMass_wFake=(wp+zJ).M();
       }
     }
     break;
   }
  } //W lepton loop
  //look for W tau
  if (!wLeptType) {
    for (size_t i=0;i<genParticles->size();i++) {
      const Candidate *p = &genParticles->at(i);
      if (abs(p->pdgId())==24) {

	if ( p->status()!=2
	    && (
	      abs(p->daughter(0)->pdgId())==15 
	      || abs(p->daughter(1)->pdgId())==15
	      )
	   )
	{
	  if (abs(p->daughter(0)->pdgId())==15)
	    wLeptType=p->daughter(0)->pdgId();
	  else
	    wLeptType=p->daughter(1)->pdgId();

          //find tau decay mode
          const Candidate *ptau = 0;
          if (p->daughter(0) && abs(p->daughter(0)->pdgId())==15) {
            ptau = p->daughter(0);
          }
          else {
            ptau = p->daughter(1);
          }
          const reco::Candidate * fsrFromTau = 0;
          intMap_["genWZ_wLeptTauDecay"]=findTauDecayMode(ptau,fsrFromTau);

              if (abs(intMap_["genWZ_wLeptTauDecay"])==13 || abs(intMap_["genWZ_wLeptTauDecay"])==11) {
                if (fsrFromTau) {
                 floatMap_["genWZ_wLeptTauDecayPt"]=fsrFromTau->pt();
                 floatMap_["genWZ_wLeptTauDecayEta"]=fsrFromTau->eta();
                 floatMap_["genWZ_wLeptTauDecayPhi"]=fsrFromTau->phi();
                }
              }

	  break;
	}
      }
    }
  }
  //find neutrino
  if (wLeptType) {
    for (size_t i=0;i<genParticles->size();i++) {
      const Candidate *p = &genParticles->at(i);
      const Candidate *mother = 0;
      int pdg=0;
      if (abs(wLeptType)+1==abs(p->pdgId()) && p->status()==1) pdg=p->pdgId();
      else continue;
      mother = p->mother();
     while (mother) {
       if (mother->pdgId()==pdg) {mother=mother->mother();continue;}
       else if (abs(mother->pdgId())==24) {
         wNeutrino=p;
       }
       //else {/* not neutrino directly from W */}
       break;
     }
   }
  }

  intMap_["genWZ_zLeptId1"]=zLeptType1;
  intMap_["genWZ_zLeptId2"]=zLeptType2;
  intMap_["genWZ_wLeptId"]=wLeptType;
  floatMap_["genWZ_FSRZMass"]=zMass;
  floatMap_["genWZ_FSRZMass_wFake"]=zMass_wFake;
  if (zLept1) {
    floatMap_["genWZ_zLept1Pt"]=zLept1->pt();
    floatMap_["genWZ_zLept1Eta"]=zLept1->eta();
    floatMap_["genWZ_zLept1Phi"]=zLept1->phi();
    floatMap_["genWZ_zLept1Energy"]=zLept1->energy();
    floatMap_["genWZ_zLept1Charge"]=zLept1->charge();
  }
  if (zLept2) {
    floatMap_["genWZ_zLept2Pt"]=zLept2->pt();
    floatMap_["genWZ_zLept2Eta"]=zLept2->eta();
    floatMap_["genWZ_zLept2Phi"]=zLept2->phi();
    floatMap_["genWZ_zLept2Energy"]=zLept2->energy();
    floatMap_["genWZ_zLept2Charge"]=zLept2->charge();
  }
  if (wLept) {
    floatMap_["genWZ_wLeptPt"]=wLept->pt();
    floatMap_["genWZ_wLeptEta"]=wLept->eta();
    floatMap_["genWZ_wLeptPhi"]=wLept->phi();
    floatMap_["genWZ_wLeptEnergy"]=wLept->energy();
    floatMap_["genWZ_wLeptCharge"]=wLept->charge();
  }
  if (wNeutrino) {
    intMap_["genWZ_wNeutrinoId"] = wNeutrino->pdgId();
    floatMap_["genWZ_wNeutrinoPt"]=wNeutrino->pt();
    floatMap_["genWZ_wNeutrinoEta"]=wNeutrino->eta();
    floatMap_["genWZ_wNeutrinoPhi"]=wNeutrino->phi();
    floatMap_["genWZ_wNeutrinoEnergy"]=wNeutrino->energy();
  }

}

int WZTreeMaker::findTauDecayMode(const reco::Candidate *ptau, const reco::Candidate *found) {

  if (ptau) {
    while (ptau) {
      if (ptau->status()==3) ptau=ptau->daughter(0);
      else break;
    }
    int dID=0;
    while (ptau && ptau->daughter(dID)) {
    dID++;
    }
    dID=0;
    while (ptau && ptau->daughter(dID)) {
      int tdpdg=ptau->daughter(dID)->pdgId();

      //retry if W found
      if (abs(tdpdg)==24) {
        ptau=ptau->daughter(dID);
        dID=0;
        continue;
      }
      dID++;
      //skip neutrinos
      if (abs(tdpdg)==12 || abs(tdpdg)==14 || abs(tdpdg)==16) continue;
      found=ptau->daughter(dID-1);
      return tdpdg;
	// break;
    }
  }
  return 0;
}

void
WZTreeMaker::storeTriggerMatchingInfo( const edm::Event *iEvent, const Handle<pat::TriggerEvent> *patTriggerEvent, Handle<ElectronV> *electrons, Handle<MuonV> *muons)
{
  size_t vSe=hltElePaths_.size()+triggersToStore_.size();
  size_t vSm=hltMuPaths_.size()+triggersToStore_.size();
  if (iEvent==0) {
    //initialization of vector maps (for up to 4 e and mu trigger objects):
    for (int i=0;i<4;i++) {
      ostringstream ost;
      ost << i;
      storeI("trigMatchInfoEle_"+ ost.str() + "Index",-1);
      storeF("trigMatchInfoEle_"+ ost.str() + "Et",0);
      storeF("trigMatchInfoEle_"+ ost.str() + "Eta",0);
      storeF("trigMatchInfoEle_"+ ost.str() + "Phi",0);

      storeI("trigMatchInfoMu_"+ ost.str() + "Index",-1);
      storeF("trigMatchInfoMu_"+ ost.str() + "Pt",0);
      storeF("trigMatchInfoMu_"+ ost.str() + "Eta",0);
      storeF("trigMatchInfoMu_"+ ost.str() + "Phi",0);
    }
    //intMap_["trigMatchInfo_definedForRun"]= false;
    return;
  }
  for (size_t j=0;j<vSe;j++) {
    vector<int> eleMatchIndexList;
    vector<float> eleMatchEtList;
    vector<float> eleMatchEtaList;
    vector<float> eleMatchPhiList;
    string *name = &hltElePathsFull_[j];
    if (j>=hltElePaths_.size()) name = &triggersToStoreFull_[j-hltElePaths_.size()];
    for (size_t i = 0; i < (*electrons)->size(); i++) {
      const pat::Electron &p= (*electrons)->at(i);
      if (p.triggerObjectMatchesByPath(*name, true, false).size()) {
        const pat::TriggerObjectStandAlone * trigRef = p.triggerObjectMatchByPath(*name, true, false);
	eleMatchIndexList.push_back(i);
	eleMatchEtList.push_back(trigRef->et());
	eleMatchEtaList.push_back(trigRef->eta());
	eleMatchPhiList.push_back(trigRef->phi());
      }
    }
    for (size_t i=0;i<eleMatchIndexList.size();i++) {
      for (size_t c=i+1;c<eleMatchIndexList.size();c++) {
	if (eleMatchEtList[i]<eleMatchEtList[c]) {
	  int itmp=eleMatchIndexList[i];
	  eleMatchIndexList[i]=eleMatchIndexList[c];
	  eleMatchIndexList[c]=itmp;
	  float ftmp = eleMatchEtList[i];
	  eleMatchEtList[i]=eleMatchEtList[c];
	  eleMatchEtList[c]=ftmp;
	  ftmp = eleMatchEtaList[i];
	  eleMatchEtaList[i]=eleMatchEtaList[c];
	  eleMatchEtaList[c]=ftmp;
	  ftmp = eleMatchPhiList[i];
	  eleMatchPhiList[i]=eleMatchPhiList[c];
	  eleMatchPhiList[c]=ftmp;
	}
      }
    }
    for (size_t i=0;i<4;i++) {
      ostringstream ost;
      ost << i;
      bool store = i<eleMatchIndexList.size();
      storeI("trigMatchInfoEle_"+ ost.str() + "Index", store ? eleMatchIndexList[i] :-1);
      storeF("trigMatchInfoEle_"+ ost.str() + "Et",    store ? eleMatchEtList[i]    :-1.);
      storeF("trigMatchInfoEle_"+ ost.str() + "Eta",   store ? eleMatchEtaList[i]   :0.);
      storeF("trigMatchInfoEle_"+ ost.str() + "Phi",   store ? eleMatchPhiList[i]   :0.);
    }
  }
  for (size_t j=0;j<vSm;j++) {
    vector<int> muMatchIndexList;
    vector<float> muMatchPtList;
    vector<float> muMatchEtaList;
    vector<float> muMatchPhiList;
    string *name = &hltMuPathsFull_[j];
    if (j>=hltMuPaths_.size()) name = &triggersToStoreFull_[j-hltMuPaths_.size()];
    for (size_t i = 0; i < (*muons)->size(); i++) {
      const pat::Muon & p = (*muons)->at(i);
      if (p.triggerObjectMatchesByPath(*name, true, false).size()) {
        const pat::TriggerObjectStandAlone * trigRef = p.triggerObjectMatchByPath(*name, true, false);
	muMatchIndexList.push_back(i);
	muMatchPtList.push_back(trigRef->pt());
	muMatchEtaList.push_back(trigRef->eta());
	muMatchPhiList.push_back(trigRef->phi());

      }
    }
    for (size_t i=0;i<muMatchIndexList.size();i++) {
      for (size_t c=i+1;c<muMatchIndexList.size();c++) {
	if (muMatchPtList[i]<muMatchPtList[c]) {
	  int itmp=muMatchIndexList[i];
	  muMatchIndexList[i]=muMatchIndexList[c];
	  muMatchIndexList[c]=itmp;
	  float ftmp = muMatchPtList[i];
	  muMatchPtList[i]=muMatchPtList[c];
	  muMatchPtList[c]=ftmp;
	  ftmp = muMatchEtaList[i];
	  muMatchEtaList[i]=muMatchEtaList[c];
	  muMatchEtaList[c]=ftmp;
	  ftmp = muMatchPhiList[i];
	  muMatchPhiList[i]=muMatchPhiList[c];
	  muMatchPhiList[c]=ftmp;
	}
      }
    }
    //fill in values
    for (size_t i=0;i<4;i++) {
      ostringstream ost;
      ost << i;
      bool store = i< muMatchIndexList.size();
      storeI("trigMatchInfoMu_"+ ost.str() + "Index", store ? muMatchIndexList[i] : -1);
      storeF("trigMatchInfoMu_"+ ost.str() + "Pt",    store ? muMatchPtList[i]  :   -1.);
      storeF("trigMatchInfoMu_"+ ost.str() + "Eta",   store ? muMatchEtaList[i] :   0.);
      storeF("trigMatchInfoMu_"+ ost.str() + "Phi",   store ? muMatchPhiList[i] :   0.);
    }
  }
  //intMap_["trigMatchInfo_definedForRun"] = int(iEvent->id().run()) <= int(hltUpdatedForRun_);
}

void
WZTreeMaker::storeTriggerSummary(const Event * iEvent, const EventSetup * iSetup, 
                                 const TriggerResults * triggerResults, const TriggerEvent * triggerEvent)
{
  if (iEvent==0) {
      storeI("trigMatchInfoMu_pass",0);
      storeI("trigMatchInfoMu_prescale",0);
      storeI("trigMatchInfoEle_pass",0);
      storeI("trigMatchInfoEle_prescale",0);
      //format names
      for (size_t j=0;j< triggersToStore_.size(); j++){
        size_t posVer=triggersToStore_[j].rfind("_v");
        if (posVer!=string::npos) trStoreNames_.push_back( triggersToStore_[j].substr(0,posVer));
        else trStoreNames_.push_back( triggersToStore_[j].substr(0,triggersToStore_[j].rfind("*")));
        if (triggersToStore_[j].rfind("*")==string::npos) {
          trMatchExactName_.push_back(1);
          triggersToStoreAux_.push_back(triggersToStore_[j]);
        }
        else  {
          trMatchExactName_.push_back(0);
          triggersToStoreAux_.push_back(triggersToStore_[j].substr(0,triggersToStore_[j].rfind("*")));
        }
      }
      for (size_t j=0;j< hltElePaths_.size();j++) {
        size_t posVer= hltElePaths_[j].rfind("_v");
        if (posVer!=string::npos) trStoreNames_.push_back( hltElePaths_[j].substr(0,posVer));
        else trStoreNames_.push_back(hltElePaths_[j].substr(0,hltElePaths_[j].rfind("*")));
        if (hltElePaths_[j].rfind("*")==string::npos) {
          trMatchExactName_.push_back(1);
          hltElePathsAux_.push_back(hltElePaths_[j]);
        }
        else  {
          trMatchExactName_.push_back(0);
          hltElePathsAux_.push_back(hltElePaths_[j].substr(0,hltElePaths_[j].rfind("*")));
        }
      }
      for (size_t j=0;j< hltMuPaths_.size();j++) {
        size_t posVer= hltMuPaths_[j].rfind("_v");
        if (posVer!=string::npos) trStoreNames_.push_back( hltMuPaths_[j].substr(0,posVer));
        else trStoreNames_.push_back(hltMuPaths_[j].substr(0,hltMuPaths_[j].rfind("*")));
        if (hltMuPaths_[j].rfind("*")==string::npos) {
          trMatchExactName_.push_back(1);
          hltMuPathsAux_.push_back(hltMuPaths_[j]);
        }
        else  {
          trMatchExactName_.push_back(0);
          hltMuPathsAux_.push_back(hltMuPaths_[j].substr(0,hltMuPaths_[j].rfind("*")));
        }
      }
  }
  for (size_t j = 0; j < trStoreNames_.size(); j++){
      intMap_["pass_" + trStoreNames_[j]] = 0;
      intMap_["prescale_" + trStoreNames_[j]] = 0;
  }
  if (triggerResults) {
    TriggerNames names = iEvent->triggerNames(* triggerResults);
    vector<string> strings = names.triggerNames();
    //save a list of available triggers once per run
    if(!trNamesShort_.size()) {
      for (size_t i = 0; i < names.size(); i++) {
        //cout << strings.at(i) << endl;
        size_t posVer=names.triggerName(i).rfind("_v");
        string subname;
        if (posVer==string::npos) subname=names.triggerName(i);
        else subname = names.triggerName(i).substr(0,posVer+2);
        trNamesShort_.push_back(subname);
      }
      //store full names for trigger matching
      hltElePathsFull_.clear();
      hltEleIsTracker_.clear();
      for (size_t j=0;j< hltElePaths_.size();j++) {
        if (hltElePaths_[j].rfind("_v*")==string::npos)
          hltElePathsFull_.push_back(hltElePaths_[j]);
        else {
          bool foundName=false;
          for (size_t i = 0; i < names.size(); i++) {
             if (trNamesShort_[i] == hltElePathsAux_[j]) {
               hltElePathsFull_.push_back(names.triggerName(i));
               foundName=true;
               break;
             }
          }
          if (!foundName) hltElePathsFull_.push_back(hltElePaths_[j]);
        }
        //find if electron path contains tracking
        if (hltElePathsFull_[j].find("TrkId")!=string::npos) hltEleIsTracker_.push_back(0);
        else hltEleIsTracker_.push_back(1);
      }
      hltMuPathsFull_.clear();
      for (size_t j=0;j< hltMuPaths_.size();j++) {
        if (hltMuPaths_[j].rfind("_v*")==string::npos)
          hltMuPathsFull_.push_back(hltMuPaths_[j]);
        else {
          bool foundName=false;
          for (size_t i = 0; i < names.size(); i++) {
             if (trNamesShort_[i] == hltMuPathsAux_[j]) {
               hltMuPathsFull_.push_back(names.triggerName(i));
               foundName=true;
               break;
             }
          }
          if (!foundName) hltMuPathsFull_.push_back(hltMuPaths_[j]);
        }
      }
      triggersToStoreFull_.clear();
      trsEleIsTracker_.clear();
      for (size_t j=0;j< triggersToStore_.size();j++) {
        if (triggersToStore_[j].rfind("_v*")==string::npos)
          triggersToStoreFull_.push_back(triggersToStore_[j]);
        else {
          bool foundName=false;
          for (size_t i = 0; i < names.size(); i++) {
             if (trNamesShort_[i] == triggersToStoreAux_[j]) {
               triggersToStoreFull_.push_back(names.triggerName(i));
               foundName=true;
               break;
             }
          }
          if (!foundName) triggersToStoreFull_.push_back(triggersToStore_[j]);
        }
        //find if path contains trigger electron or cluster
        if (triggersToStoreFull_[j].find("TrkId")!=string::npos) trsEleIsTracker_.push_back(0);
        else trsEleIsTracker_.push_back(1);
      }
    }
    //reset trigger list
    vector<int> triggers_pass(triggersToStore_.size()+hltElePaths_.size()+hltMuPaths_.size());
    vector<int> triggers_prescale(triggersToStore_.size()+hltElePaths_.size()+hltMuPaths_.size());
    trigMatchInfoEle_pass_ = vector<int>(hltElePaths_.size());
    trigMatchInfoEle_prescale_=vector<int>(hltElePaths_.size());
    trigMatchInfoMu_pass_=vector<int>(hltMuPaths_.size());
    trigMatchInfoMu_prescale_=vector<int>(hltMuPaths_.size());

    for (size_t i = 0; i < names.size(); i++) {
      for (size_t j = 0; j < triggersToStore_.size(); j++) {
	bool pass = false;
	int prescale = 0;
	bool matchExact =(bool) trMatchExactName_[j];
	if ( (matchExact && names.triggerName(i) == triggersToStoreAux_[j])
	    || (!matchExact && trNamesShort_[i] == triggersToStoreAux_[j]) ) {
	  bool accept = triggerResults->accept(i);
	  if (accept) pass = true;
	  bool prescale_set = hltConfigProvider_.prescaleSet(*iEvent, *iSetup);
	  if (prescale_set != -1){
	    prescale = hltConfigProvider_.prescaleValue(*iEvent, *iSetup, names.triggerName(i));
	  }
          triggers_pass[j] = pass;
          triggers_prescale[j] = prescale;
	} // match name
      }
      size_t offset = triggersToStore_.size();
      for (size_t j=0;j<hltElePaths_.size();j++) {
        bool matchExact =(bool) trMatchExactName_[j+offset];
        if ( (matchExact &&  names.triggerName(i)== hltElePathsAux_[j])
            || (!matchExact && trNamesShort_[i] == hltElePathsAux_[j]) ) {
          trigMatchInfoEle_pass_[j]=triggerResults->accept(i);
          bool prescale_set = hltConfigProvider_.prescaleSet(*iEvent, *iSetup);
          if (prescale_set != -1) {
              int prescale = hltConfigProvider_.prescaleValue(*iEvent, *iSetup, names.triggerName(i));
	      trigMatchInfoEle_prescale_[j]=prescale;
	  }
          else trigMatchInfoEle_prescale_[j]=0;
	}
        triggers_pass[j+offset] = trigMatchInfoEle_pass_[j];
        triggers_prescale[j+offset] = trigMatchInfoEle_prescale_[j];
      }
      offset += hltElePaths_.size();
      for (size_t j=0;j<hltMuPaths_.size();j++) {
        bool matchExact =(bool) trMatchExactName_[j+offset];
        if ( (matchExact &&  names.triggerName(i)== hltMuPathsAux_[j])
            || (!matchExact && trNamesShort_[i] == hltMuPathsAux_[j]) ) {
          trigMatchInfoMu_pass_[j]=triggerResults->accept(i);
          bool prescale_set = hltConfigProvider_.prescaleSet(*iEvent, *iSetup);
          if (prescale_set != -1) {
              int prescale = hltConfigProvider_.prescaleValue(*iEvent, *iSetup, names.triggerName(i));
              trigMatchInfoMu_prescale_[j]=prescale;
          }
          else trigMatchInfoMu_prescale_[j]=0;
        }
        triggers_pass[j+offset] = trigMatchInfoMu_pass_[j];
        triggers_prescale[j+offset] = trigMatchInfoMu_prescale_[j];
      }
    }
    for (size_t j=0;j<triggersToStore_.size()+hltElePaths_.size()+hltMuPaths_.size();j++) {
      intMap_["pass_" + trStoreNames_[j]] = triggers_pass[j];
      intMap_["prescale_" + trStoreNames_[j]] = triggers_prescale[j];
    }
    for (size_t j=0;j<hltElePaths_.size();j++) {
      storeI("trigMatchInfoEle_pass",trigMatchInfoEle_pass_[j]);
      storeI("trigMatchInfoEle_prescale",trigMatchInfoEle_prescale_[j]);
    }
    for (size_t j=0;j<hltMuPaths_.size();j++) {
      storeI("trigMatchInfoMu_pass",trigMatchInfoMu_pass_[j]);
      storeI("trigMatchInfoMu_prescale",trigMatchInfoMu_prescale_[j]);
    }
  }
}

void
WZTreeMaker::storePileup(const std::vector<PileupSummaryInfo> * PupInfo, const Event * iEvent)
{

  floatMap_["PU_TrueNumInteractions"] = 0;
  intMap_["PU_BX_Minus1"] = 0;
  intMap_["PU_BX_Zero"] = 0;
  intMap_["PU_BX_Plus1"] = 0;
//  floatMap_["PU_Weight1BX"] = 0;
//  floatMap_["PU_Weight3BX"] = 0;
//  floatMap_["PU_Weight3D"] = 0;
  floatMap_["PU_WeightTruth"] = 0;
  
  if(PupInfo){
    float sum_nvtx = 0;
    float true_nvtx = 0;
    int nm1 = -1; int n0 = -1; int np1 = -1;
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      float npv = PVI->getPU_NumInteractions();
      sum_nvtx += float(npv);
      int BX = PVI->getBunchCrossing();
      if(BX == 0){//in time PU
        floatMap_["PU_TrueNumInteractions"] = true_nvtx = PVI->getTrueNumInteractions();
//        floatMap_["PU_Weight1BX"] = LumiWeights_.weight(npv);
      }
      if(BX == -1) intMap_["PU_BX_Minus1"] = nm1 = PVI->getPU_NumInteractions();
      if(BX == 0) intMap_["PU_BX_Zero"] = n0 = PVI->getPU_NumInteractions();
      if(BX == 1) intMap_["PU_BX_Plus1"] = np1 = PVI->getPU_NumInteractions();
    }
    float avg_nvtx = sum_nvtx/3.;//+1, 0, -1 BX
    floatMap_["PU_NumInteractions3BX"] = avg_nvtx;
//    floatMap_["PU_Weight3BX"] = LumiWeights_.weight3BX(avg_nvtx);
//    floatMap_["PU_Weight3D"] = (float) LumiWeights3D_.weight3D(nm1,n0,np1);
    //floatMap_["PU_WeightTruth"] = LumiWeightsTruth_.weight(true_nvtx);
  }
  if (iEvent) {
    Handle<double> userWeightCount_;
    iEvent->getByLabel( weightProducerName_,userWeightCount_);
    if (userWeightCount_.isValid()) {
      floatMap_["PU_WeightTruth"] = *userWeightCount_;
    }
  }


}


void 
WZTreeMaker::beginJob()
{

  if (debug_) cout << ">>>>>>>>>>>>>>>>>>> In the beginJob" << endl;

  outputFile_ = new TFile(outputFileName_.c_str(), "recreate");
  outputTree_ = new TTree("WZ", "0");

  const unsigned int nNames = numEventsNames_.size();
  numEventsHist_ = new TH1F("numEvents", "", nNames, 0, nNames);
  for (size_t i = 0; i < nNames; i++)
    numEventsHist_->GetXaxis()->SetBinLabel(i + 1, numEventsNames_[i].c_str());

  userWeightHist_ = new TH1D("userWeight", "", 2,0,2);
  userWeightHist_->GetXaxis()->SetBinLabel(1, "userWeightSum");
 
  initializeTree();

}



void 
WZTreeMaker::beginRun(const Run & run, const EventSetup & setup)
{

  bool changed;
  if(!hltConfigProvider_.init(run, setup, hltProcessName_, changed))
    cerr<<"Something went wrong with HLTConfigProvider init!!!!\n\n\n";

}



void 
WZTreeMaker::endJob() {

  if (debug_) cout << ">>>>>>>>>>>>>>>>>>> In the endJob" << endl;

  cout << "Total events processed: "
       << numEventsHist_->GetBinContent(1) << endl;

  outputFile_->cd();
  outputTree_->Write();
  numEventsHist_->Write();
  userWeightHist_->Write();

  outputFile_->Close();

}



void 
WZTreeMaker::endLuminosityBlock(const edm::LuminosityBlock & lumi, 
                               const EventSetup & setup)
{

  if (debug_) cout << ">>>>>>>>>>>>>>>>>>> In endLuminosityBlock" << endl;

  for (size_t i = 0; i < numEventsNames_.size(); i++) {

    string name = numEventsNames_[i];
    Handle<MergeableCounter> numEventsCounter;
    lumi.getByLabel(name, numEventsCounter);

    if (numEventsCounter.isValid()) {

      numEventsHist_->AddBinContent(i + 1, numEventsCounter->value);

      if (i == 0) {
        if (debug_) cout << "Adding " << numEventsCounter->value 
                         << " to total events" << endl;
      }

    }

  }

  Handle<double> userWeightCount_;
  lumi.getByLabel( weightProducerName_,userWeightCount_);
  if (userWeightCount_.isValid()) {
    userWeightHist_->AddBinContent(1,*userWeightCount_);
  }
}



void 
WZTreeMaker::endRun(const Run & run, const EventSetup & setup)
{
  trNamesShort_.clear();
  Handle<GenRunInfoProduct> genInfo;
  run.getByLabel("generator", genInfo);
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(WZTreeMaker);
