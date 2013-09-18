// -*- C++ -*-

#ifndef UserAnalysis_WZProducts_h
#define UserAnalysis_WZProducts_h

#include <vector>
#include <string>
#include <iostream>

#include <TROOT.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TFile.h>
#include <TSystem.h>

#include "DataFormats/Common/interface/MergeableCounter.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "PhysicsTools/SelectorUtils/interface/Selector.h"
// #include "FWCore/ParameterSet/interface/ProcessDesc.h"
// #include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/Framework/interface/ESHandle.h"


namespace wz {

  //// Type Definitions //////////////////////////////////////////////////////

  typedef std::vector<int> vint;
  typedef std::vector<std::string> vstring;
  typedef edm::ParameterSet PSet;
  typedef edm::MergeableCounter Counter;
  typedef math::XYZPoint Point;
  typedef math::XYZTLorentzVector LorentzVector;
  typedef std::vector<reco::Track      > TrackV;
  typedef std::vector<reco::Candidate  > CandV;
  typedef std::vector<reco::GenParticle> GenParticleV;
  typedef std::vector< edm::InputTag   > VInputTag;
  typedef std::vector< pat::Electron   > ElectronV;
  typedef std::vector< pat::Muon       > MuonV;
  typedef std::vector< pat::Jet        > JetV;
  typedef std::vector< pat::MET        > METV;



  //// Global constants //////////////////////////////////////////////////////

  const float ZMASS = 91.188;
  const float WMASS = 80.398;



  //////// Comparators (for sorting vectors) /////////////////////////////////

  struct closestToZMass {
    bool operator() (const reco::Candidate & a, const reco::Candidate & b) {
      return fabs(ZMASS - a.mass()) < fabs(ZMASS - b.mass());
    }
  };
  struct highestPt {
    bool operator() (const reco::Candidate * a, const reco::Candidate * b) {
      return a->pt() > b->pt();
    }
  };
  struct highestPtLepton {
    bool operator() (const reco::CompositeCandidate a, 
                     reco::CompositeCandidate b) {
      return a.daughter(0)->pt() > b.daughter(0)->pt();
    }
  };



  //////// Useful Algorithms /////////////////////////////////////////////////

  inline bool areIdentical(const reco::Candidate & p1, 
                           const reco::Candidate & p2)
  {
    float tolerance = 0.0001;
    if (p1.pdgId() == p2.pdgId() &&
        fabs(p1.eta() - p2.eta()) < tolerance &&
        fabs(p1.phi() - p2.phi()) < tolerance &&
        fabs(p1.pt () - p2.pt ()) < tolerance)
      return true;
    return false;
  }



  inline bool areOverlapping(const reco::Candidate & p1, 
                             const reco::Candidate & p2)
  {
    if (p1.numberOfDaughters() == 0 && p2.numberOfDaughters() == 0)
      return areIdentical(p1, p2);
    for (size_t i = 0; i < p1.numberOfDaughters(); i++)
      if (areOverlapping(p2, * p1.daughter(i)))
        return true;
    for (size_t i = 0; i < p2.numberOfDaughters(); i++)
      if (areOverlapping(p1, * p2.daughter(i)))
        return true;
    return false;
  }



  inline const reco::Candidate * findMother(const reco::Candidate * p)
  {
    const reco::Candidate * mother = p ? p->mother(0) : 0;
    if (mother) {
      if (mother->pdgId() == p->pdgId()) return findMother(mother);
      else return mother;
    }
    else return 0;
  }



  inline const reco::Candidate * findAncestor(const reco::Candidate *p)
  {
    const reco::Candidate * mother = p ? p->mother(0) : 0;
    if (!mother) return 0;
    do {
      const reco::Candidate * motherNext=mother->mother(0);
      if (motherNext) mother=motherNext;
      else break;
    }
    while (1);
    return mother;
  }


  /// Functions to get products from the event
  template <class T, class P>
    T getProduct(const P & event, std::string productName) {
    edm::Handle<T> handle;
    event.getByLabel(edm::InputTag(productName), handle);
    return * handle;
  }

  template <class T, class P>
    T getProduct(const P & event, std::string productName, T defaultVal) {
    edm::Handle<T> handle;
    event.getByLabel(edm::InputTag(productName), handle);
    if (handle.isValid()) return * handle;
    return defaultVal;
  }

  template <class T, class P>
    const T * getPointer(const P & event, std::string productName) {
    edm::Handle<T> handle;
    event.getByLabel(edm::InputTag(productName), handle);
    if (handle.isValid()) return handle.product();
    return 0;
  }

  //////// Classes ///////////////////////////////////////////////////////////

  class BosonCandidate : public reco::CompositeCandidate {
  public:
    int flavor() const {
      if (numberOfDaughters() > 0)
        return abs(daughter(0)->pdgId());
      return 0;
    }
    operator bool() const {
      return (numberOfDaughters() > 0);
    }
  protected:
    void addDaughters(const reco::Candidate & p1, const reco::Candidate & p2) {
      addDaughter(p1);
      addDaughter(p2);
      AddFourMomenta addP4;
      addP4.set(* this);
    }
    int findGenMotherId(const reco::Candidate * p) const {
      const reco::Candidate * mom = findMother(p);
      if (mom) return mom->pdgId();
      return 0;
    }
  };

  class ZCandidate : public BosonCandidate {
  public:
    ZCandidate() {genLepton1_ = genLepton2_ = 0;}
    ZCandidate(const pat::Electron & p1, const pat::Electron & p2) {
      genLepton1_ = p1.genLepton();
      genLepton2_ = p2.genLepton();
      addDaughters(p1, p2);
    }
    ZCandidate(const pat::Muon & p1, const pat::Muon & p2) {
      genLepton1_ = p1.genLepton();
      genLepton2_ = p2.genLepton();
      addDaughters(p1, p2);
    }
    const reco::Candidate * genLepton1() const {return genLepton1_;}  
    const reco::Candidate * genLepton2() const {return genLepton2_;}
    int genMotherId1() const {return findGenMotherId(genLepton1_);}
    int genMotherId2() const {return findGenMotherId(genLepton2_);}
  private:
    const reco::Candidate * genLepton1_;
    const reco::Candidate * genLepton2_;
  };
  typedef std::vector<ZCandidate > ZCandV;



  class WCandidate : public BosonCandidate {
  public:
    WCandidate() {genLepton_ = 0;}
    WCandidate(const pat::Electron & lepton, const reco::Candidate & met) {
      genLepton_ = lepton.genLepton();
      addDaughters(lepton, met);
    }
    WCandidate(const pat::Muon & lepton, const reco::Candidate & met) {
      genLepton_ = lepton.genLepton();
      addDaughters(lepton, met);
    }
    const reco::Candidate * lepton() const {return daughter(0);}
    const reco::Candidate * met() const {return daughter(1);}
    const reco::Candidate * genLepton() const {return genLepton_;}
    int genMotherId() const {return findGenMotherId(genLepton_);}
    double mt() const {
      double dphi = 1 - cos(reco::deltaPhi(lepton()->phi(), met()->phi()));
    return sqrt(2 * met()->et() * lepton()->et() * dphi);
    }
  private:
    const reco::Candidate * genLepton_;
  };
  typedef std::vector<WCandidate > WCandV;



  class WZCandidate {

  public:

    WZCandidate() {initialize_();}

    WZCandidate(const ZCandidate & Z, const WCandidate & W) {

      initialize_();
      wcand_ = & W;
      zcand_ = & Z;
      const reco::Candidate * wLep  = W.daughter(0);
      const reco::Candidate * met   = W.daughter(1);
      const reco::Candidate * zLep1 = Z.daughter(0);
      const reco::Candidate * zLep2 = Z.daughter(1);

      LorentzVector trilepP4 = zLep1->p4() + zLep2->p4() + wLep->p4();
      double term1 = sqrt(trilepP4.M2() + pow(trilepP4.pt(),2)) + met->pt();
      double term2 = (trilepP4 + met->p4()).pt();
      transMass_ = sqrt(pow(term1, 2) - pow(term2, 2));

      double dPhi         = deltaPhi(wLep->phi(), met->phi());
      double g            = (WMASS * WMASS / 2. + 
                             wLep->pt() * met->pt() * cos(dPhi));
      double a            = - pow(wLep->pt(), 2);
      double b            = 2 * g * wLep->pz();
      double c            = pow(g, 2) - pow(wLep->p(), 2) * pow(met->pt(), 2);
      double discriminant = (b * b) - (4 * a * c);
    
      if (discriminant > 0) {
    
        TVector3 leptonP3(wLep->px(), wLep->py(), wLep->pz());
    
        double   pz1 = -b/(2*a) + sqrt(discriminant)/(2*a);
        TVector3 p1  = TVector3(met->px(), met->py(), pz1);
        double   dr1 = p1.DeltaR(leptonP3);
      
        double   pz2 = -b/(2*a) - sqrt(discriminant)/(2*a);
        TVector3 p2  = TVector3(met->px(), met->py(), pz2);
        double   dr2 = p2.DeltaR(leptonP3);

        neutrinoPz_[0] = (dr1 < dr2) ? pz1 : pz2;
        neutrinoPz_[1] = (dr1 < dr2) ? pz2 : pz1;
        neutrinoPz_[2] = (fabs(pz1) > fabs(pz2)) ? pz1 : pz2;
        neutrinoPz_[3] = (fabs(pz1) > fabs(pz2)) ? pz2 : pz1;

        for (size_t i = 0; i < neutrinoPz_.size(); i++) {
          double        px = met->px();
          double        py = met->py();
          double        pz = neutrinoPz_[i];
          double        E  = sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
          LorentzVector p4 = LorentzVector(px, py, pz, E);
          invariantMass_[i] = (Z.p4() + wLep->p4() + p4).M();
        }

      }

    }

    double ht() {
      if (!zcand_ || !wcand_) return 0.;
      return (zcand_->daughter(0)->pt() + 
              zcand_->daughter(1)->pt() + 
              wcand_->lepton()->pt());
    }

    double mt() {return transMass_;}
    double transMass() {return transMass_;}

    double neutrinoPz(std::string type) {
      return neutrinoPz_[index_(type)];
    }

    double mass(std::string type) {
      return invariantMass_[index_(type)];
    }

  private:

    double transMass_;
    std::vector<double> neutrinoPz_;
    std::vector<double> invariantMass_;
    const WCandidate * wcand_;
    const ZCandidate * zcand_;

    void initialize_() {
      transMass_ = 0.;
      neutrinoPz_ = std::vector<double>(4, 0.);
      invariantMass_ = std::vector<double>(4, 0.);
      wcand_ = 0;
      zcand_ = 0;
    }

    size_t index_(std::string type) {
      if (type == "minAngle") return 0;
      if (type == "maxAngle") return 1;
      if (type == "maxPz"   ) return 2;
      if (type == "minPz"   ) return 3;
      printf("Used a non-existent index in WZCandidate\n");
      return 9;
    }

  };
  typedef std::vector<WZCandidate> WZCandV;



  //////// Selectors /////////////////////////////////////////////////////////

  /// Class derived from Selector that implements some convenience functions
  template<class T>
    class MinMaxSelector : public Selector<T> {
  public:
    bool cutOnMin(std::string param) {
      bool useMin = param.find("min") != std::string::npos;
      bool useMax = param.find("max") != std::string::npos;
      if (!useMin && !useMax) {
        std::cout << param 
                  << " has neither 'min' nor 'max' in its name!" 
                  << std::endl;
        exit(1);
      }
      return useMin;
    }
    template<class C>
      void loadFromPset(PSet params, std::string param, bool shouldSet = true) {
      C defaultValue = 0;
      if (!this->cutOnMin(param))
        defaultValue = std::numeric_limits<C>::max();
      C val = params.getUntrackedParameter<C>(param, defaultValue);
      this->push_back(param, val);
      this->set(param, shouldSet);
    }
    template<class C>
      void setPassCut(std::string param, C value, pat::strbitset & ret) {
      bool useMin = this->cutOnMin(param);
      bool passMin = useMin && value >= this->cut(param, C());
      bool passMax = !useMin && value <= this->cut(param, C());
      if (passMin || passMax || this->ignoreCut(param))
        this->passCut(ret, param);
    }
  };



  /// Selector for electrons (either barrel or endcap) based on params
  class ElectronSelectorBase : public MinMaxSelector<pat::Electron> {
  public:
    ElectronSelectorBase() {}
    ElectronSelectorBase(PSet const params) {
      // Set the last parameter to false to turn off the cut
      loadFromPset<double>(params, "minPt", true);
      loadFromPset<double>(params, "maxDeltaEta", true);
      loadFromPset<double>(params, "maxDeltaPhi", true);
      loadFromPset<double>(params, "maxSigmaEtaEta", true);
      loadFromPset<double>(params, "maxSigmaIEtaIEta", true);
      loadFromPset<double>(params, "maxCombRelIso", true);
      loadFromPset<double>(params, "minConvDist", true);
      loadFromPset<double>(params, "minConvDcot", true);
      loadFromPset<double>(params, "minEoverP", true);
      loadFromPset<double>(params, "maxHoverE", true);
      loadFromPset<double>(params, "maxCaloIso", true);
      loadFromPset<double>(params, "maxTrackIso", true);
      loadFromPset<double>(params, "maxSwissCross", true);
      loadFromPset<int>(params, "maxMissingHits", true);
    }
    virtual bool operator()(const pat::Electron & p, pat::strbitset & ret) {
      ret.set(false);
      setPassCut("minPt", p.pt(), ret);
      setPassCut("maxDeltaEta", fabs(p.deltaEtaSuperClusterTrackAtVtx()), ret);
      setPassCut("maxDeltaPhi", fabs(p.deltaPhiSuperClusterTrackAtVtx()), ret);
      setPassCut("maxSigmaEtaEta", p.scSigmaEtaEta(), ret);
      setPassCut("maxSigmaIEtaIEta", p.sigmaIetaIeta(), ret);
      setPassCut("maxCombRelIso", p.userFloat("combRelIso"), ret);
      setPassCut("minConvDist", fabs(p.convDist()), ret);
      setPassCut("minConvDcot", fabs(p.convDcot()), ret);
      setPassCut("minEoverP", p.eSuperClusterOverP(), ret);
      setPassCut("maxHoverE", p.hadronicOverEm(), ret);
      setPassCut("maxCaloIso", p.userFloat("relCaloIso"), ret);
      setPassCut("maxTrackIso", p.userFloat("relTrackIso"), ret);
      setPassCut("maxSwissCross", p.userFloat("swissCross"), ret);
      setPassCut("maxMissingHits", 
                 p.gsfTrack()->trackerExpectedHitsInner().numberOfHits(), ret);
      setIgnored(ret);
      return (bool) ret;
    }
  };



  /// A wrapper to handle the barrel/endcap split for electrons
  class ElectronSelector {
  public:
    ElectronSelector() {}
    ElectronSelector(PSet pset, std::string selectorName) {
      PSet const params = pset.getParameter<PSet>(selectorName);
      PSet paramsEB(params);
      paramsEB.augment(params.getParameter<PSet>("barrel"));
      barrelSelector_ = ElectronSelectorBase(paramsEB);
      PSet paramsEE(params);
      paramsEE.augment(params.getParameter<PSet>("endcap"));
      endcapSelector_ = ElectronSelectorBase(paramsEE);
    }
    bool operator()(const pat::Electron & p, pat::strbitset & ret) {
      if (p.isEB()) return barrelSelector_(p, ret);
      if (p.isEE()) return endcapSelector_(p, ret);
      return false;
    }
    pat::strbitset getBitTemplate() { return barrelSelector_.getBitTemplate(); }
  private:
    ElectronSelectorBase barrelSelector_;
    ElectronSelectorBase endcapSelector_;
  };



  /// Selector for muons based on params
  class MuonSelector : public MinMaxSelector<pat::Muon> {
  public:
    MuonSelector() {}
    MuonSelector(PSet pset, std::string selectorName) {
      PSet const params = pset.getParameter<PSet>(selectorName);
      // Set the last parameter to false to turn off the cut
      loadFromPset<double>(params, "minPt", true);
      loadFromPset<double>(params, "maxEta", true);
      loadFromPset<double>(params, "maxSip", true);
      loadFromPset<double>(params, "maxDxy", true);
      //loadFromPset<double>(params, "maxCombRelIso", true);
      loadFromPset<double>(params, "maxCombRelPFIso", true);
      loadFromPset<double>(params, "maxNormalizedChi2", true);
      loadFromPset<int>(params, "minIsGlobal", true);
      loadFromPset<int>(params, "minIsTracker", true);
      //loadFromPset<int>(params, "minNTrackerHits", true);
      loadFromPset<int>(params, "minNTrackerLayers", true);
      loadFromPset<int>(params, "minNPixelHits", true);
      loadFromPset<int>(params, "minNMuonHits", true);
      loadFromPset<int>(params, "minNMatches", true);
      loadFromPset<int>(params, "minNSegments", true);
    }
    virtual bool operator()(const pat::Muon & p, pat::strbitset & ret) {
      ret.set(false);
      reco::TrackRef inner = p.innerTrack();
      reco::TrackRef global = p.globalTrack();
      setPassCut("minPt", p.pt(), ret);
      setPassCut("maxEta", fabs(p.eta()), ret);
      setPassCut("maxSip", p.userFloat("sip"), ret);
      setPassCut("maxDxy", fabs(p.userFloat("d0")), ret);
      //setPassCut("maxCombRelIso", p.userFloat("combRelIso"), ret);
      //setPassCut("maxCombRelPFIso", p.userFloat("combRelPFIsoCorr"), ret);
      setPassCut("maxCombRelPFIso", p.userFloat("combRelPFIso"), ret);
      setPassCut("maxNormalizedChi2", global.isNull() ? 0. : 
                 global->normalizedChi2(), ret);
      setPassCut("minIsGlobal", p.isGlobalMuon(), ret);
      setPassCut("minIsTracker", p.isTrackerMuon(), ret);
      //setPassCut("minNTrackerHits", inner.isNull() ? 0 :
      //           inner->hitPattern().numberOfValidTrackerHits(), ret);
      setPassCut("minNTrackerLayers", inner.isNull() ? 0 :
                 inner->hitPattern().trackerLayersWithMeasurement(), ret);
      setPassCut("minNPixelHits", inner.isNull() ? 0 :
                 inner->hitPattern().numberOfValidPixelHits(), ret);
      setPassCut("minNMuonHits", global.isNull() ? 0 :
                 global->hitPattern().numberOfValidMuonHits(), ret);
      setPassCut("minNMatches", p.numberOfMatches(), ret);
      setPassCut("minNSegments", p.userInt("nSegments"), ret);
      setIgnored(ret);
      return (bool) ret;
    }
  };



  //////// Functions /////////////////////////////////////////////////////////

  /// Return a vector of non-overlapping Z candidates
  inline ZCandV getZCands(const ElectronV & electrons, 
                          const MuonV & muons, 
                          float minMass = 60., 
                          float maxMass = 120.)
  {
    ZCandV zCands;

    // Get opposite-charge pairs of electrons
    for (size_t i = 0; i < electrons.size(); i++)
      for (size_t j = i + 1; j < electrons.size(); j++)
        if (electrons[i].charge() != electrons[j].charge())
          zCands.push_back(ZCandidate(electrons[i], electrons[j]));

    // Get opposite-charge pairs of muons
    for (size_t i = 0; i < muons.size(); i++)
      for (size_t j = i + 1; j < muons.size(); j++)
        if (muons[i].charge() != muons[j].charge())
          zCands.push_back(ZCandidate(muons[i], muons[j]));
  
    // Order by difference from Z mass
    std::sort(zCands.begin(), zCands.end(), closestToZMass());

    // Get rid of overlapping candidates
    for (ZCandV::iterator i = zCands.begin(); i != zCands.end(); ++i)
      for (ZCandV::iterator j = i + 1; j != zCands.end(); ++j)
        if (areOverlapping(*i, *j)) {
          zCands.erase(j);
          j = i;
        }

    return zCands;

  }

  /// Return a WCandidate using the highest-pT non-Z lepton
  inline WCandidate getWCand(const ElectronV & electrons,
                             const MuonV & muons, 
                             const pat::MET & met,
                             const ZCandidate & zCand,
                             double minDeltaR = 0.)
  {
    if (!zCand) return WCandidate();

    WCandV wCands;

    for (ElectronV::const_iterator i = electrons.begin(); 
         i != electrons.end(); ++i)
      if (!areOverlapping(* i, * zCand.daughter(0)) &&
          !areOverlapping(* i, * zCand.daughter(1)) &&
          reco::deltaR(* i, * zCand.daughter(0)) > minDeltaR &&
          reco::deltaR(* i, * zCand.daughter(1)) > minDeltaR) {
        wCands.push_back(WCandidate(* i, met));
      }

    for (MuonV::const_iterator i = muons.begin(); 
         i != muons.end(); ++i)
      if (!areOverlapping(* i, * zCand.daughter(0)) &&
          !areOverlapping(* i, * zCand.daughter(1)) &&
          reco::deltaR(* i, * zCand.daughter(0)) > minDeltaR &&
          reco::deltaR(* i, * zCand.daughter(1)) > minDeltaR) {
        wCands.push_back(WCandidate(* i, met));
      }

    std::sort(wCands.begin(), wCands.end(), highestPtLepton());

    if (wCands.size()) return wCands[0];

    return WCandidate();
  }

  // inline const ElectronSelector getElectronSelector(const std::string type)
  // {
  //   // PythonProcessDesc builder("params.py");
  //   // edm::ParameterSet pset = builder.processDesc()->getProcessPSet()->
  //   //   getParameter<edm::ParameterSet>("FWLiteParams");
  //   // PSet const eSelectorPset = pset.getParameter<PSet>("electronSelectors");
  //   // return ElectronSelector(eSelectorPset, type);
  //   return ElectronSelector();
  // }

  inline float combIso(float rho,
                       const pat::Muon & muon,
                       const pat::Muon & m2 = pat::Muon())
  {
    float sumIso = (muon.isolationR03().sumPt +
                    muon.isolationR03().emEt +
                    muon.isolationR03().hadEt);
    if (fabs(muon.eta()) < 1.05)
      sumIso -= rho * 0.1057;
    else
      sumIso -= rho * 0.0769;
    if (m2.pt() == 0.) // m2 is a default object.
      return sumIso;
    sumIso -= m2.innerTrack()->pt();
    if (deltaR(muon.momentum(), m2.calEnergy().ecal_position) < 0.3)
      sumIso -= m2.isolationR03().emVetoEt;
    if (deltaR(muon.momentum(), m2.calEnergy().hcal_position) < 0.3)
      sumIso -= m2.isolationR03().hadVetoEt;
    return sumIso;
  }

  inline float combIso(float rho, const pat::Electron & electron)
  {
    reco::SuperClusterRef scRef = electron.superCluster();
    float R = TMath::Sqrt(scRef->x() * scRef->x() + 
                          scRef->y() * scRef->y() + 
                          scRef->z() * scRef->z());
    float Rt = TMath::Sqrt(scRef->x() * scRef->x() + 
                           scRef->y() * scRef->y());
    float sumIso = (electron.dr03TkSumPt() + 
                    electron.dr03HcalTowerSumEt() +
                    (electron.hadronicOverEm() * 
                     electron.energy() * (Rt/R)));
    if (electron.isEB()) {
      sumIso += std::max(0., electron.dr03EcalRecHitSumEt() - 1.);
      sumIso -= rho * 0.0997;
    }
    if (electron.isEE()) {
      sumIso += electron.dr03EcalRecHitSumEt();
      sumIso -= rho * 0.1123;
    }
    return sumIso;
  }

  class Algorithms {
  public:
    static const reco::Candidate * findMother(const reco::Candidate * p)
    {
      return wz::findMother(p);
    }
    static const reco::Candidate * findAncestor(const reco::Candidate *p)
    {
      return wz::findAncestor(p);
    }
    static ZCandV getZCands(const ElectronV & electrons, 
                            const MuonV & muons, 
                            float minMass = 60.,
                            float maxMass = 120.)
    {
      return wz::getZCands(electrons, muons, minMass, maxMass);
    }
    static WCandidate getWCand(const ElectronV & electrons,
                               const MuonV & muons, 
                               const pat::MET & met,
                               const ZCandidate & zCand,
                               double minDeltaR = 0.)
    {
      return wz::getWCand(electrons, muons, met, zCand, minDeltaR);
    }
    static float combIsoMuon(float rho, const pat::Muon & muon, 
                             const pat::Muon & m2 = pat::Muon())
    {
      return wz::combIso(rho, muon, m2);
    }
    static float combIsoElectron(float rho, const pat::Electron & electron)
    {
      return wz::combIso(rho, electron);
    }
    // static ElectronSelector getElectronSelector(const std::string type)
    // {
    //   return wz::getElectronSelector(type);
    // }
  };



}


#endif
