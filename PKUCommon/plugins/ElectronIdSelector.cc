
/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   Chayanit Asawatangtrakuldee chayanit@cern.ch 
 *
 * Description:
 *   - Selects "loose" and "tight" electrons needed for V-boson analysis.
 *   - Saves collection of the reference vectors of electrons passing the 
 *     required electron ID.
 * History:
 *   
 *
 *****************************************************************************/
////////////////////////////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////////////////////////////
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class ElectronIdSelector : public edm::EDProducer
{
public:
  // construction/destruction
  ElectronIdSelector(const edm::ParameterSet& iConfig);
  virtual ~ElectronIdSelector();
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();

private:  
  // member data
  edm::InputTag  src_;
  std::string    moduleLabel_;
  std::string    idLabel_;  
  bool           useDetectorIsolation_;
  bool           applyTightID_;
  bool           applyMediumID_;
  bool           applyLooseID_;
  bool           applyVetoID_;
  unsigned int nTot_;
  unsigned int nPassed_;
};



////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
ElectronIdSelector::ElectronIdSelector(const edm::ParameterSet& iConfig)
  : src_    (iConfig.getParameter<edm::InputTag>     ("src"))
  , moduleLabel_(iConfig.getParameter<std::string>   ("@module_label"))
  , idLabel_(iConfig.existsAs<std::string>("idLabel") ? iConfig.getParameter<std::string>("idLabel") : "loose")
  , useDetectorIsolation_(iConfig.existsAs<bool>("useDetectorIsolation") ? iConfig.getParameter<bool>("useDetectorIsolation") : false)
  , nTot_(0)
  , nPassed_(0)
{
  produces<std::vector<pat::Electron> >();

  /// ------- Decode the ID criteria --------
  applyTightID_ = false;
  applyMediumID_ = false;
  applyLooseID_ = false;
  applyVetoID_ = false;

  if( (idLabel_.compare("tight")==0) || 
      (idLabel_.compare("Tight")==0) || 
      (idLabel_.compare("TIGHT")==0) ||
      (idLabel_.compare("WP70")==0) ||
      (idLabel_.compare("wp70")==0) )  
    applyTightID_ = true;
  else if( (idLabel_.compare("medium")==0) ||
      (idLabel_.compare("Medium")==0) ||
      (idLabel_.compare("MEDIUM")==0) ||
      (idLabel_.compare("WP80")==0) ||
      (idLabel_.compare("wp80")==0) )  applyMediumID_ = true;
  else if( (idLabel_.compare("loose")==0) || 
      (idLabel_.compare("Loose")==0) || 
      (idLabel_.compare("LOOSE")==0) ||
      (idLabel_.compare("WP90")==0) ||
      (idLabel_.compare("wp90")==0) )  applyLooseID_ = true;
  else if( (idLabel_.compare("veto")==0) || 
      (idLabel_.compare("Veto")==0) || 
      (idLabel_.compare("VETO")==0) ||
      (idLabel_.compare("VETOid")==0) ||
      (idLabel_.compare("VetoId")==0) )  applyVetoID_ = true;
}

 
//______________________________________________________________________________
ElectronIdSelector::~ElectronIdSelector(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
 
//______________________________________________________________________________
void ElectronIdSelector::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtxs);
 
//  edm::Handle<reco::ConversionCollection> conversions;
//  iEvent.getByLabel("allConversions", conversions);

//  edm::Handle<reco::BeamSpot> beamspot_h;
//  iEvent.getByLabel("offlineBeamSpot", beamspot_h);
//  const reco::BeamSpot &beamspot = *(beamspot_h.product());

  std::auto_ptr<std::vector<pat::Electron> > passingElectrons(new std::vector<pat::Electron >);

  edm::Handle<edm::View<pat::Electron> > electrons;
  iEvent.getByLabel(src_,electrons);
  
  bool* isPassing = new bool[electrons->size()];

  double rhoVal_;
  rhoVal_=-99.;
  edm::Handle<double> rho;
  iEvent.getByLabel("fixedGridRhoFastjetAll",rho);
  rhoVal_ = *rho;

  for(unsigned int iElec=0; iElec<electrons->size(); iElec++) { 

    isPassing[iElec]=false;

    const pat::Electron& ele = electrons->at(iElec);

    // -------- Make sure that the electron is within acceptance ------
    float eta = ele.superCluster()->eta();
    bool isEB = ele.isEB() && fabs(eta) < 1.479;
    bool isEE = ele.isEE() && fabs(eta) > 1.479 && fabs(eta) < 2.5;
    //bool inAcceptance = (isEB || isEE) && (ele.ecalDrivenSeed()==1);
    float pt  = ele.pt();

    // -------- Compute Detector isolation ------
    const double PI = 4.0*atan(1.0);
    float detector_isolation = (ele.dr03TkSumPt() + 
			       std::max(0.,ele.dr03EcalRecHitSumEt()-1.0) + 
			       ele.dr03HcalTowerSumEt() - 
			       PI*0.3*0.3*rhoVal_) / pt;

//    float EffArea = ElectronEffectiveArea::GetElectronEffectiveArea( ElectronEffectiveArea::kEleGammaAndNeutralHadronIso03 , eta , ElectronEffectiveArea::kEleEAData2012);
//    float pfIso03EA = (ele.chargedHadronIso() + std::max(0.,ele.neutralHadronIso() + ele.photonIso() - EffArea*rhoVal_)) / pt;

    float isolation = 100.;
    isolation = detector_isolation;
//    if(useDetectorIsolation_) isolation = detector_isolation;
//    else isolation = pfIso03EA;

    // -------- Compute ID ------
    double sigmaIEtaIEta   = ele.sigmaIetaIeta();
    double dPhiIn    = fabs(ele.deltaPhiSuperClusterTrackAtVtx());
    double dEtaIn    = fabs(ele.deltaEtaSuperClusterTrackAtVtx());
    double hoe     = ele.hadronicOverEm();
    double ooemoop = fabs((1.0/ele.ecalEnergy() - ele.eSuperClusterOverP()/ele.ecalEnergy()));

    // impact parameter variables
    float d0vtx         = 0.0;
    float dzvtx         = 0.0;
    if (vtxs->size() > 0) {
        reco::VertexRef vtx(vtxs, 0);    
        d0vtx = ele.gsfTrack()->dxy(vtx->position());
        dzvtx = ele.gsfTrack()->dz(vtx->position());
    } else {
        d0vtx = ele.gsfTrack()->dxy();
        dzvtx = ele.gsfTrack()->dz();
    }

    // conversion rejection variables
//    bool vtxFitConversion = ConversionTools::hasMatchedConversion(ele, conversions, beamspot.position());
      bool vtxFitConversion = !(ele.passConversionVeto()==1);
//    float mHits = ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits();
    float mHits=ele.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);  

    bool isTight  = false;  /////// <--- equivalent to WP70
    bool isMedium = false;  /////// <--- equivalent to WP80
    bool isLoose  = false;  /////// <--- equivalent to WP90
    bool isVeto    = false;  /////// <--- the loosest cut for veto

    // ---------- cut-based ID -----------------

   isTight = (pt>20.) && (mHits<=1) &&
        (!vtxFitConversion) && ((isEB && isolation<0.074355 && sigmaIEtaIEta<0.010181 && dPhiIn<0.022868 && dEtaIn<0.006574 && hoe<0.037553 && ooemoop<0.131191 && fabs(d0vtx)<0.009924 && fabs(dzvtx)<0.01531) || (isEE && isolation<0.090185 && sigmaIEtaIEta<0.028766 && dPhiIn<0.032046 && dEtaIn<0.005681 && hoe<0.081902 && ooemoop<0.106055 && fabs(d0vtx)<0.027261 && fabs(dzvtx)<0.147154));

    isMedium = (pt>20.) && (mHits<=1) &&
        (!vtxFitConversion) &&
        ((isEB && isolation<0.097213 && sigmaIEtaIEta<0.010399 && dPhiIn<0.032643 && dEtaIn<0.007641 && hoe<0.060662 && ooemoop<0.153897 && fabs(d0vtx)<0.011811 && fabs(dzvtx)<0.070775) ||
         (isEE && isolation<0.116708 && sigmaIEtaIEta<0.029524 && dPhiIn<0.042447 && dEtaIn<0.009285 && hoe<0.104263 && ooemoop<0.137468 && fabs(d0vtx)<0.051682 && fabs(dzvtx)<0.180720));

    isLoose = (pt>20.) && (mHits<=1) &&
         (!vtxFitConversion) &&
        ((isEB && isolation<0.120026 && sigmaIEtaIEta<0.120026 && dPhiIn<0.072624  && dEtaIn<0.012442  && hoe<0.121476 && ooemoop<0.221803 && fabs(d0vtx)<0.022664 && fabs(dzvtx)<0.173670) ||
         (isEE && isolation<0.162914 && sigmaIEtaIEta<0.032602  && dPhiIn<0.145129 && dEtaIn<0.010654 && hoe<0.131862 && ooemoop<0.142283 && fabs(d0vtx)<0.097358 && fabs(dzvtx)<0.198444));

    isVeto = (pt>10.) && (mHits<=999) &&
         (!vtxFitConversion) &&
        ((isEB && isolation<0.164369 && sigmaIEtaIEta<0.011100 && dPhiIn<0.252044 && dEtaIn<0.016315 && hoe<0.345843 && ooemoop<0.248070  && fabs(d0vtx)<0.060279 && fabs(dzvtx)<0.800538) ||
         (isEE && isolation<0.212604 && sigmaIEtaIEta<0.033987 && dPhiIn<0.245263 && dEtaIn<0.010671 && hoe<0.134691 && ooemoop<0.157160  && fabs(d0vtx)<0.273097 && fabs(dzvtx)<0.885860));

    /// ------- Finally apply selection --------
    if(applyTightID_ && isTight)   isPassing[iElec]= true;
    if(applyMediumID_ && isMedium) isPassing[iElec]= true;
    if(applyLooseID_ && isLoose)   isPassing[iElec]= true;
    if(applyVetoID_ && isVeto) isPassing[iElec]= true;
    
 }
   
  unsigned int counter=0;
  edm::View<pat::Electron>::const_iterator tIt, endcands = electrons->end();
  for (tIt = electrons->begin(); tIt != endcands; ++tIt, ++counter) {
    if(isPassing[counter]) passingElectrons->push_back( *tIt );  
  }

  nTot_  +=electrons->size();
  nPassed_+=passingElectrons->size();

  delete [] isPassing;  
  iEvent.put(passingElectrons);
}

 
//______________________________________________________________________________
void ElectronIdSelector::endJob()
{
  std::stringstream ss;
  ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    <<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   <<"\n"<<moduleLabel_<<"(ElectronIdSelector) SUMMARY:\n"<<ss.str()
	   <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef ElectronIdSelector   			    PATElectronIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATElectronIdSelector);
