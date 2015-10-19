
/*****************************************************************************
 * Project: CMS detector at the CERN
 *
 * Package: ElectroWeakAnalysis/VPlusJets
 *
 *
 * Authors:
 *
 *   PKU 
 *
 * Description:
 *   - Selects "loose" and "tight" muons needed for V-boson analysis.
 *   - Saves collection of the reference vectors of muons passing the 
 *     required muon ID.
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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include <memory>
#include <vector>
#include <sstream>
#include <cmath>

////////////////////////////////////////////////////////////////////////////////
// class definition
////////////////////////////////////////////////////////////////////////////////
class MuonIdSelector : public edm::EDProducer
{
public:
  // construction/destruction
  MuonIdSelector(const edm::ParameterSet& iConfig);
  virtual ~MuonIdSelector();
  
  // member functions
  void produce(edm::Event& iEvent,const edm::EventSetup& iSetup);
  void endJob();

private:  
  // member data
  edm::InputTag  src_;
  std::string    moduleLabel_;
  std::string    idLabel_;  
  bool           applyTightID_;
  bool           applyLooseID_;

  unsigned int nTot_;
  unsigned int nPassed_;
};



////////////////////////////////////////////////////////////////////////////////
// construction/destruction
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
MuonIdSelector::MuonIdSelector(const edm::ParameterSet& iConfig)
  : src_    (iConfig.getParameter<edm::InputTag>     ("src"))
  , moduleLabel_(iConfig.getParameter<std::string>   ("@module_label"))
  , idLabel_(iConfig.existsAs<std::string>("idLabel") ? iConfig.getParameter<std::string>("idLabel") : "loose")
  , nTot_(0)
  , nPassed_(0)
{
  produces<std::vector<pat::Muon> >();

  /// ------- Decode the ID criteria --------
  applyTightID_ = false;
  applyLooseID_ = false;

  if( (idLabel_.compare("tight")==0) || 
      (idLabel_.compare("Tight")==0) || 
      (idLabel_.compare("TIGHT")==0) ||
      (idLabel_.compare("WP70")==0) ||
      (idLabel_.compare("wp70")==0) )  
    applyTightID_ = true;
  else if( (idLabel_.compare("loose")==0) || 
      (idLabel_.compare("Loose")==0) || 
      (idLabel_.compare("LOOSE")==0) ||
      (idLabel_.compare("WP90")==0) ||
      (idLabel_.compare("wp90")==0) )  applyLooseID_ = true;
}

 
//______________________________________________________________________________
MuonIdSelector::~MuonIdSelector(){}


////////////////////////////////////////////////////////////////////////////////
// implementation of member functions
////////////////////////////////////////////////////////////////////////////////
 
//______________________________________________________________________________
void MuonIdSelector::produce(edm::Event& iEvent,const edm::EventSetup& iSetup)
{

  /////// Pileup density "rho" in the event from fastJet pileup calculation /////
 
  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByLabel("offlineSlimmedPrimaryVertices", vtxs);
 
  std::auto_ptr<std::vector<pat::Muon> > passingMuons(new std::vector<pat::Muon >);

  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByLabel(src_,muons);
  
  bool* isPassing = new bool[muons->size()];

  for(unsigned int iMu=0; iMu<muons->size(); iMu++) { 

    isPassing[iMu]=false;

    const pat::Muon& mu1 = muons->at(iMu);

  
 
 
    float isolation = 100.;
    isolation =  (mu1.pfIsolationR04().sumChargedHadronPt+ std::max(0.0,mu1.pfIsolationR04().sumNeutralHadronEt+mu1.pfIsolationR04().sumPhotonEt-0.5*mu1.pfIsolationR04().sumPUPt))/mu1.pt();;
 
  
    // impact parameter variables
    float d0vtx         = 0.0;
    float dzvtx         = 0.0;
    if (vtxs->size() > 0) {
        reco::VertexRef vtx(vtxs, 0);    
        d0vtx = mu1.muonBestTrack()->dxy(vtx->position());
        dzvtx = mu1.muonBestTrack()->dz(vtx->position());
    } else {
        d0vtx = mu1.muonBestTrack()->dxy();
        dzvtx = mu1.muonBestTrack()->dz();
    }
 
 

    bool isTight  = false;  /////// <--- equivalent to WP70
    bool isLoose  = false;  /////// <--- equivalent to WP90


  if(mu1.pt()>20 && abs(mu1.eta())<2.1 && mu1.isGlobalMuon() && mu1.isPFMuon() && (mu1.globalTrack()->normalizedChi2())<10 && (mu1.globalTrack()->hitPattern().numberOfValidMuonHits())>0 && (mu1.numberOfMatchedStations())>1 && abs(d0vtx)<0.2 && abs(dzvtx)<0.5 && (mu1.numberOfMatchedStations())>1 && (mu1.globalTrack()->hitPattern().trackerLayersWithMeasurement())>5 && abs(isolation)<0.12) { isTight = true;}
 

  if(mu1.pt()>20 && abs(mu1.eta())<2.1  && (mu1.isGlobalMuon() || mu1.isTrackerMuon()) && mu1.isPFMuon() && abs(isolation)<0.20) { isLoose = true;}

    /// ------- Finally apply selection --------
    if(applyTightID_ && isTight)   isPassing[iMu]= true;
    if(applyLooseID_ && isLoose)   isPassing[iMu]= true;
 
    
 }
   
  unsigned int counter=0;
  edm::View<pat::Muon>::const_iterator tIt, endcands = muons->end();
  for (tIt = muons->begin(); tIt != endcands; ++tIt, ++counter) {
    if(isPassing[counter]) passingMuons->push_back( *tIt );  
  }

  nTot_  +=muons->size();
  nPassed_+=passingMuons->size();

  delete [] isPassing;  
  iEvent.put(passingMuons);
}

 
//______________________________________________________________________________
void MuonIdSelector::endJob()
{
  std::stringstream ss;
  ss<<"nTot="<<nTot_<<" nPassed="<<nPassed_
    <<" effPassed="<<100.*(nPassed_/(double)nTot_)<<"%\n";
  std::cout<<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   <<"\n"<<moduleLabel_<<"(MuonIdSelector) SUMMARY:\n"<<ss.str()
	   <<"++++++++++++++++++++++++++++++++++++++++++++++++++"
	   << std::endl;
}


////////////////////////////////////////////////////////////////////////////////
// plugin definition
////////////////////////////////////////////////////////////////////////////////
typedef MuonIdSelector    PATMuonIdSelector;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PATMuonIdSelector);
