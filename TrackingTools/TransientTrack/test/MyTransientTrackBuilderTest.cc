#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include <memory>
#include <iostream>
#include <string>

using namespace edm;

class MyTransientTrackBuilderTest : public edm::one::EDAnalyzer<> {
public:
  MyTransientTrackBuilderTest(const edm::ParameterSet& pset)
      : ttkToken_(esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"}) ) {

    trkToken = consumes<reco::TrackCollection>(pset.getParameter<edm::InputTag>("TrackLabel"));
    bsToken = consumes<reco::BeamSpot>(pset.getParameter<edm::InputTag>("beamSpotLabel"));

}

  ~MyTransientTrackBuilderTest() = default;

  virtual void analyze(const edm::Event& event, const edm::EventSetup& setup) {
    using namespace std;
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  event.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }
  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  
  event.getByToken(trkToken, tks);

  // interface RECO tracks to vertex reconstruction
  const auto& theB = &setup.getData(ttkToken_);
  std::vector<reco::TransientTrack> t_tks;

    t_tks = (*theB).build(tks, beamSpot);

    edm::LogPrint("TrackerTrackBuilderTest")
        << " Asking for the TransientTrackBuilder with name TransientTrackBuilder\n";
    //const TransientTrackBuilder* theB = &setup.getData(ttkToken_);

    edm::LogPrint("TrackerTrackBuilderTest") << " Got a " << typeid(*theB).name() << endl;
    edm::LogPrint("TrackerTrackBuilderTest")
        << "Field at origin (in Testla): " << (*theB).field()->inTesla(GlobalPoint(0., 0., 0.)) << endl;
  }

private:
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;

};
