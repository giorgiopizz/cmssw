// -*- C++ -*-
//
// Package:    TrackingTools/TTProdNew
// Class:      TTProdNew
//
/**\class TTProdNew TTProdNew.cc TrackingTools/TTProdNew/plugins/TTProdNew.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Giorgio Pizzati
//         Created:  Thu, 03 Mar 2022 07:51:36 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
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


//
// class declaration
//

class TTProdNew : public edm::stream::EDProducer<> {
public:
  explicit TTProdNew(const edm::ParameterSet&);
  ~TTProdNew();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;
};

//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
TTProdNew::TTProdNew(const edm::ParameterSet& iConfig) : ttkToken_(esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"}) )
{
    trkToken = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackLabel"));
    bsToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotLabel"));
  //register your products
  //produces<std::vector<reco::TransientTrack>>("TransientTrack");
/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/
  //now do what ever other initialization is needed
}

TTProdNew::~TTProdNew() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void TTProdNew::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from iEventiSetup";
  }
  // get RECO tracks from the iEvent
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  
  iEvent.getByToken(trkToken, tks);

  // interface RECO tracks to vertex reconstruction
  const auto& theB = &iSetup.getData(ttkToken_);
  std::vector<reco::TransientTrack> t_tks;

  t_tks = (*theB).build(tks, beamSpot);
    edm::LogPrint("TrackerTrackBuilderTest")
        << " Asking for the TransientTrackBuilder with name TransientTrackBuilder\n";
    //const TransientTrackBuilder* theB = &iSetup.getData(ttkToken_);

    edm::LogPrint("TrackerTrackBuilderTest") << " Got a " << typeid(*theB).name() << endl;
    edm::LogPrint("TrackerTrackBuilderTest")
        << "Field at origin (in Testla): " << (*theB).field()->inTesla(GlobalPoint(0., 0., 0.)) << endl;
   // do something with TransientTracks: call kernel ecc...
   
/* NO MORE NEEDED
  // insert transient track in event
  //
  auto result = std::make_unique<std::vector<reco::TransientTrack>>();
  result->push_back(t_tks.at(0));
  iEvent.put(std::move(result), "TransientTrack");
*/

/* This is an event example
  //Read 'ExampleData' from the Event
  ExampleData const& in = iEvent.get(inToken_);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  iEvent.put(std::make_unique<ExampleData2>(in));
*/

/* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
*/
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void TTProdNew::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void TTProdNew::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
TTProdNew::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
TTProdNew::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
TTProdNew::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
TTProdNew::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TTProdNew::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TTProdNew);
