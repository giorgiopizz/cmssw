// -*- C++ -*-
//
// Class:      gpuTest
//
/**\class gpuTest gpuTest.cc TrackingTools/gpuTest/plugins/gpuTest.cc

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

#include <cuda_runtime.h>

#include "CUDADataFormats/Common/interface/Product.h"
#include "CUDADataFormats/Common/interface/HostProduct.h"
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

#include <cuda_runtime.h>

#include "CUDADataFormats/Common/interface/Product.h"

#include "RecoVertex/PrimaryVertexProducer/test/gpuKernel.h"
#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"

//
// class declaration
//

class gpuTest : public edm::stream::EDProducer<edm::ExternalWork> {
public:
  explicit gpuTest(const edm::ParameterSet&);
  ~gpuTest();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void acquire(edm::Event&, const edm::EventSetup&,edm::WaitingTaskWithArenaHolder waitingTaskHolder);
  void produce(edm::Event&, const edm::EventSetup&) override;

private:
//  void beginStream(edm::StreamID) override;
//  void produce(edm::Event&, const edm::EventSetup&) override;
//  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  
  const gpuKernel::Producer gpuAlgo;

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttkToken_;
  edm::EDGetTokenT<reco::BeamSpot> bsToken;
  edm::EDGetTokenT<reco::TrackCollection> trkToken;

//  edm::EDPutTokenT<cms::cuda::Product<reco::VertexCollection>>();
  cms::cuda::ContextState ctxState_;
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
gpuTest::gpuTest(const edm::ParameterSet& iConfig) : ttkToken_(esConsumes(edm::ESInputTag{"", "TransientTrackBuilder"}) )
{
    trkToken = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TrackLabel"));
    bsToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotLabel"));
  //register your products
//    produces<cms::cuda::Product<reco::VertexCollection>>();

/* Examples
  produces<ExampleData2>();

  //if do put with a label
  produces<ExampleData2>("label");
 
  //if you want to put into the Run
  produces<ExampleData2,InRun>();
*/
  //now do what ever other initialization is needed
}

gpuTest::~gpuTest() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

void gpuTest::acquire(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::WaitingTaskWithArenaHolder waitingTaskHolder) {
    cms::cuda::ScopedContextAcquire ctx{iEvent.streamID(), std::move(waitingTaskHolder), ctxState_};

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
    
  gpuAlgo.makeAsync(ctx.stream(), t_tks);

}

// ------------ method called to produce the data  ------------
void gpuTest::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  cms::cuda::ScopedContextProduce ctx{ctxState_};
    
  // now we should get results from gpuAlgo  

}

// ------------ method called when starting to processes a run  ------------
/*
void
gpuTest::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
// ------------ method called when ending the processing of a run  ------------
/*
void
gpuTest::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
gpuTest::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
gpuTest::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void gpuTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(gpuTest);
