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
//#include <memory>

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



#include "HeterogeneousCore/CUDACore/interface/ScopedContext.h"

#include "gpuKernel.cu"

//
// class declaration
//

class gpuTest : public edm::stream::EDProducer<edm::ExternalWork> {
public:
  explicit gpuTest(const edm::ParameterSet&);
  ~gpuTest() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, edm::WaitingTaskWithArenaHolder waitingTaskHolder) override;
  void produce(edm::Event& iEvent,  edm::EventSetup const& iSetup) override;

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

//
// member functions
//


/*
struct track_t {
    double z;                        // z-coordinate at point of closest approach to the beamline
    double dz2;                      // square of the error of z(pca)
    const reco::TransientTrack *tt;  // a pointer to the Transient Track
    double Z;                        // Z[i] for DA clustering
    double pi;                       // track weight
};

*/
void gpuTest::acquire(edm::Event const& iEvent, edm::EventSetup const& iSetup, edm::WaitingTaskWithArenaHolder waitingTaskHolder) {
  cms::cuda::ScopedContextAcquire ctx{iEvent.streamID(), std::move(waitingTaskHolder), ctxState_};

//  using namespace edm;
//  using namespace std;
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
/*
  edm::LogPrint("TrackerTrackBuilderTest")
      << " Asking for the TransientTrackBuilder with name TransientTrackBuilder\n";
  //const TransientTrackBuilder* theB = &iSetup.getData(ttkToken_);

  edm::LogPrint("TrackerTrackBuilderTest") << " Got a " << typeid(*theB).name() << std::endl;
  edm::LogPrint("TrackerTrackBuilderTest")
      << "Field at origin (in Testla): " << (*theB).field()->inTesla(GlobalPoint(0., 0., 0.)) << std::endl;
*/

  gpuKernel::track_SoA tks_SoA;
  double vertexSize_ = 0.006;
  double d0CutOff_ = 3;
  int i_t = 0;
  for (std::vector<reco::TransientTrack>::const_iterator it = t_tks.begin(); it!= t_tks.end(); it++){
         
    //gpuKernel::track_t t;
    tks_SoA.z[i_t] = ((*it).stateAtBeamLine().trackStateAtPCA()).position().z();
    double tantheta = tan(((*it).stateAtBeamLine().trackStateAtPCA()).momentum().theta());
    double phi = ((*it).stateAtBeamLine().trackStateAtPCA()).momentum().phi();
    //  get the beam-spot
    reco::BeamSpot beamspot = (it->stateAtBeamLine()).beamSpot();
    tks_SoA.dz2[i_t] = pow((*it).track().dzError(), 2)  // track errror
            + (pow(beamspot.BeamWidthX() * cos(phi), 2) + pow(beamspot.BeamWidthY() * sin(phi), 2)) /
                  pow(tantheta, 2)  // beam-width induced
            + pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
    if (d0CutOff_ > 0) {
      Measurement1D IP = (*it).stateAtBeamLine().transverseImpactParameter();       // error constains beamspot
      tks_SoA.pi[i_t] = 1. / (1. + exp(pow(IP.value() / IP.error(), 2) - pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
    } else {
      tks_SoA.pi[i_t] = 1.;
    }
    //t.tt = &(*it);
    tks_SoA.Z[i_t] = 1.;
  }
 
  gpuAlgo.makeAsync(ctx.stream(), tks_SoA);

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
