#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducerDumbFitter.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "RecoVertex/VertexTools/interface/GeometricAnnealing.h"

#include "RecoVertex/LinearizationPointFinders/interface/CrossingPtBasedLinearizationPointFinder.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include <cmath>
#include <algorithm>

namespace {
  inline GlobalPoint operator-(const GlobalPoint& a, const GlobalPoint& b) {
    return GlobalPoint(a.x() - b.x(), a.y() - b.y(), a.z() - b.z());
  }

  inline GlobalPoint operator+(const GlobalPoint& a, const GlobalPoint& b) {
    return GlobalPoint(a.x() + b.x(), a.y() + b.y(), a.z() + b.z());
  }

  inline GlobalPoint operator/(const GlobalPoint& a, const double b) {
    return GlobalPoint(a.x() / b, a.y() / b, a.z() / b);
  }

#ifndef __clang__
  inline GlobalPoint operator*(const GlobalPoint& a, const double b) {
    return GlobalPoint(a.x() * b, a.y() * b, a.z() * b);
  }

  inline GlobalPoint operator*(const double b, const GlobalPoint& a) {
    return GlobalPoint(a.x() * b, a.y() * b, a.z() * b);
  }
#endif
}


struct CompareTwoTracks {
    int operator()(const reco::TransientTrack &a, const reco::TransientTrack &b) {
      return a.initialFreeState().momentum().mag() > b.initialFreeState().momentum().mag();
      //       return a.p() > b.p();
    };
};

PrimaryVertexProducerDumbFitter::PrimaryVertexProducerDumbFitter(const edm::ParameterSet& conf)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  fVerbose = conf.getUntrackedParameter<bool>("verbose", false);

  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  f4D = false;

  // select and configure the track selection
  std::string trackSelectionAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
  if (trackSelectionAlgorithm == "filter") {
    theTrackFilter = new TrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else if (trackSelectionAlgorithm == "filterWithThreshold") {
    theTrackFilter = new HITrackFilterForPVFinding(conf.getParameter<edm::ParameterSet>("TkFilterParameters"));
  } else {
    throw VertexException("PrimaryVertexProducer: unknown track selection algorithm: " + trackSelectionAlgorithm);
  }

  // select and configure the track clusterizer
  std::string clusteringAlgorithm =
      conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<std::string>("algorithm");
  if (clusteringAlgorithm == "gap") {
    theTrackClusterizer = new GapClusterizerInZ(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkGapClusParameters"));
  } else if (clusteringAlgorithm == "DA") {
    theTrackClusterizer = new DAClusterizerInZ(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  }
  // provide the vectorized version of the clusterizer, if supported by the build
  else if (clusteringAlgorithm == "DA_vect") {
    theTrackClusterizer = new DAClusterizerInZ_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
  } else if (clusteringAlgorithm == "DA2D_vect") {
    theTrackClusterizer = new DAClusterizerInZT_vect(
        conf.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
    f4D = true;
  }

  else {
    throw VertexException("PrimaryVertexProducer: unknown clustering algorithm: " + clusteringAlgorithm);
  }

  if (f4D) {
    trkTimesToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes<edm::ValueMap<float> >(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
  }

  // select and configure the vertex fitters
  if (conf.exists("vertexCollections")) {
    std::vector<edm::ParameterSet> vertexCollections =
        conf.getParameter<std::vector<edm::ParameterSet> >("vertexCollections");

    for (std::vector<edm::ParameterSet>::const_iterator algoconf = vertexCollections.begin();
         algoconf != vertexCollections.end();
         algoconf++) {
      algo algorithm;
      std::string fitterAlgorithm = algoconf->getParameter<std::string>("algorithm");
      if (fitterAlgorithm == "KalmanVertexFitter") {
        algorithm.fitter = new KalmanVertexFitter();
      } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
        algorithm.fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
      } else {
        throw VertexException("PrimaryVertexProducer: unknown algorithm: " + fitterAlgorithm);
      }
      algorithm.label = algoconf->getParameter<std::string>("label");
      algorithm.minNdof = algoconf->getParameter<double>("minNdof");
      algorithm.useBeamConstraint = algoconf->getParameter<bool>("useBeamConstraint");
      algorithm.vertexSelector =
          new VertexCompatibleWithBeam(VertexDistanceXY(), algoconf->getParameter<double>("maxDistanceToBeam"));
      algorithms.push_back(algorithm);

      produces<reco::VertexCollection>(algorithm.label);
    }
  } else {
    edm::LogWarning("MisConfiguration")
        << "this module's configuration has changed, please update to have a vertexCollections=cms.VPSet parameter.";

    algo algorithm;
    std::string fitterAlgorithm = conf.getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.fitter = new KalmanVertexFitter();
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      algorithm.fitter = new AdaptiveVertexFitter();
    } else {
      throw VertexException("PrimaryVertexProducerAlgorithm: unknown algorithm: " + fitterAlgorithm);
    }
    algorithm.label = "";
    algorithm.minNdof = conf.getParameter<double>("minNdof");
    algorithm.useBeamConstraint = conf.getParameter<bool>("useBeamConstraint");

    algorithm.vertexSelector = new VertexCompatibleWithBeam(
        VertexDistanceXY(),
        conf.getParameter<edm::ParameterSet>("PVSelParameters").getParameter<double>("maxDistanceToBeam"));

    algorithms.push_back(algorithm);
    produces<reco::VertexCollection>(algorithm.label);
  }

  //check if this is a recovery iteration
  fRecoveryIteration = conf.getParameter<bool>("isRecoveryIteration");
  if (fRecoveryIteration) {
    if (algorithms.empty()) {
      throw VertexException("PrimaryVertexProducer: No algorithm specified. ");
    } else if (algorithms.size() > 1) {
      throw VertexException(
          "PrimaryVertexProducer: Running in Recovery mode and more than one algorithm specified.  Please "
          "only one algorithm.");
    }
    recoveryVtxToken = consumes<reco::VertexCollection>(conf.getParameter<edm::InputTag>("recoveryVtxCollection"));
  }
}
TransientVertex weightedMean(const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0, chi2=0;
     unsigned int ndof = 0;
     //bool isValid = false;
     for (const auto& p : points){ 
     //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

            wxy = pow(p.second.x(),2);
            wxy = wxy == 0. ? 1. : 1. / wxy;
            x += p.first.x() * wxy;
            y += p.first.y() * wxy;
            s_wxy += wxy;

            wz = pow(p.second.z(),2);
            wz = wz == 0. ? 1. : 1. / wz;
            z += p.first.z() * wz;
            s_wz += wz;
     }

     //auto err = v.positionError().matrix();
     AlgebraicSymMatrix33 err;
     if ( s_wxy != 0. && s_wz != 0. ){
         x /= s_wxy;     
         y /= s_wxy;     
         z /= s_wz;     
         float err_x=0, err_y=0, err_z=0;
         for (const auto& p : points){ 
            wxy = pow(p.second.x(),2);
            wxy = wxy == 0. ? 1. : 1. / wxy;

            wz = pow(p.second.z(),2);
            wz = wz == 0. ? 1. : 1. / wz;

            err_x += wxy * pow(p.first.x() - x, 2);
            err_y += wxy * pow(p.first.y() - y, 2);
            err_z += wz * pow(p.first.z() - z, 2);
         }
         err(0,0) = std::sqrt( err_x / (s_wxy * (points.size() -1) ));
         err(1,1) = std::sqrt( err_y / (s_wxy * (points.size() -1) ));
         err(2,2) = std::sqrt( err_z / (s_wz *  (points.size() -1) ));
         //isValid = true;
         float dist = 0; 
         //float varikk = 0;
         //float var_x = 0, var_y = 0, var_z = 0;
         //std::vector<float> weights;
         for (const auto& p : points){ 
              //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            wxy = p.second.x();
            wz = p.second.z();
            dist = std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
            dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
            dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
            //weights.push_back(dist)

            /*
            dist = std::pow(p.first.x() - x, 2) / p.first.x();
            dist += std::pow(p.first.y() - y, 2) / p.first.y();
            dist += std::pow(p.first.z() - z, 2) / p.first.z(); 
            */
            chi2 += dist;
            ndof ++;
            
           //var_x += wxy * pow(p.first.x() - x, 2)
         }
         TransientVertex v(GlobalPoint(x,y,z), err, *iclus, chi2, ndof);
         return v;
    }
    err(0,0) = 10000;
    err(1,1) = 10000;
    err(2,2) = 10000;
    return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
}

TransientVertex weightedMean1Iter(const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0, chi2=0;
     unsigned int ndof = 0;
     //bool isValid = false;
     for (const auto& p : points){ 
            //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

            wxy = pow(p.second.x(),2);
            wxy = wxy == 0. ? 1. : 1. / wxy;
            x += p.first.x() * wxy;
            y += p.first.y() * wxy;
            s_wxy += wxy;

            wz = pow(p.second.z(),2);
            wz = wz == 0. ? 1. : 1. / wz;
            z += p.first.z() * wz;
            s_wz += wz;
     }

     //auto err = v.positionError().matrix();
     AlgebraicSymMatrix33 err;
     if ( s_wxy != 0. && s_wz != 0. ){
         x /= s_wxy;     
         y /= s_wxy;     
         z /= s_wz;     
         float err_x=0, err_y=0, err_z=0;
         for (const auto& p : points){ 
            wxy = pow(p.second.x(),2);
            wxy = wxy == 0. ? 1. : 1. / wxy;

            wz = pow(p.second.z(),2);
            wz = wz == 0. ? 1. : 1. / wz;

            err_x += wxy * pow(p.first.x() - x, 2);
            err_y += wxy * pow(p.first.y() - y, 2);
            err_z += wz  * pow(p.first.z() - z, 2);
         }
//         err(0,0) = std::sqrt( err_x / (s_wxy * (points.size() -1) ));
//         err(1,1) = std::sqrt( err_y / (s_wxy * (points.size() -1) ));
//         err(2,2) = std::sqrt( err_z / (s_wz *  (points.size() -1) ));
         //isValid = true;
         float dist = 0; 
         //float varikk = 0;
         //float var_x = 0, var_y = 0, var_z = 0;
         //std::vector<float> weights;
         s_wxy = 0; s_wz = 0;
         float old_x = x, old_y = y, old_z = z;
         x=0; y=0; z=0;
         for (const auto& p : points){ 
              //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();

//            wxy = p.second.x();
//            wz = p.second.z();
//            dist = std::pow(p.first.x() - old_x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
//            dist += std::pow(p.first.y() - old_y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
//            dist += std::pow(p.first.z() - old_z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
//
            dist =  std::pow(p.first.x() - old_x, 2); 
            dist += std::pow(p.first.y() - old_y, 2); 
            dist += std::pow(p.first.z() - old_z, 2);
            //weights.push_back(dist)
            dist = dist == 0 ? 10000 : 1. / dist;

//            wxy = pow(p.second.x(),2);
//            wxy = wxy == 0. ? 1. : 1. / wxy;
//
//            wz = pow(p.second.z(),2);
//            wz = wz == 0. ? 1. : 1. / wz;

            x += p.first.x() * dist;
            y += p.first.y() * dist;
            z += p.first.z() * dist;

            s_wxy += dist ;
            s_wz  += dist ;

            /*
            dist = std::pow(p.first.x() - x, 2) / p.first.x();
            dist += std::pow(p.first.y() - y, 2) / p.first.y();
            dist += std::pow(p.first.z() - z, 2) / p.first.z(); 
            */
            //chi2 += dist;
            //ndof ++;
            
           //var_x += wxy * pow(p.first.x() - x, 2)
         }
         if ( s_wxy != 0. && s_wz != 0. ){
             x /= s_wxy;     
             y /= s_wxy;     
             z /= s_wz;     
             err_x=0; err_y=0; err_z=0;
             //for (const auto& p : points){ 
             for (int i=0; i< (int)points.size(); i++){ 
                std::pair<GlobalPoint,GlobalPoint> p = points[i];
                TrajectoryStateClosestToPoint tscp = (*iclus)[i].trajectoryStateClosestToPoint(GlobalPoint(x,y,z));
                GlobalPoint trackPos = tscp.theState().position();
                dist = std::pow(p.first.x() - old_x, 2); 
                dist += std::pow(p.first.y() - old_y, 2);
                dist += std::pow(p.first.z() - old_z, 2);
                dist = dist == 0 ? 10000 : 1. / dist;

//                wxy = pow(p.second.x(),2);
//                wxy = wxy == 0. ? 1. : 1. / wxy;
//
//                wz = pow(p.second.z(),2);
//                wz = wz == 0. ? 1. : 1. / wz;

//                err_x += dist * pow(p.first.x() - x, 2);
//                err_y += dist * pow(p.first.y() - y, 2);
//                err_z += dist * pow(p.first.z() - z, 2);
                err_x += dist * pow(trackPos.x() - x, 2);
                err_y += dist * pow(trackPos.y() - y, 2);
                err_z += dist * pow(trackPos.z() - z, 2);
             }
             err(0,0) = std::sqrt( err_x / (s_wxy * (points.size() -1) )) / 1.;
             err(1,1) = std::sqrt( err_y / (s_wxy * (points.size() -1) )) / 1.;
             err(2,2) = std::sqrt( err_z / (s_wz *  (points.size() -1) )) / 1.;
             //isValid = true;
             for (const auto& p : points){
                  //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
                //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                wxy = p.second.x();
                wz = p.second.z();
                dist = std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
                dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
                dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
                //weights.push_back(dist)
                /*
                x += p.first.x() * dist * p.second.x();
                y += p.first.y() * dist * p.second.x();
                z += p.first.z() * dist * p.second.z();
                s_wxy += dist * p.second.x();
                s_wz +=  dist * p.second.z();
                */
                /*
                dist = std::pow(p.first.x() - x, 2) / p.first.x();
                dist += std::pow(p.first.y() - y, 2) / p.first.y();
                dist += std::pow(p.first.z() - z, 2) / p.first.z(); 
                */
                chi2 += dist;
                ndof ++;
                
               //var_x += wxy * pow(p.first.x() - x, 2)
             }
            TransientVertex v(GlobalPoint(x,y,z), err, *iclus, chi2, ndof);
            return v;
        }
    }
    err(0,0) = 10000;
    err(1,1) = 10000;
    err(2,2) = 10000;
    return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
}

TransientVertex weightedMeanManyIter(const std::vector<std::pair<GlobalPoint, GlobalPoint>>& points, std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0, chi2=0;
     unsigned int ndof = 0;
     //bool isValid = false;
     for (const auto& p : points){ 
            //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
            //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

            wxy = pow(p.second.x(),2);
            wxy = wxy == 0. ? 1. : 1. / wxy;
            x += p.first.x() * wxy;
            y += p.first.y() * wxy;
            s_wxy += wxy;

            wz = pow(p.second.z(),2);
            wz = wz == 0. ? 1. : 1. / wz;
            z += p.first.z() * wz;
            s_wz += wz;
     }

     //auto err = v.positionError().matrix();
     AlgebraicSymMatrix33 err;
     if ( s_wxy != 0. && s_wz != 0. ){
         x /= s_wxy;     
         y /= s_wxy;     
         z /= s_wz;     
         float err_x=0, err_y=0, err_z=0;
         for (const auto& p : points){ 
            wxy = pow(p.second.x(),2);
            wxy = wxy == 0. ? 1. : 1. / wxy;

            wz = pow(p.second.z(),2);
            wz = wz == 0. ? 1. : 1. / wz;

            err_x += wxy * pow(p.first.x() - x, 2);
            err_y += wxy * pow(p.first.y() - y, 2);
            err_z += wz  * pow(p.first.z() - z, 2);
         }
//         err(0,0) = std::sqrt( err_x / (s_wxy * (points.size() -1) ));
//         err(1,1) = std::sqrt( err_y / (s_wxy * (points.size() -1) ));
//         err(2,2) = std::sqrt( err_z / (s_wz *  (points.size() -1) ));
         //isValid = true;
         float dist = 0; 
         //float varikk = 0;
         //float var_x = 0, var_y = 0, var_z = 0;
         //std::vector<float> weights;
         s_wxy = 0; s_wz = 0;
         float old_x = x, old_y = y, old_z = z;
         x=0; y=0; z=0;
         for (const auto& p : points){ 
              //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
            //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();

//            wxy = p.second.x();
//            wz = p.second.z();
//            dist = std::pow(p.first.x() - old_x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
//            dist += std::pow(p.first.y() - old_y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
//            dist += std::pow(p.first.z() - old_z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
//
            dist =  std::pow(p.first.x() - old_x, 2); 
            dist += std::pow(p.first.y() - old_y, 2); 
            dist += std::pow(p.first.z() - old_z, 2);
            //weights.push_back(dist)
            dist = dist == 0 ? 10000 : 1. / dist;

//            wxy = pow(p.second.x(),2);
//            wxy = wxy == 0. ? 1. : 1. / wxy;
//
//            wz = pow(p.second.z(),2);
//            wz = wz == 0. ? 1. : 1. / wz;

            x += p.first.x() * dist;
            y += p.first.y() * dist;
            z += p.first.z() * dist;

            s_wxy += dist ;
            s_wz  += dist ;

            /*
            dist = std::pow(p.first.x() - x, 2) / p.first.x();
            dist += std::pow(p.first.y() - y, 2) / p.first.y();
            dist += std::pow(p.first.z() - z, 2) / p.first.z(); 
            */
            //chi2 += dist;
            //ndof ++;
            
           //var_x += wxy * pow(p.first.x() - x, 2)
         }
         if ( s_wxy != 0. && s_wz != 0. ){
             x /= s_wxy;     
             y /= s_wxy;     
             z /= s_wz;     
             err_x=0; err_y=0; err_z=0;
             for (const auto& p : points){ 

                dist = std::pow(p.first.x() - old_x, 2); 
                dist += std::pow(p.first.y() - old_y, 2);
                dist += std::pow(p.first.z() - old_z, 2);
                dist = dist == 0 ? 10000 : 1. / dist;

//                wxy = pow(p.second.x(),2);
//                wxy = wxy == 0. ? 1. : 1. / wxy;
//
//                wz = pow(p.second.z(),2);
//                wz = wz == 0. ? 1. : 1. / wz;

                err_x += dist * pow(p.first.x() - x, 2);
                err_y += dist * pow(p.first.y() - y, 2);
                err_z += dist * pow(p.first.z() - z, 2);
             }
             err(0,0) = std::sqrt( err_x / (s_wxy * (points.size() -1) )) / 1;
             err(1,1) = std::sqrt( err_y / (s_wxy * (points.size() -1) )) / 1;
             err(2,2) = std::sqrt( err_z / (s_wz *  (points.size() -1) )) / 1;
             //isValid = true;
             for (const auto& p : points){
                  //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
                //GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                wxy = p.second.x();
                wz = p.second.z();
                dist = std::pow(p.first.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
                dist += std::pow(p.first.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
                dist += std::pow(p.first.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
                //weights.push_back(dist)
                /*
                x += p.first.x() * dist * p.second.x();
                y += p.first.y() * dist * p.second.x();
                z += p.first.z() * dist * p.second.z();
                s_wxy += dist * p.second.x();
                s_wz +=  dist * p.second.z();
                */
                /*
                dist = std::pow(p.first.x() - x, 2) / p.first.x();
                dist += std::pow(p.first.y() - y, 2) / p.first.y();
                dist += std::pow(p.first.z() - z, 2) / p.first.z(); 
                */
                chi2 += dist;
                ndof ++;
                
               //var_x += wxy * pow(p.first.x() - x, 2)
             }
            TransientVertex v(GlobalPoint(x,y,z), err, *iclus, chi2, ndof);
            return v;
        }
    }
    err(0,0) = 10000;
    err(1,1) = 10000;
    err(2,2) = 10000;
    return TransientVertex(GlobalPoint(0,0,0), err, *iclus, 0, 0);
}


PrimaryVertexProducerDumbFitter::~PrimaryVertexProducerDumbFitter() {
  if (theTrackFilter)
    delete theTrackFilter;
  if (theTrackClusterizer)
    delete theTrackClusterizer;
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    if (algorithm->fitter)
      delete algorithm->fitter;
    if (algorithm->vertexSelector)
      delete algorithm->vertexSelector;
  }
}

void PrimaryVertexProducerDumbFitter::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // get the BeamSpot, it will always be needed, even when not used as a constraint
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(bsToken, recoBeamSpotHandle);
  if (recoBeamSpotHandle.isValid()) {
    beamSpot = *recoBeamSpotHandle;
  } else {
    edm::LogError("UnusableBeamSpot") << "No beam spot available from EventSetup";
  }

  bool validBS = true;
  VertexState beamVertexState(beamSpot);
  if ((beamVertexState.error().cxx() <= 0.) || (beamVertexState.error().cyy() <= 0.) ||
      (beamVertexState.error().czz() <= 0.)) {
    validBS = false;
    edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors " << beamVertexState.error().matrix();
  }

  //if this is a recovery iteration, check if we already have a valid PV
  if (fRecoveryIteration) {
    auto const& oldVertices = iEvent.get(recoveryVtxToken);
    //look for the first valid (not-BeamSpot) vertex
    for (auto const& old : oldVertices) {
      if (!(old.isFake())) {
        //found a valid vertex, write the first one to the collection and return
        //otherwise continue with regular vertexing procedure
        auto result = std::make_unique<reco::VertexCollection>();
        result->push_back(old);
        iEvent.put(std::move(result), algorithms.begin()->label);
        return;
      }
    }
  }

  // get RECO tracks from the event
  // `tks` can be used as a ptr to a reco::TrackCollection
  edm::Handle<reco::TrackCollection> tks;
  iEvent.getByToken(trkToken, tks);

  // interface RECO tracks to vertex reconstruction
  const auto& theB = &iSetup.getData(theTTBToken);
  std::vector<reco::TransientTrack> t_tks;

  if (f4D) {
    edm::Handle<edm::ValueMap<float> > trackTimesH;
    edm::Handle<edm::ValueMap<float> > trackTimeResosH;
    iEvent.getByToken(trkTimesToken, trackTimesH);
    iEvent.getByToken(trkTimeResosToken, trackTimeResosH);
    t_tks = (*theB).build(tks, beamSpot, *(trackTimesH.product()), *(trackTimeResosH.product()));
  } else {
    t_tks = (*theB).build(tks, beamSpot);
  }
  if (fVerbose) {
    std::cout << "RecoVertex/PrimaryVertexProducer"
              << "Found: " << t_tks.size() << " reconstructed tracks"
              << "\n";
  }

  // select tracks
  std::vector<reco::TransientTrack>&& seltks = theTrackFilter->select(t_tks);
//  if (seltks.size() > 4096) {
//    std::cout << "RecoVertex/PrimaryVertexProducer"
//              << "Found: " << seltks.size() << " reconstructed tracks"
//              << "\n";
//    }

  // clusterize tracks in Z
  std::vector<std::vector<reco::TransientTrack> >&& clusters = theTrackClusterizer->clusterize(seltks);

  if (fVerbose) {
    std::cout << " clustering returned  " << clusters.size() << " clusters  from " << seltks.size()
              << " selected tracks" << std::endl;
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);

    std::vector<TransientVertex> pvs;
    for (std::vector<std::vector<reco::TransientTrack> >::const_iterator iclus = clusters.begin();
         iclus != clusters.end();
         iclus++) {
         /*
         GlobalPoint p;
         GlobalError e;
         if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
            p = itrack.impacPointState().globalPosition();
            e = itrack.impacPointState().globalPosition();
         }
         else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {


         }
         */
         //if (validBS && !validBS) std::cout << "ciao" << std::endl;
         std::vector<std::pair<GlobalPoint, GlobalPoint>> points;
         //std::vector<GlobalPoint> errors;
         if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
             for (const auto& itrack : *iclus){ 
                    GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                    GlobalPoint err(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.track().dzError());
                    std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                    points.push_back(p2);
                    //errors.push_back(err);
             }
            TransientVertex v = weightedMean1Iter(points, iclus);
            if ((v.positionError().matrix())(0,0) != 10000) pvs.push_back(v);
            //pvs.push_back(v);
         } 
         else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {
            for (const auto& itrack : *iclus){ 
                    GlobalPoint p = itrack.impactPointState().globalPosition();
                    GlobalPoint err(itrack.track().dxyError(), itrack.track().dxyError(), itrack.track().dzError());
                    std::pair<GlobalPoint, GlobalPoint> p2(p, err);
                    points.push_back(p2);
            }
            //std::vector<PointAndDistance> vgp;
            //for (auto i=0; i < int( iclus->size() / 2) ; i++){
            //for (auto i=0; i < int( iclus->size() - 1) ; i++){
            /*
            unsigned int theNPairs = 30; 
            unsigned int n_tracks = (2 * (unsigned int)(theNPairs)) < iclus->size() ? 2 * theNPairs : iclus->size();

            std::vector<reco::TransientTrack> goodtracks(*iclus);
            std::partial_sort(goodtracks.begin(), goodtracks.begin() + n_tracks, goodtracks.end(), CompareTwoTracks());
            goodtracks.erase(goodtracks.begin() + n_tracks, goodtracks.end());
            unsigned int t_first = 0;
            unsigned int t_interval = goodtracks.size() / 2;
            unsigned int lim = goodtracks.size() - 1;

            // the 'direction' false: intervals will expand
            //   // true: intervals will shrink
            bool dir = false;
            while (points.size() < theNPairs){
                unsigned int i1 = t_first;
                unsigned int i2 = t_first + t_interval;
                TwoTrackMinimumDistance ttmd;
                reco::TransientTrack tt1 = goodtracks[i1];
                reco::TransientTrack tt2 = goodtracks[i2];
                //bool status = ttmd.calculate()[i1].impactPointState(), (*iclus)[i2].impactPointState());
                bool status = ttmd.calculate(tt1.impactPointState(), tt2.impactPointState());
                if (status) {
                    std::pair<GlobalPoint, GlobalPoint> pts = ttmd.points();
                    float x = (pts.second.x() + pts.first.x()) / 2.;
                    float y = (pts.second.y() + pts.first.y()) / 2.;
                    float z = (pts.second.z() + pts.first.z()) / 2.;
                    GlobalPoint p(x, y, z);

                    float xerr = std::sqrt( (pow(pts.first.x() - p.x(), 2) + pow(pts.second.x() - p.x(), 2) ) / 2);
                    float yerr = std::sqrt( (pow(pts.first.y() - p.y(), 2) + pow(pts.second.y() - p.y(), 2) ) / 2);
                    float zerr = std::sqrt( (pow(pts.first.z() - p.z(), 2) + pow(pts.second.z() - p.z(), 2) ) / 2);
                    //float xerr=0, yerr=0, zerr=0;
                    points.push_back(std::pair<GlobalPoint, GlobalPoint>(p, GlobalPoint(xerr,yerr,zerr)));
                    //PointAndDistance v((pts.second + pts.first) / 2., (pts.second - pts.first).mag());
                    //vgp.push_back(v); 
                }
                if ((t_first + t_interval) < lim) {
                  t_first++;
                } else if (dir) {
                  t_first = 0;
                  t_interval--;
                  if (t_interval == 0) {
                    break;
                  }
                } else {
                  t_first = 0;
                  t_interval++;
                  if (t_interval == goodtracks.size()) {
                    dir = true;
                    t_interval = goodtracks.size() / 2 - 1;
                  }
                }
            }
            */
            TransientVertex v = weightedMean1Iter(points, iclus);
            if ((v.positionError().matrix())(0,0) != 10000) pvs.push_back(v);

         }

           /* 

             float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0, chi2=0;
             unsigned int ndof = 0;
             //bool isValid = false;
             for (const auto& itrack : *iclus){ 
             //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
                    GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                    //GlobalError e(itrack.stateAtBeamLine().transverseImpactParameter().error(), itrack.stateAtBeamLine().transverseImpactParameter().error(),  itrack.track().dzError());

                    wxy = pow(itrack.stateAtBeamLine().transverseImpactParameter().error(),2);
                    wxy = wxy == 0. ? 1. : 1. / wxy;
                    x += p.x() * wxy;
                    y += p.y() * wxy;
                    s_wxy += wxy;

                    wz = pow(itrack.track().dzError(),2);
                    wz = wz == 0. ? 1. : 1. / wz;
                    z += p.z() * wz;
                    s_wz += wz;
             }

             //auto err = v.positionError().matrix();
             AlgebraicSymMatrix33 err;
             if ( s_wxy != 0. && s_wz != 0. ){
                 x /= s_wxy;     
                 y /= s_wxy;     
                 err(0,0) = std::sqrt(1. / s_wxy);
                 err(1,1) = std::sqrt(1. / s_wxy);
                 z /= s_wz;     
                 err(2,2) = std::sqrt(1. / s_wz);
                 //isValid = true;
                 float dist = 0; 
                 for (const auto& itrack : *iclus){ 
                      //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
                    GlobalPoint p =  itrack.stateAtBeamLine().trackStateAtPCA().position();
                    wxy = itrack.stateAtBeamLine().transverseImpactParameter().error();
                    wz = itrack.track().dzError();

                      dist = std::pow(p.x() - x, 2) / ( std::pow(wxy, 2) + std::pow( err(0,0),2) );
                      dist += std::pow(p.y() - y, 2) / ( std::pow(wxy, 2) + std::pow( err(1,1),2) );
                      dist += std::pow(p.z() - z, 2) / ( std::pow(wz, 2) + std::pow( err(2,2),2) ); 
                      chi2 += dist;
                      ndof ++;
                 }
                TransientVertex v(GlobalPoint(x,y,z), err, (*iclus), chi2, ndof);
                pvs.push_back(v);
                
             }
         
         }
         else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {

             float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0, chi2=0;
             unsigned int ndof = 0;
             //bool isValid = false;
             for (const auto& itrack : *iclus){ 
             //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
                    wxy = pow(itrack.track().dxyError(),2);
                    wxy = wxy == 0. ? 1. : 1. / wxy;
                    x += itrack.impactPointState().globalPosition().x() * wxy;
                    y += itrack.impactPointState().globalPosition().y() * wxy;
                    s_wxy += wxy;

                    wz = pow(itrack.track().dzError(),2);
                    wz = wz == 0. ? 1. : 1. / wz;
                    z += itrack.impactPointState().globalPosition().z() * wz;
                    s_wz += wz;
             }

             //auto err = v.positionError().matrix();
             AlgebraicSymMatrix33 err;
             if ( s_wxy != 0. && s_wz != 0. ){
                 x /= s_wxy;     
                 y /= s_wxy;     
                 err(0,0) = std::sqrt(1. / s_wxy);
                 err(1,1) = std::sqrt(1. / s_wxy);
                 z /= s_wz;     
                 err(2,2) = std::sqrt(1. / s_wz);
                 //isValid = true;
                 float dist = 0; 
                 for (const auto& itrack : *iclus){ 
                      //dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) + std::pow(itrack.impactPointState().globalPosition().y() - y, 2) + std::pow(itrack.impactPointState().globalPosition().z() - z, 2); 
                      dist = std::pow(itrack.impactPointState().globalPosition().x() - x, 2) / ( std::pow(itrack.track().dxyError(), 2) + std::pow( err(0,0),2) );
                      dist += std::pow(itrack.impactPointState().globalPosition().y() - y, 2) / ( std::pow(itrack.track().dxyError(), 2) + std::pow( err(1,1),2) );
                      dist += std::pow(itrack.impactPointState().globalPosition().z() - z, 2) / ( std::pow(itrack.track().dzError(), 2) + std::pow( err(2,2),2) ); 
                      chi2 += dist;
                      ndof ++;
                 }
                TransientVertex v(GlobalPoint(x,y,z), err, (*iclus), chi2, ndof);
                pvs.push_back(v);
                
             }
         }
        */
         
//         if (iclus->size()>1){
//             x /= iclus->size();     
//             y /= iclus->size();     
//             z /= iclus->size();     
//         }
         //TransientVertex v2(GlobalPoint(x,y,z), v.positionError(), (*iclus), v.totalChiSquared());
      /*
      double sumwt = 0.;
      double sumwt2 = 0.;
      double sumw = 0.;
      double meantime = 0.;
      double vartime = 0.;
      if (f4D) {
        for (const auto& tk : *iclus) {
          const double time = tk.timeExt();
          const double err = tk.dtErrorExt();
          const double inverr = err > 0. ? 1.0 / err : 0.;
          const double w = inverr * inverr;
          sumwt += w * time;
          sumwt2 += w * time * time;
          sumw += w;
        }
        meantime = sumwt / sumw;
        double sumsq = sumwt2 - sumwt * sumwt / sumw;
        double chisq = iclus->size() > 1 ? sumsq / double(iclus->size() - 1) : sumsq / double(iclus->size());
        vartime = chisq / sumw;
      }

      TransientVertex v;
      if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
        v = algorithm->fitter->vertex(*iclus, beamSpot);
      } else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {
        v = algorithm->fitter->vertex(*iclus);
      }  // else: no fit ==> v.isValid()=False
     if (v.isValid()){
     float x=0, y=0, z=0, s_wxy=0, s_wz=0, wxy=0, wz=0;
         for (const auto& itrack : *iclus){ 
         //std::vector<reco::TransientTrack>::const_iterator itrack = iclus.begin(); itrack!= iclus.end(); itrack++) {
                wxy = pow(itrack.track().dxyError(),2);
                wxy = wxy == 0. ? 1. : 1. / wxy;
                x += itrack.impactPointState().globalPosition().x() * wxy;
                y += itrack.impactPointState().globalPosition().y() * wxy;
                s_wxy += wxy;

                wz = pow(itrack.track().dzError(),2);
                wz = wz == 0. ? 1. : 1. / wz;
                z += itrack.impactPointState().globalPosition().z() * wz;
                s_wz += wz;
         }

         auto err = v.positionError().matrix();

         if (s_wxy != 0.){
             x /= s_wxy;     
             y /= s_wxy;     
             err(0,0) = std::sqrt(1. / s_wxy);
             err(1,1) = std::sqrt(1. / s_wxy);
         }
         if (s_wz != 0.) {
             z /= s_wz;     
             err(2,2) = std::sqrt(1. / s_wz);
         }
//         if (iclus->size()>1){
//             x /= iclus->size();     
//             y /= iclus->size();     
//             z /= iclus->size();     
//         }
         //TransientVertex v2(GlobalPoint(x,y,z), v.positionError(), (*iclus), v.totalChiSquared());
         v = TransientVertex(GlobalPoint(x,y,z), err, v.originalTracks(), v.totalChiSquared());
     }
      // 4D vertices: add timing information
      if (f4D and v.isValid()) {
        auto err = v.positionError().matrix4D();
        err(3, 3) = vartime;
        auto trkWeightMap3d = v.weightMap();  // copy the 3d-fit weights
        v = TransientVertex(v.position(), meantime, err, v.originalTracks(), v.totalChiSquared(), v.degreesOfFreedom());
        v.weightMap(trkWeightMap3d);
      }

      if (fVerbose) {
        if (v.isValid()) {
          std::cout << "x,y,z";
          if (f4D)
            std::cout << ",t";
          std::cout << "=" << v.position().x() << " " << v.position().y() << " " << v.position().z();
          if (f4D)
            std::cout << " " << v.time();
          std::cout << " cluster size = " << (*iclus).size() << std::endl;
        } else {
          std::cout << "Invalid fitted vertex,  cluster size=" << (*iclus).size() << std::endl;
        }
      }

      if (v.isValid() && (v.degreesOfFreedom() >= algorithm->minNdof) &&
          (!validBS || (*(algorithm->vertexSelector))(v, beamVertexState)))
        pvs.push_back(v);
      */
    }  // end of cluster loop

    if (fVerbose) {
      std::cout << "PrimaryVertexProducerAlgorithm::vertices  candidates =" << pvs.size() << std::endl;
    }

    if (clusters.size() > 2 && clusters.size() > 2 * pvs.size())
      edm::LogWarning("PrimaryVertexProducer")
          << "more than half of candidate vertices lost " << pvs.size() << ' ' << clusters.size();

    if (pvs.empty() && seltks.size() > 5)
      edm::LogWarning("PrimaryVertexProducer")
          << "no vertex found with " << seltks.size() << " tracks and " << clusters.size() << " vertex-candidates";

    // sort vertices by pt**2  vertex (aka signal vertex tagging)
    if (pvs.size() > 1) {
      sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
    }

    // convert transient vertices returned by the theAlgo to (reco) vertices
    for (std::vector<TransientVertex>::const_iterator iv = pvs.begin(); iv != pvs.end(); iv++) {
      reco::Vertex v = *iv;
      vColl.push_back(v);
    }

    if (vColl.empty()) {
      GlobalError bse(beamSpot.rotatedCovariance3D());
      if ((bse.cxx() <= 0.) || (bse.cyy() <= 0.) || (bse.czz() <= 0.)) {
        AlgebraicSymMatrix33 we;
        we(0, 0) = 10000;
        we(1, 1) = 10000;
        we(2, 2) = 10000;
        vColl.push_back(reco::Vertex(beamSpot.position(), we, 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << "Beamspot with invalid errors " << bse.matrix() << std::endl;
          std::cout << "Will put Vertex derived from dummy-fake BeamSpot into Event.\n";
        }
      } else {
        vColl.push_back(reco::Vertex(beamSpot.position(), beamSpot.rotatedCovariance3D(), 0., 0., 0));
        if (fVerbose) {
          std::cout << "RecoVertex/PrimaryVertexProducer: "
                    << " will put Vertex derived from BeamSpot into Event.\n";
        }
      }
    }

    if (fVerbose) {
      int ivtx = 0;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        std::cout << "recvtx " << ivtx++ << "#trk " << std::setw(3) << v->tracksSize() << " chi2 " << std::setw(4)
                  << v->chi2() << " ndof " << std::setw(3) << v->ndof() << " x " << std::setw(6) << v->position().x()
                  << " dx " << std::setw(6) << v->xError() << " y " << std::setw(6) << v->position().y() << " dy "
                  << std::setw(6) << v->yError() << " z " << std::setw(6) << v->position().z() << " dz " << std::setw(6)
                  << v->zError();
        if (f4D) {
          std::cout << " t " << std::setw(6) << v->t() << " dt " << std::setw(6) << v->tError();
        }
        std::cout << std::endl;
      }
    }
      int ivtx = 0;
        std::cout << "recvtx,#trk,chi2,ndof,x,dx,y,dy,z,dz" << std::endl;
      for (reco::VertexCollection::const_iterator v = vColl.begin(); v != vColl.end(); ++v) {
        std::cout << ivtx++ << "," << v->tracksSize() << "," << v->chi2() << ","  << v->ndof() << ","  << v->position().x()
                  << ","  << v->xError() << ","  << v->position().y() << ","
                   << v->yError() << ","  << v->position().z() << "," 
                  << v->zError();
        std::cout << std::endl;

      }
    iEvent.put(std::move(result), algorithm->label);
  }
}


void PrimaryVertexProducerDumbFitter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // offlinePrimaryVertices
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<double>("maxDistanceToBeam", 1.0);
    vpsd1.add<std::string>("algorithm", "AdaptiveVertexFitter");
    vpsd1.add<bool>("useBeamConstraint", false);
    vpsd1.add<std::string>("label", "");
    vpsd1.add<double>("chi2cutoff", 2.5);
    vpsd1.add<double>("minNdof", 0.0);
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(2);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", false);
      temp2.addParameter<std::string>("label", "");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("minNdof", 0.0);
      temp1.push_back(temp2);
    }
    {
      edm::ParameterSet temp2;
      temp2.addParameter<double>("maxDistanceToBeam", 1.0);
      temp2.addParameter<std::string>("algorithm", "AdaptiveVertexFitter");
      temp2.addParameter<bool>("useBeamConstraint", true);
      temp2.addParameter<std::string>("label", "WithBS");
      temp2.addParameter<double>("chi2cutoff", 2.5);
      temp2.addParameter<double>("minNdof", 2.0);
      temp1.push_back(temp2);
    }
    desc.addVPSet("vertexCollections", vpsd1, temp1);
  }
  desc.addUntracked<bool>("verbose", false);
  {
    edm::ParameterSetDescription psd0;
    TrackFilterForPVFinding::fillPSetDescription(psd0);
    psd0.add<int>("numTracksThreshold", 0);  // HI only
    desc.add<edm::ParameterSetDescription>("TkFilterParameters", psd0);
  }
  desc.add<edm::InputTag>("beamSpotLabel", edm::InputTag("offlineBeamSpot"));
  desc.add<edm::InputTag>("TrackLabel", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("TrackTimeResosLabel", edm::InputTag("dummy_default"));  // 4D only
  desc.add<edm::InputTag>("TrackTimesLabel", edm::InputTag("dummy_default"));      // 4D only

  {
    edm::ParameterSetDescription psd0;
    {
      edm::ParameterSetDescription psd1;
      DAClusterizerInZT_vect::fillPSetDescription(psd1);
      psd0.add<edm::ParameterSetDescription>("TkDAClusParameters", psd1);

      edm::ParameterSetDescription psd2;
      GapClusterizerInZ::fillPSetDescription(psd2);
      psd0.add<edm::ParameterSetDescription>("TkGapClusParameters", psd2);
    }
    psd0.add<std::string>("algorithm", "DA_vect");
    desc.add<edm::ParameterSetDescription>("TkClusParameters", psd0);
  }

  desc.add<bool>("isRecoveryIteration", false);
  desc.add<edm::InputTag>("recoveryVtxCollection", {""});

  descriptions.add("primaryVertexProducerDumbFitter", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexProducerDumbFitter);
