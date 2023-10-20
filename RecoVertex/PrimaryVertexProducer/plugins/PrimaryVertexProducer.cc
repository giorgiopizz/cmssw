#include "RecoVertex/PrimaryVertexProducer/interface/PrimaryVertexProducer.h"
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
#include <algorithm>
#include <numeric>

constexpr double vertexSize_ = 1 * 0.006;  //0.006;
constexpr uint neighboursCore = 2;
constexpr double tmin = 2.0;

double myDistance(double z1, double err1, double z2, double err2) {
  return (std::pow(z1 - z2, 2) / (std::pow(err1, 2) + std::pow(err2, 2)));
}
bool passDistance(double z1, double err1, double z2, double err2) { return myDistance(z1, err1, z2, err2) < tmin; }

std::vector<uint> find_neighbours(
    std::vector<double>& zs, std::vector<double>& errs, std::vector<bool>& core, uint i, bool onlyCore) {
  std::vector<uint> ret;
  for (auto j = 0U; j < zs.size(); j++) {
    if (i == j)
      continue;
    if (onlyCore && !core[j])
      continue;
    if (passDistance(zs[i], errs[i], zs[j], errs[j]))
      ret.push_back(j);
  }
  return ret;
}

std::tuple<double, double> computeClusterCenter(std::vector<uint>& cluster, std::vector<double>& zs, std::vector<double>& errs) {
  double center = 0;
  double weights = 0;
  for (auto i = 0U; i < cluster.size(); i++) {
    uint ind = cluster[i];
    center += zs[ind] / std::pow(errs[ind], 2);
    weights += 1.0 / std::pow(errs[ind], 2);
  }
  return {center / weights, 1.0 / std::sqrt(weights)};
}

std::vector<uint> updateCluster(std::vector<uint>&& cluster,
                                std::vector<bool>& assigned,
                                std::vector<double>& zs,
                                std::vector<double>& errs,
                                std::vector<bool>& core,
                                uint i,
                                bool onlyCore,
                                double z_cluster,
                                double err_cluster) {
  // for (auto i = 0; i < cluster.size())

  auto neighbours = find_neighbours(zs, errs, core, i, onlyCore);  // find all neighbours that are core
  for (auto j = 0U; j < neighbours.size(); j++) {
    // if (std::find(cluster.begin(), cluster.end(), neighbours[j]) == cluster.end()) continue;
    uint ind = neighbours[j];
    if (assigned[ind])
      continue;
    if (!passDistance(z_cluster, err_cluster, zs[ind], errs[ind]))
      continue;  // do not add tracks not compatible with cluster
    assigned[ind] = true;
    cluster.push_back(neighbours[j]);
    auto [ center, err ] = computeClusterCenter(cluster, zs, errs);
    // cluster = updateCluster(std::move(cluster), assigned, zs, errs, core, neighbours[j], onlyCore, center, err);
    cluster = updateCluster(std::move(cluster), assigned, zs, errs, core, neighbours[j], onlyCore, center, err);
  }

  return cluster;
}

template <typename T>
std::vector<T> reorder(std::vector<T>&& vec, std::vector<uint>& ind) {
  std::vector<T> ret(vec.size());
  for (auto i = 0U; i < vec.size(); i++) {
    ret[i] = vec[ind[i]];
  }
  return ret;
}

std::vector<std::vector<reco::TransientTrack>> clusterize(std::vector<reco::TransientTrack>& tracks) {
  // prepare track data for clustering
  // track_t tks;
  // double sumtkwt = 0.;
  // auto const& t_mom = (*it).stateAtBeamLine().trackStateAtPCA().momentum();
  std::vector<double> zs, errs;
  std::vector<const reco::TransientTrack*> tt;
  for (auto it = tracks.begin(); it != tracks.end(); it++) {
    if (!(*it).isValid())
      continue;
    double t_tkwt = 1.;
    double t_z = ((*it).stateAtBeamLine().trackStateAtPCA()).position().z();
    if (std::fabs(t_z) > 1000.)
      continue;
    auto const& t_mom = (*it).stateAtBeamLine().trackStateAtPCA().momentum();
    //  get the beam-spot
    reco::BeamSpot beamspot = (it->stateAtBeamLine()).beamSpot();
    double t_dz2 = std::pow((*it).track().dzError(), 2)  // track errror
                   + (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
                         std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
                   + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
    if (edm::isNotFinite(t_dz2) || t_dz2 < std::numeric_limits<double>::min())
      continue;
    zs.push_back(t_z);
    errs.push_back(std::sqrt(t_dz2));
    tt.push_back(&(*it));
  }
  assert(zs.size() == errs.size());

  //sort tracks
  std::vector<uint> indices(zs.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](uint A, uint B) -> bool { return zs[A] < zs[B]; });
  zs = reorder<double>(std::move(zs), indices);
  errs = reorder<double>(std::move(errs), indices);
  tt = reorder<const reco::TransientTrack*>(std::move(tt), indices);
  // for (auto i = 0U; i < zs.size(); i++){
  //   std::cout << "Track "<< i << " "<< zs[i] << "\n";
  // }

  std::vector<bool> core(zs.size());
  std::vector<bool> assigned(zs.size(), false);
  std::vector<double> neighbours(zs.size(), 0.0);
  for (auto i = 0U; i < zs.size(); i++) {
    // uint neighbours = 0;
    for (auto j = 0U; j < zs.size(); j++) {
      if (i == j)
        continue;
      if (passDistance(zs[i], errs[i], zs[j], errs[j]))
        neighbours[i]++;
    }
    core[i] = (neighbours[i] >= neighboursCore);
    neighbours[i] *= 1 / std::pow(errs[i], 1/4);
    // core.push_back(neighbours >= neighboursCore);
    // assigned.push_back(false);
  }

  // indices.clear();
  // indices.resize(zs.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(), [&](uint A, uint B) -> bool { return neighbours[A] > neighbours[B]; });
  std::vector<std::vector<uint>> clusters;
  // add only core elements to clusters
  for (auto i = 0U; i < zs.size(); i++) {
    uint ind = indices[i];
    if (!core[ind])
      continue;
    if (assigned[ind])
      continue;
    assigned[ind] = true;
    std::vector<uint> cluster;
    cluster.push_back(ind);
    // auto [ center, err ] = computeClusterCenter(cluster, zs, errs);
    cluster = updateCluster(std::move(cluster), assigned, zs, errs, core, ind, true, zs[ind], errs[ind]);

    clusters.push_back(cluster);
  }

  std::vector<double> z_clusters(clusters.size());
  std::vector<double> err_clusters(clusters.size());
  for (auto i = 0U; i < z_clusters.size(); i++) {
    // double center = 0;
    // double weights = 0;
    // for (auto j = 0U; j < clusters[i].size(); j++) {
    //   auto idx = clusters[i][j];
    //   center += zs[idx] / std::pow(errs[idx], 2);
    //   weights += 1.0 / std::pow(errs[idx], 2);
    // }
    auto [center, err] = computeClusterCenter(clusters[i], zs, errs);
    z_clusters[i] = center;
    err_clusters[i] = err;
//    z_clusters[i] = center / weights;
//    err_clusters[i] = 1.0 / std::sqrt(weights);
  }

  // add non core elements
  for (auto i = 0U; i < zs.size(); i++) {
    if (core[i])
      continue;
    if (assigned[i])
      continue;
    double minDist = 10000;
    uint minClust = 0;
    for (auto clustIdx = 0U; clustIdx < clusters.size(); clustIdx++) {
      double currDist = myDistance(zs[i], errs[i], z_clusters[clustIdx], err_clusters[clustIdx]);
      if (currDist < minDist) {
        minDist = currDist;
        minClust = clustIdx;
      }
    }
    if (minDist < tmin) {
      clusters[minClust].push_back(i);
      assigned[i] = true;
    }
  }


  // // split vertices
  // for (auto i = 0U; i < z_clusters.size(); i++) {


  /*
  // add non core elements
  for (auto clustIdx = 0U; clustIdx < clusters.size(); clustIdx++) {
    uint clustSize = clusters[clustIdx].size();
    for (auto i = 0U; i < clustSize; i++) {
      uint idx = clusters[clustIdx][i];
      if (!core[idx])
        continue;                                                     // should I really check?
      auto neighbours = find_neighbours(zs, errs, core, idx, false);  // find all neighbours that are core
      for (auto newIdx : neighbours) {
        if (assigned[newIdx])
          continue;
        clusters[clustIdx].push_back(newIdx);
        assigned[newIdx] = true;
      }
    }
  }
  */
  // for (auto i = 0U; i < zs.size(); i++){
  //   if (assigned[i]) continue;

  //   auto neighbours = find_neighbours(zs, errs, core, i, onlyCore); // find all neighbours that are core

  //   cluster = updateCluster(std::move(cluster), assigned, zs, errs, core, i, true);

  //   clusters.push_back(cluster);
  // }

  std::vector<std::vector<reco::TransientTrack>> ret;
  for (auto clustIdx = 0U; clustIdx < clusters.size(); clustIdx++) {
    std::vector<reco::TransientTrack> cluster;
    // std::cout << "\nBegin of cluster " << clustIdx << "\n";
    for (auto i = 0U; i < clusters[clustIdx].size(); i++) {
      uint idx = clusters[clustIdx][i];
      // cluster.push_back(tracks[idx]);
      cluster.push_back(*tt[idx]);
      // std::cout << "Track " << zs[idx] << " " << errs[idx] << "\n";
    }
    // if (clusters[clustIdx].size()>=2){
    // std::cout << "Track 0-1 distance" << zs[] << " " << errs[i] << "\n";
    // }
    ret.push_back(cluster);
  }

  return ret;
}

PrimaryVertexProducer::PrimaryVertexProducer(const edm::ParameterSet& conf)
    : theTTBToken(esConsumes(edm::ESInputTag("", "TransientTrackBuilder"))), theConfig(conf) {
  fVerbose = conf.getUntrackedParameter<bool>("verbose", false);

  trkToken = consumes<reco::TrackCollection>(conf.getParameter<edm::InputTag>("TrackLabel"));
  bsToken = consumes<reco::BeamSpot>(conf.getParameter<edm::InputTag>("beamSpotLabel"));
  f4D = false;
  weightFit = false;

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
    trkTimesToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("TrackTimesLabel"));
    trkTimeResosToken = consumes<edm::ValueMap<float>>(conf.getParameter<edm::InputTag>("TrackTimeResosLabel"));
  }

  // select and configure the vertex fitters
  std::vector<edm::ParameterSet> vertexCollections =
      conf.getParameter<std::vector<edm::ParameterSet>>("vertexCollections");

  for (std::vector<edm::ParameterSet>::const_iterator algoconf = vertexCollections.begin();
       algoconf != vertexCollections.end();
       algoconf++) {
    algo algorithm;
    std::string fitterAlgorithm = algoconf->getParameter<std::string>("algorithm");
    if (fitterAlgorithm == "KalmanVertexFitter") {
      algorithm.fitter = new KalmanVertexFitter();
    } else if (fitterAlgorithm == "AdaptiveVertexFitter") {
      algorithm.fitter = new AdaptiveVertexFitter(GeometricAnnealing(algoconf->getParameter<double>("chi2cutoff")));
    } else if (fitterAlgorithm == "WeightedMeanFitter") {
      algorithm.fitter = nullptr;
      weightFit = true;
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

PrimaryVertexProducer::~PrimaryVertexProducer() {
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

void PrimaryVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
    edm::Handle<edm::ValueMap<float>> trackTimesH;
    edm::Handle<edm::ValueMap<float>> trackTimeResosH;
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

  // clusterize tracks in Z
  std::vector<std::vector<reco::TransientTrack>>&& clusters = theTrackClusterizer->clusterize(seltks);

  //std::vector<std::vector<reco::TransientTrack>>&& clusters = clusterize(seltks);

  // for (auto clustIdx = 0U; clustIdx < clusters.size(); clustIdx++) {
  //   // std::vector<reco::TransientTrack> cluster;
  //   std::cout << "\nBegin of cluster " << clustIdx << "\n";
  //   for (auto i = 0U; i < clusters[clustIdx].size(); i++) {
  //     auto track = clusters[clustIdx][i];
  //     auto z = ((track).stateAtBeamLine().trackStateAtPCA()).position().z();
  //     auto const& t_mom = track.stateAtBeamLine().trackStateAtPCA().momentum();
  //     reco::BeamSpot beamspot = (track.stateAtBeamLine()).beamSpot();
  //     double t_dz2 = std::pow(track.track().dzError(), 2)  // track errror
  //                    +
  //                    (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
  //                        std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
  //                    + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
  //     std::cout << "Track " << z << " " << t_dz2 << "\n";
  //   }
  // }
  // std::cout << "\n\nCorrect clustering\n";
  // clusters = theTrackClusterizer->clusterize(seltks);
  // for (auto clustIdx = 0U; clustIdx < clusters.size(); clustIdx++) {
  //   // std::vector<reco::TransientTrack> cluster;
  //   std::cout << "\nBegin of cluster " << clustIdx << "\n";
  //   for (auto i = 0U; i < clusters[clustIdx].size(); i++) {
  //     auto track = clusters[clustIdx][i];
  //     auto z = ((track).stateAtBeamLine().trackStateAtPCA()).position().z();
  //     auto const& t_mom = track.stateAtBeamLine().trackStateAtPCA().momentum();
  //     reco::BeamSpot beamspot = (track.stateAtBeamLine()).beamSpot();
  //     double t_dz2 = std::pow(track.track().dzError(), 2)  // track errror
  //                    +
  //                    (std::pow(beamspot.BeamWidthX() * t_mom.x(), 2) + std::pow(beamspot.BeamWidthY() * t_mom.y(), 2)) *
  //                        std::pow(t_mom.z(), 2) / std::pow(t_mom.perp2(), 2)  // beam spot width
  //                    + std::pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
  //     std::cout << "Track " << z << " " << t_dz2 << "\n";
  //   }
  // }

  if (fVerbose) {
    std::cout << " clustering returned  " << clusters.size() << " clusters  from " << seltks.size()
              << " selected tracks" << std::endl;
  }

  // vertex fits
  for (std::vector<algo>::const_iterator algorithm = algorithms.begin(); algorithm != algorithms.end(); algorithm++) {
    auto result = std::make_unique<reco::VertexCollection>();
    reco::VertexCollection& vColl = (*result);

    std::vector<TransientVertex> pvs;
    for (std::vector<std::vector<reco::TransientTrack>>::const_iterator iclus = clusters.begin();
         iclus != clusters.end();
         iclus++) {
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
      if (algorithm->fitter) {
        if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
          v = algorithm->fitter->vertex(*iclus, beamSpot);
        } else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {
          v = algorithm->fitter->vertex(*iclus);
        }  // else: no fit ==> v.isValid()=False
      } else if (weightFit) {
        std::vector<std::pair<GlobalPoint, GlobalPoint>> points;
        if (algorithm->useBeamConstraint && validBS && (iclus->size() > 1)) {
          for (const auto& itrack : *iclus) {
            GlobalPoint p = itrack.stateAtBeamLine().trackStateAtPCA().position();
            GlobalPoint err(itrack.stateAtBeamLine().transverseImpactParameter().error(),
                            itrack.stateAtBeamLine().transverseImpactParameter().error(),
                            itrack.track().dzError());
            std::pair<GlobalPoint, GlobalPoint> p2(p, err);
            points.push_back(p2);
          }

          v = WeightedMeanFitter::weightedMeanOutlierRejectionBeamSpot(points, *iclus, beamSpot);
          if ((v.positionError().matrix())(2, 2) != (WeightedMeanFitter::startError * WeightedMeanFitter::startError))
            pvs.push_back(v);
        } else if (!(algorithm->useBeamConstraint) && (iclus->size() > 1)) {
          for (const auto& itrack : *iclus) {
            GlobalPoint p = itrack.impactPointState().globalPosition();
            GlobalPoint err(itrack.track().dxyError(), itrack.track().dxyError(), itrack.track().dzError());
            std::pair<GlobalPoint, GlobalPoint> p2(p, err);
            points.push_back(p2);
          }

          v = WeightedMeanFitter::weightedMeanOutlierRejection(points, *iclus);
          if ((v.positionError().matrix())(2, 2) != (WeightedMeanFitter::startError * WeightedMeanFitter::startError))
            pvs.push_back(v);  //FIX with constants
        }
      } else
        throw VertexException(
            "PrimaryVertexProducer: Something went wrong. You are not using the weighted mean fit and no algorithm was "
            "selected.");

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

      //for weightFit we have already pushed it above (no timing infomration anyway)
      if (v.isValid() && not weightFit && (v.degreesOfFreedom() >= algorithm->minNdof) &&
          (!validBS || (*(algorithm->vertexSelector))(v, beamVertexState)))
        pvs.push_back(v);
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

    iEvent.put(std::move(result), algorithm->label);
  }
}

void PrimaryVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
    psd0.add<int>("numTracksThreshold", 0);            // HI only
    psd0.add<int>("maxNumTracksThreshold", 10000000);  // HI only
    psd0.add<double>("minPtTight", 0.0);               // HI only
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

  descriptions.add("primaryVertexProducer", desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PrimaryVertexProducer);
