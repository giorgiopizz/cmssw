#include "SimTracker/VertexAssociation/interface/calculateVertexSharedTracks.h"

bool matchSimTrack (const SimTrack& tk1, const SimTrack& tk2){
  return tk1.momentum() == tk2.momentum() && tk1.charge() == tk2.charge();
}

bool matchSimTracks (const std::vector<SimTrack>& tk1, const std::vector<SimTrack>& tk2){
  if (tk1.size() != tk2.size()) return false;
  for (auto i = 0U; i < tk1.size() ; i++){
    if (! matchSimTrack(tk1[i], tk2[i])) return false;
  }
  return true;
}

// bool SimTrack::operator == (const SimTrack& other){
//   return this->momentum() == other.momentum() && this->charge() == other.charge();
// }

bool matchSim(const TrackingParticle & simTrack, const TrackingParticle& daughterTrack){
              // bool ret = true;
              // for 
              return (
                matchSimTracks(simTrack.g4Tracks(), daughterTrack.g4Tracks()) &&
                simTrack.parentVertex() == daughterTrack.parentVertex()
              );
}


unsigned int calculateVertexSharedTracks(const reco::Vertex &recoV,
                                         const TrackingVertex &simV,
                                         const reco::RecoToSimCollection &trackRecoToSimAssociation) {
  unsigned int sharedTracks = 0;
  for (auto iTrack = recoV.tracks_begin(); iTrack != recoV.tracks_end(); ++iTrack) {
    auto found = trackRecoToSimAssociation.find(*iTrack);

    if (found == trackRecoToSimAssociation.end())
      continue;

    // matched TP equal to any TP of sim vertex => increase counter
    for (const auto &tp : found->val) {
      if (std::find_if(simV.daughterTracks_begin(), simV.daughterTracks_end(), [&](const TrackingParticleRef &vtp) {
            // return tp.first == vtp;
            return ((tp.first.get() != nullptr) && (vtp.get() != nullptr)) ? matchSim(*(tp.first.get()), *(vtp.get())) : false;

          }) != simV.daughterTracks_end()) {
        sharedTracks += 1;
        break;
      }
    }
  }
  // std::cout << "Shared tracks reco -> sim: " << sharedTracks << "\n";
  return sharedTracks;
}

// bool matchVal (auto a, auto b){
//   if (abs(b) > 1e-7){
//     return abs(a-b)/b <= 1e-7;
//   }
//   return abs(a-b) <= 1e-8;
// }

// bool matchReco(auto& recoTrack, auto& vtk){
//               // return (matchVal(recoTrack.first.get()->pt(), vtk->pt()) &&
//               // matchVal(recoTrack.first.get()->eta() , vtk->eta()) &&
//               // matchVal(recoTrack.first.get()->phi() , vtk->phi()));
//               return (recoTrack.innerMomentum() == vtk.innerMomentum()) && (recoTrack.innerPosition() == vtk.innerPosition());
//               // return (
//               //   recoTrack.momentum() == vtk.momentum() &&
//               //   recoTrack.vertex() == vtk.vertex() &&
//               //   recoTrack.covariance() == vtk.covariance()
//               // );
              // return (matchVal(recoTrack.first.get()->pt(), vtk->pt()) &&
              // matchVal(recoTrack.first.get()->eta() , vtk->eta()) &&
              // matchVal(recoTrack.first.get()->phi() , vtk->phi()));
              // return (recoTrack.innerMomentum() == vtk.innerMomentum()) && (recoTrack.innerPosition() == vtk.innerPosition());
// }

bool matchReco(const reco::Track & recoTrack, const reco::Track& vtk){
              return (
                recoTrack.momentum() == vtk.momentum() &&
                recoTrack.vertex() == vtk.vertex() &&
                recoTrack.covariance() == vtk.covariance()
              );
}

unsigned int calculateVertexSharedTracks(const TrackingVertex &simV,
                                         const reco::Vertex &recoV,
                                         const reco::SimToRecoCollection &trackSimToRecoAssociation) {
  unsigned int sharedTracks = 0;
  // for (auto& el: trackSimToRecoAssociation){
  //   for (auto & recoTrack: el.val){
  // //     // std::cout << static_cast<edm::RefToBaseVector<reco::Track>>(recoTrack.first).pt() << "\n";
  // //     // std::cout << recoTrack.first.get()->pt() << "\t" << recoTrack.second << "\n";
  //     if (std::find_if(recoV.tracks_begin(), recoV.tracks_end(), [&](const reco::TrackBaseRef &vtk) {
  //           // return recoTrack.first.get()->stateAtBeamLine().trackStateAtPCA().position().z() == vtk->stateAtBeamLine().trackStateAtPCA().position().z();
  //           // return (abs(recoTrack.first.get()->pt() - vtk->pt())/vtk->pt() < 0.001 &&
  //           //   abs(recoTrack.first.get()->eta() - vtk->eta())/vtk->eta() < 0.001 &&
  //           //   abs(recoTrack.first.get()->phi() - vtk->phi())/vtk->phi() < 0.001);
  //           return matchReco(recoTrack, vtk);
  //         }) != recoV.tracks_end()) {
  //           std::cout << "Matched with reco\n";
  //         }
  //   }
  // }
  uint notRecoTracks = 0;
  uint notMatched = 0;
  // uint simV_id = 0;
  for (auto iTP = simV.daughterTracks_begin(); iTP != simV.daughterTracks_end(); ++iTP) {
    auto found = trackSimToRecoAssociation.find(*iTP);
    // auto found = std::find_if(trackSimToRecoAssociation.begin(), trackSimToRecoAssociation.end(), [&](const auto & oiTP) {
    //         // return tk.first == vtk;
    //         return *iTP == oiTP.key;
    //       });

    // << iTP
    if (found == trackSimToRecoAssociation.end()){
      // std::cout << "Did not find the sim track in map\n";
      // std::cout << "Sim track from sim vertex with no tracks\n";
      notRecoTracks ++;
      continue;
    } else{
      // std::cout << "Sim track from sim vertex with matched track\n";

    }

    // matched track equal to any track of reco vertex => increase counter
    for (const auto &tk : found->val) {
      if (std::find_if(recoV.tracks_begin(), recoV.tracks_end(), [&](const reco::TrackBaseRef &vtk) {
            return ((tk.first.get() != nullptr) && (vtk.get() != nullptr)) ? matchReco(*(tk.first.get()), *(vtk.get())) : false;
            //return tk.first == vtk;
          }) != recoV.tracks_end()) {
        sharedTracks += 1;
        break;
      }
      else{
        notMatched ++;
      }
    }
  }
  // std::cout << "Shared tracks sim -> reco: " << sharedTracks << "\n";
  // std::cout << "Not reco sim tracks: " << notRecoTracks << "\n";
  // std::cout << "Not matched with reco vertex: " << notMatched << "\n";
  // std::cout << "Sum: " << sharedTracks + notRecoTracks + notMatched << "\n";
  // std::cout << "Total number of sim tracks for vertex: " << simV.nDaughterTracks() << "\n";
  // std::cout << "Total number of reco tracks for vertex: " << recoV.tracksSize() << "\n";
  return sharedTracks;
}
