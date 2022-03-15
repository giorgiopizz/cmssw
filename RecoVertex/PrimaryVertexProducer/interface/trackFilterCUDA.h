#ifndef trackFilterCUDA_h
#define trackFilterCUDA_h
#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include <cstddef>
#include <cstdint>
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"

namespace trackFilterCUDA {
  struct filterParameters {
    double maxSignificance;
    double maxdxyError;
    double maxdzError;
    double minpAtIP;
    double maxetaAtIP;
    double maxchi2;
    int minpixelHits;
    int mintrackerHits;
  };
  void filterWrapper(TrackForPV::TrackForPVSoA* const itk, bool* obool, int size, filterParameters params, cudaStream_t stream);
  std::vector<reco::TransientTrack> filterWrapperOnCPU(std::vector<reco::TransientTrack> itk);
}

#endif
