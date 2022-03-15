#ifndef CUDADataFormats_Track_TrackForPVSoAT_H
#define CUDADataFormats_Track_TrackForPVSoAT_H

#include <string>
#include <algorithm>

#include "HeterogeneousCore/CUDAUtilities/interface/HistoContainer.h"

#include "CUDADataFormats/Common/interface/HeterogeneousSoA.h"
#include "HeterogeneousCore/CUDAUtilities/interface/eigenSoA.h"
template <int32_t S>
class TrackForPVSoAHeterogeneousT {
public:
  static constexpr int32_t stride() { return S; }

public:
  // Track properties needed for the PV selection + fitting

  eigenSoA::ScalarSoA<double, S> significance;
  eigenSoA::ScalarSoA<double, S> dxyError;
  eigenSoA::ScalarSoA<double, S> dzError;

  // For now, we can consider saving the full 4-momentum?
  eigenSoA::ScalarSoA<double, S> pAtIP;
  eigenSoA::ScalarSoA<double, S> etaAtIP;

  eigenSoA::ScalarSoA<double, S> chi2;

  eigenSoA::ScalarSoA<int8_t, S> nPixelHits;
  eigenSoA::ScalarSoA<int8_t, S> nTrackerHits;

};

namespace TrackForPV {

#ifdef GPU_SMALL_EVENTS
  // kept for testing and debugging
  constexpr uint32_t maxNumber() { return 2 * 1024; }
#else
  // tested on MC events with 55-75 pileup events
  constexpr uint32_t maxNumber() { return 32 * 1024; }
#endif

  using TrackForPVSoA = TrackForPVSoAHeterogeneousT<maxNumber()>;

}  
#endif
