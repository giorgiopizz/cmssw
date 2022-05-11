// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoVertex/PrimaryVertexProducer/interface/fitterCUDA.h"
#include <stdio.h>
namespace clusterizerCUDA {
__global__ void vertex(TrackForPV::TrackForFitSoA * tracks, TrackForPV::VertexForFitSoA) {
        // here we should fit tracks to a vertex (that we will have to find)
        // chi2 should be included in the vertex
        
        // linearizePointFinder
    }
}