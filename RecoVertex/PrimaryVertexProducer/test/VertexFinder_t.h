nclude <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/requireDevices.h"
#include "HeterogeneousCore/CUDAUtilities/interface/launch.h"
#ifdef USE_DBSCAN
#include "RecoPixelVertexing/PixelVertexFinding/plugins/gpuClusterTracksDBSCAN.h"
#define CLUSTERIZE gpuVertexFinder::clusterTracksDBSCAN
#elif USE_ITERATIVE
#include "RecoPixelVertexing/PixelVertexFinding/plugins/gpuClusterTracksIterative.h"
#define CLUSTERIZE gpuVertexFinder::clusterTracksIterative
#else
#include "RecoPixelVertexing/PixelVertexFinding/plugins/gpuClusterTracksByDensity.h"
#define CLUSTERIZE gpuVertexFinder::clusterTracksByDensityKernel
#endif
#include "RecoPixelVertexing/PixelVertexFinding/plugins/gpuFitVertices.h"
#include "RecoPixelVertexing/PixelVertexFinding/plugins/gpuSortByPt2.h"
#include "RecoPixelVertexing/PixelVertexFinding/plugins/gpuSplitVertices.h"



