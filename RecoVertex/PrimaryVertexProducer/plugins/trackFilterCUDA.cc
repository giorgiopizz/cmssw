// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoVertex/PrimaryVertexProducer/interface/trackFilterCUDA.h"
#include <stdio.h>
namespace trackFilterCUDA {

  __global__ void filterKernel(TrackForPV::TrackForPVSoA* const itk, bool* obool, int size, filterParameters params){
    //printf("Filter kernel %i.%i starts", threadIdx.x, blockIdx.x);
    auto &tk = *itk;
    //printf("Conversion done");
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x;
    size_t gridSize = blockDim.x * gridDim.x;
    //printf("Grid size got");
    for (int i = firstElement; i < size; i += gridSize) {
      obool[i] = false;
      if (tk.significance(i) < params.maxSignificance){
        if (tk.dxyError(i) < params.maxdxyError){
          if (tk.dzError(i) < params.maxdzError){
            if (tk.pAtIP(i) > params.minpAtIP){
              if (tk.etaAtIP(i) < params.maxetaAtIP){
                if (tk.chi2(i) < params.maxchi2){
                  if (tk.nPixelHits(i) >= params.minpixelHits){
                    if (tk.nTrackerHits(i) >= params.mintrackerHits){
                      obool[i] = true;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
 
#ifdef __CUDACC__
  // Only on GPUs, of course...
  void filterWrapper(TrackForPV::TrackForPVSoA* const itk, bool* obool, int size, filterParameters params, cudaStream_t stream){
    //std::cout << "Will filter tracks in GPU" << std::endl;
    //std::cout << "Size: "<< size << std::endl;

    int blockSize = 32;
    int gridSize  = (size + blockSize - 1)/blockSize; 
    //std::cout << gridSize << "; "<< blockSize << "; "<< size << std::endl;
    filterKernel<<<gridSize, blockSize, 0, stream>>>(itk, obool, size, params);
    cudaCheck(cudaGetLastError());
}
#else
  // Not really used, TODO: think whether to move all the CPU filter stuff from main plugin to here, makes a bit more sense
  std::vector<reco::TransientTrack> filterWrapperOnCPU(std::vector<reco::TransientTrack> itk){
    //std::cout << "Will filter tracks in CPU" << std::endl;
    return itk;
  }
#endif
}
