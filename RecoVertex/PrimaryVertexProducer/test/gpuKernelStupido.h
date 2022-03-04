#include <cmath>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>


#include <cstddef>
#include <cstdint>

//#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
/*
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "HeterogeneousCore/CUDAUtilities/interface/requireDevices.h"
#include "HeterogeneousCore/CUDAUtilities/interface/launch.h"
*/


/*
struct Event {
  std::vector<float> zvert;
  std::vector<uint16_t> itrack;
  std::vector<float> ztrack;
  std::vector<float> eztrack;
  std::vector<float> pttrack;
  std::vector<uint16_t> ivert;
};

//__device__ void kernel1(
*/
/*

struct vertex_t {
    double z;   //           z coordinate
    double pk;  //           vertex weight for "constrained" clustering
    // --- temporary numbers, used during update
    double ei;
    double sw;
    double swz;
    double se;
    // ---for Tc
    double swE;
    double Tc;
  };


struct Event {
  std::vector<float> zvert;
  std::vector<uint16_t> itrack;
  std::vector<float> ztrack;
  std::vector<float> eztrack;
  std::vector<float> pttrack;
  std::vector<uint16_t> ivert;
};
*/
//using namespace std;





/*
vector<track_t> fill(const const vector<reco::TransientTrack>& tracks) const {
  // prepare track data for clustering
  vector<track_t> tks;
  for (vector<reco::TransientTrack>::const_iterator it = tracks.begin(); it != tracks.end(); it++) {
    track_t t;
    t.z = ((*it).stateAtBeamLine().trackStateAtPCA()).position().z();
    double tantheta = tan(((*it).stateAtBeamLine().trackStateAtPCA()).momentum().theta());
    double phi = ((*it).stateAtBeamLine().trackStateAtPCA()).momentum().phi();
    //  get the beam-spot
    reco::BeamSpot beamspot = (it->stateAtBeamLine()).beamSpot();
    t.dz2 = pow((*it).track().dzError(), 2)  // track errror
            + (pow(beamspot.BeamWidthX() * cos(phi), 2) + pow(beamspot.BeamWidthY() * sin(phi), 2)) /
                  pow(tantheta, 2)  // beam-width induced
            + pow(vertexSize_, 2);  // intrinsic vertex size, safer for outliers and short lived decays
    if (d0CutOff_ > 0) {
      Measurement1D IP = (*it).stateAtBeamLine().transverseImpactParameter();       // error constains beamspot
      t.pi = 1. / (1. + exp(pow(IP.value() / IP.error(), 2) - pow(d0CutOff_, 2)));  // reduce weight for high ip tracks
    } else {
      t.pi = 1.;
    }
    t.tt = &(*it);
    t.Z = 1.;
    tks.push_back(t);
  }
  return tks;
}


*/


#define N   1024
#define RADIUS 3
#define BLOCK_SIZE 16

namespace gpuKernelStupido {
    
    static constexpr int32_t MAXTRACKS = 1000;

    struct track_SoA {
        double z[MAXTRACKS];                        // z-coordinate at point of closest approach to the beamline
        double dz2[MAXTRACKS];                      // square of the error of z(pca)
        //const reco::TransientTrack *tt;  // a pointer to the Transient Track
        double Z[MAXTRACKS];                        // Z[i] for DA clustering
        double pi[MAXTRACKS];                       // track weight
    };
    // void fill_ints(int * x, int n);
    // __global__ void myKernel(int *in, int *out);

    
    void makeAsync(track_SoA tks_SoA);
/*
        void makeAsync(cudaStream_t stream, track_SoA tks_SoA) const{
            std::cout << "Ciao" << std::endl;         
            int *in, * out;
            int *d_in, *d_out;
            int size = (N+2*RADIUS) * sizeof(int);
            
            in = (int * ) malloc(size); fill_ints(in, N+ 2*RADIUS);
            out = (int *) malloc(size); fill_ints(out, N+ 2*RADIUS);
            
            cudaMalloc((void **) &d_in, size);
            cudaMalloc((void **) &d_out, size);
            
            cudaMemcpy(d_in, in, size, cudaMemcpyHostToDevice);
            cudaMemcpy(d_out, out, size, cudaMemcpyHostToDevice);

            myKernel<<<N/BLOCK_SIZE, BLOCK_SIZE>>>(d_in + RADIUS, d_out + RADIUS);

            cudaMemcpy(out, d_out,size, cudaMemcpyDeviceToHost);
            
            free(in); free(out);
            cudaFree(d_in); cudaFree(d_out);
        }
*/            

}
/*
        int *in, * out;
        int *d_in, *d_out;
        int size = (N+2*RADIUS) * sizeof(int);
        
        in = (int * ) malloc(size); fill_ints(in, N+ 2*RADIUS);
        out = (int *) malloc(size); fill_ints(out, N+ 2*RADIUS);
        
        cudaMalloc((void **) &d_in, size);
        cudaMalloc((void **) &d_out, size);
        
        cudaMemcpy(d_in, in, size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_out, out, size, cudaMemcpyHostToDevice);

        myKernel<<<N/BLOCK_SIZE, BLOCK_SIZE>>>(d_in + RADIUS, d_out + RADIUS);

        cudaMemcpy(out, d_out,size, cudaMemcpyDeviceToHost);
        
        free(in); free(out);
        cudaFree(d_in); cudaFree(d_out);
*/
        
