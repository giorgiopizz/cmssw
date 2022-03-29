#ifndef clusterizerCUDA_h
#define clusterizerCUDA_h
#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "DataFormats/Math/interface/Error.h"
#include <cstddef>
#include <cstdint>
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>

namespace clusterizerCUDA {
  struct clusterParameters {
    double Tmin;
    double Tpurge;
    double Tstop;
    double vertexSize;
    double coolingFactor;
    double d0CutOff;
    double dzCutOff;
    double uniquetrkweight;
    double uniquetrkminp;
    double zmerge;
    double sel_zrange;
    int convergence_mode;
    double delta_lowT;
    double delta_highT;
  };



  __device__ __forceinline__ void update(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* osumtkwt, double* beta, double rho0, bool updateTc){
    // printf("Update start:\n");
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x; // This is going to be the track index
    size_t gridSize = blockDim.x * gridDim.x;
    double Z_init = rho0 * exp(-(*beta) * params.dzCutOff * params.dzCutOff); // TODO::99% of the time rho0 is going to be 0, so maybe an if here can save some time. Exponentials are not cheap
    // printf("Temps set\n");

    //D if (0==threadIdx.x && 0 == blockIdx.x) printf("Update running, nt=%i, nv=%i, Zinit=%1.6f, rho0=%1.6f, updateTc=%i\n", ntracks, vertices->nTrueVertex, Z_init, rho0, updateTc);
//    clock_t start = clock();
 //   clock_t stop = clock();
    // Initiliaze stuff to 0
    for (unsigned int itrack = firstElement; itrack < ntracks ; itrack += gridSize){
      // printf("iTrack: %i:\n", itrack);

      if (not(tracks->isGood(itrack))) continue;
      for (unsigned int ivertexO = 0 ; ivertexO < vertices->nTrueVertex ; ++ivertexO){ //Only init over really existing ones
        unsigned int ivertex = vertices->order(ivertexO); // ivertex translates from ordered vertex to real vertex positions
        // printf("ivertex %i, ivertexo %i \n", ivertex, ivertexO);
        tracks->vert_sw(itrack)(ivertex) = 0.;
        // printf("--sw \n");
        tracks->vert_se(itrack)(ivertex) = 0.;
        // printf("--se \n");
        tracks->vert_swz(itrack)(ivertex) = 0.;
        // printf("--swz \n");
        if (updateTc) tracks->vert_swE(itrack)(ivertex) = 0.;
        // printf("--swE \n");
        tracks->vert_exp(itrack)(ivertex) = 0.;
        // printf("--exp \n");
        tracks->vert_exparg(itrack)(ivertex) = 0.;
        // printf("--exparg \n");
      }
    }
    //D if (0==threadIdx.x && 0==blockIdx.x) printf("Params for first vertex, before update, se=%1.10f, sw=%1.10f, swz=%1.10f, swE=%1.10f, z=%1.10f, rho=%1.10f\n", vertices->se(0), vertices->sw(0), vertices->swz(0), vertices->swE(0), vertices->z(0), vertices->rho(0));
    // printf("Everything at 0\n");
    __syncthreads();
 //   if (0==threadIdx.x && 0==blockIdx.x) printf("update stop 1: %i\n\n", (int) (clock() - stop));
 //   stop = clock();    
    // Now the monster thing about updating
    for (unsigned int itrack = firstElement; itrack < ntracks ; itrack += gridSize){
      // First, update vertex stuff
      if (not(tracks->isGood(itrack))) continue;
      double botrack_dz2 = -(*beta) * tracks->dz2(itrack);
      tracks->sum_Z(itrack) = Z_init;
      // First, let's get the partition function per track
      // printf("Track %i, kmin %i, kmax %i\n", itrack, tracks->kmin(itrack), tracks->kmax(itrack));
      for (unsigned int ivertexO = tracks->kmin(itrack) ; ivertexO < tracks->kmax(itrack) ; ++ivertexO){ // ivertexO loops over ordered vertex
        unsigned int ivertex = vertices->order(ivertexO); // ivertex translates from ordered vertex to real vertex positions
        double mult_res = tracks->z(itrack) - vertices->z(ivertex);
        // printf("Track %i, vertex %i, track_z %1.10f, track_dz2 %1.10f, track_sum_Z: %1.10f, vertex_z: %1.10f, beta: %1.10f\n", itrack, ivertexO, tracks->z(itrack), tracks->dz2(itrack), tracks->sum_Z(itrack), vertices->z(ivertexO), *beta);
        tracks->vert_exparg(itrack)(ivertex) = botrack_dz2* (mult_res * mult_res);
        tracks->vert_exp(itrack)(ivertex)    = exp(tracks->vert_exparg(itrack)(ivertex)); // exp is defined as device function in cuda
        tracks->sum_Z(itrack) += vertices->rho(ivertex)*tracks->vert_exp(itrack)(ivertex);
      }
      if(not(std::isfinite(tracks->sum_Z(itrack)))) tracks->sum_Z(itrack) = 0; // Just in case something diverges
      if(tracks->sum_Z(itrack) > 0){ // If partition > 0, then it is non-trivially assigned to a vertex and we need to compute stuff
        double sumw = tracks->weight(itrack)/tracks->sum_Z(itrack);
        for (unsigned int ivertexO = tracks->kmin(itrack) ; ivertexO < tracks->kmax(itrack) ; ++ivertexO){
          unsigned int ivertex = vertices->order(ivertexO);
          tracks->vert_se(itrack)(ivertex) = tracks->vert_exp(itrack)(ivertex) * sumw;
          double w                   = vertices->rho(ivertex) * tracks->vert_exp(itrack)(ivertex) * sumw * tracks->dz2(itrack);
          tracks->vert_sw(itrack)(ivertex)  = w;
          tracks->vert_swz(itrack)(ivertex) = w * tracks->z(itrack);
          if (updateTc) tracks->vert_swE(itrack)(ivertex) = -w * tracks->vert_exparg(itrack)(ivertex)/(*beta); // Only need it when changing the Tc
          
        }
      }
    }
    __syncthreads(); // Need to synchronize, as now we have to add across vertexes
 //   if (threadIdx.x == 0 && blockIdx.x == 0)  printf("update stop 2: %i\n\n", (int) (clock()-stop));
 //   stop = clock();
    for (unsigned int ivertexO = firstElement; ivertexO < vertices->nTrueVertex; ivertexO+=gridSize){
      unsigned int ivertex    = vertices->order(ivertexO); // ivertex translates from ordered vertex to real vertex positions
      vertices->se(ivertex)   = 0.;
      vertices->sw(ivertex)   = 0.;
      vertices->swz(ivertex)  = 0.;
      vertices->aux1(ivertex) = 0.; // Aux here is delta, the position variation in this update loop
      if (updateTc) vertices->swE(ivertex) = 0.;
    }
    __syncthreads(); //Just to be extremely careful
    
      for (unsigned int itrack = firstElement ; itrack < ntracks ; itrack+=gridSize){ //skip 0, as that is already in place
        if (not(tracks->isGood(itrack))) continue;
        for (unsigned int ivertexO = 0; ivertexO < vertices->nTrueVertex; ivertexO++){
          unsigned int ivertex    = vertices->order(ivertexO); // ivertex translates from ordered vertex to real vertex positions
        atomicAdd(&vertices->se(ivertex), tracks->vert_se(itrack)(ivertex));
        atomicAdd(&vertices->sw(ivertex) , tracks->vert_sw(itrack)(ivertex));
        atomicAdd(&vertices->swz(ivertex) , tracks->vert_swz(itrack)(ivertex));
        if (updateTc) atomicAdd(&vertices->swE(ivertex) , tracks->vert_swE(itrack)(ivertex));
      }
      }
      __syncthreads();
    for (unsigned int ivertexO = firstElement; ivertexO < vertices->nTrueVertex; ivertexO+=gridSize){
          unsigned int ivertex    = vertices->order(ivertexO); // ivertex translates from ordered vertex to real vertex positions
      if (vertices->sw(ivertex) > 0){ //The vertex position is updated
        double znew    = vertices->swz(ivertex)/vertices->sw(ivertex);
        vertices->aux1(ivertex) = abs(znew-vertices->z(ivertex));
        vertices->z(ivertex)    = znew; 
      }
      vertices->rho(ivertex)    = vertices->rho(ivertex) * vertices->se(ivertex) * (*osumtkwt);  // The relative vertex weight is updated
      // if (ivertexO==0) printf("Params for first vertex, after update, se=%1.10f, sw=%1.10f, swz=%1.10f, swE=%1.10f, z=%1.10f, rho=%1.10f\n", vertices->se(ivertexO), vertices->sw(ivertexO), vertices->swz(ivertexO), vertices->swE(ivertexO), vertices->z(ivertexO), vertices->rho(ivertexO));
    }
    __syncthreads(); //Just to be extremely careful
    // printf("Finished updating!\n");
 //   if (0==threadIdx.x && 0==blockIdx.x) printf("update stop 3: %i\n\n", (int) (clock() - stop));
 //   stop = clock();    
 //   if (0==threadIdx.x && 0==blockIdx.x) printf("update total: %i\n\n", (int) (stop - start));
    __syncthreads(); //Just to be extremely careful
  }

  __device__ __forceinline__ void set_vtx_range(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* osumtkwt, double* beta){
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x; // set_vtx_range is parallelized on tracks
    size_t gridSize = blockDim.x * gridDim.x;
    double zrange_min_= 0.1; //TODO:: put it as a param, currently hard coded as in CPU
    for (unsigned int itrack = firstElement; itrack < ntracks ; itrack+=gridSize){
      if (not(tracks->isGood(itrack))) continue;
      // printf("%i vtx_range 1\n", threadIdx.x); 
      double zrange     = std::max(params.sel_zrange/ sqrt((*beta) * tracks->dz2(itrack)), zrange_min_);
      // printf("%i vtx_range 1.1, %p\n", threadIdx.x, (void*)&zrange);
      double zmin       = tracks->z(itrack) - zrange;
      // printf("%i vtx_range 1.2, %p\n", threadIdx.x, (void*)&zmin);
      unsigned int kmin = std::min(vertices->nTrueVertex - 1,  tracks->kmin(itrack)); //We might have deleted a vertex, this might be complicated
      // printf("%i vtx_range 2, %p\n", threadIdx.x, (void*)&kmin);

      if (vertices->z(vertices->order(kmin)) > zmin){ // vertex properties always accessed through vertices->order
        while ((kmin > 0) && (vertices->z(vertices->order(std::max(kmin - 1, 0U))) > zmin)) { // i.e., while we find another vertex within range that is before the previous initial step
          kmin--;
        }
      }
      else {
        while ((kmin < (vertices->nTrueVertex - 1)) && (vertices->z(vertices->order(kmin)) < zmin)) { // Or it might happen that we have to take out vertices from the thing
          kmin++;
        }
      }
      // printf("%i vtx_range 3\n", threadIdx.x);

      // Now the same for the upper bound
      double zmax       = tracks->z(itrack) + zrange;
      unsigned int kmax = std::min(vertices->nTrueVertex - 1, tracks->kmax(itrack) - 1);
      if (vertices->z(vertices->order(kmax)) < zmax) {
        while ((kmax < (vertices->nTrueVertex  - 1)) && (vertices->z(vertices->order(kmax + 1)) < zmax)) { // As long as we have more vertex above kmax but within z range, we can add them to the collection, keep going
          kmax++;
        }
      }
      else { //Or maybe we have to restrict it
        while ((kmax > 0) && (vertices->z(vertices->order(kmax)) > zmax)) {
          kmax--;
        }
      }
      // printf("%i vtx_range 4\n", threadIdx.x);
      // Here kmin is the minimal index, kmax is the maximal
      if (kmin <= kmax) {
        tracks->kmin(itrack) = kmin;
        tracks->kmax(itrack) = kmax + 1; //always looping to tracks->kmax(i) - 1
      }
      else { // If it is here, the whole vertex are under
        tracks->kmin(itrack) = std::max(0U, std::min(kmin, kmax));
        tracks->kmax(itrack) = std::min(vertices->nTrueVertex, std::max(kmin, kmax) + 1);
      }
      // printf("%i vtx_range is finished\n", threadIdx.x);
    }
    // printf("%i device vtx_range is finished\n", threadIdx.x);
  }

   __device__ __forceinline__ void thermalize(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double * osumtkwt, double* beta, double delta_max0, double rho0){
    // We are going to be doing the same operations on all threads here
    size_t gridSize = blockDim.x * gridDim.x;
    int niter    = 0;
    double  zrange_min_ = 0.01; // Hard coded and double defined as in the CPU code, sigh....
    double delta_max = params.delta_lowT;
    ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("Start thermalizing\n");
    __syncthreads();
    if (params.convergence_mode == 0) {
      delta_max = delta_max0;
    } else if (params.convergence_mode == 1) {
      delta_max = params.delta_lowT / sqrt(std::max((*beta), 1.0));
    }
    int maxIterations_ = 1000; // TODO:: Set as Param, in the CPU version it is hard coded as well, though. Rarely goes beyond 10-20
    ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("vtx_range start\n"); 
    __syncthreads();
    set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta); // TODO::Probably want to cleanup the input a bit
    __syncthreads();
    ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("vtx_range done\n"); 
    // Accumulator of variations
    double delta_sum_range = 0;
 
    while (niter++ < maxIterations_){
      // if (threadIdx.x == 0 && blockIdx.x == 0) printf("-----Iter %i start\n", niter);
      update(ntracks, tracks, vertices, params, osumtkwt, beta, rho0, false); // Thermalizing never updates the critical T
      __syncthreads();
      ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("--------Update done\n");

      // At this stage, we have the delta per vertex in vertices->aux1
      double dmax = 0;
      for (unsigned int ivertexO = 0 ; ivertexO < vertices->nTrueVertex; ivertexO++){ // TODO::Currently we are doing this in all threads in parallel, might be optimized using shared memory?
        unsigned int ivertex = vertices->order(ivertexO);
        if (vertices->aux1(ivertex) >= dmax) dmax = vertices->aux1(ivertex);
      }
      delta_sum_range += dmax;
      __syncthreads();
      ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("--------Max delta done\n");
      if (delta_sum_range > zrange_min_) { // Check if any vertex moved a lot
        set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
        delta_sum_range = 0;
      }
      /*
        for (unsigned int ivertexO = 0 ; ivertexO < vertices->nTrueVertex; ivertexO+=gridSize){ // TODO::Currently we are doing this in all threads in parallel, might be optimized using shared memory?
          unsigned int ivertex = vertices->order(ivertexO);
          if (vertices->aux1(ivertex) > zrange_min_){ // If any moved enough, redo the track-vertex range association
            ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("--------Set vtx run for %i\n", ivertexO);
            set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
          }
        }
      }
      */
      __syncthreads();
      if (dmax < delta_max){ // at the end delta_max acts as a delta_min, below which we break the loop
        break;
      }
      ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("--------Iter done\n");
    }
    __syncthreads();
  }

  __device__ __forceinline__ void merge(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double * osumtkwt, double* beta){
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x; // This is going to be the vertex index
    size_t gridSize = blockDim.x * gridDim.x;
    unsigned int nvprev = vertices->nTrueVertex; 
    if (nvprev < 2) return; // Can't merge anything if we only have one vertex
    for (unsigned int ivertexO = firstElement; ivertexO < nvprev - 1 ; ivertexO += gridSize){
      unsigned int ivertex     = vertices->order(ivertexO);
      unsigned int ivertexnext = vertices->order(ivertexO+1);
      vertices->aux1(ivertex) = std::abs(vertices->z(ivertex)-vertices->z(ivertexnext)); // Differences between one vertex and the next
    }
    __syncthreads();

    // TODO:: MAJOR OPTIMIZATION POINT //
    // Very annoyingly, we need to sort here per thread. TODO:: Can we do this with something like thrust but better? Seems we would need to atomize the kernel significantly
    // Yes, you can't create an array of non-fixed size in cuda...
    double critical_dist[512];
    unsigned int critical_index[512];


    unsigned int ncritical = 0;
    for (unsigned int ivertexO = 0; ivertexO < nvprev - 1 ; ivertexO++){
      unsigned int ivertex     = vertices->order(ivertexO);
      if (vertices->aux1(ivertex) < params.zmerge){
        critical_dist[ncritical]  = std::fabs(vertices->aux1(ivertex));
        critical_index[ncritical] = ivertexO;
        ncritical++;
      }
    }
    if (ncritical == 0) return;
    __syncthreads();
    // Yep, this is a very bogus sorting algorithm, not even quicksort, but the size of critical shouldn't be > 10
      for (unsigned int sortO = 0; sortO < ncritical ; ++sortO){//This we might be able to parallelize more. TODO
        unsigned int ikO = 0;
        double minVal = 999999.;
        for (unsigned int sort1 = 0; sort1 < ncritical; ++sort1){
          if (critical_dist[sort1] < minVal){
            minVal = critical_dist[sort1];
            ikO    = sort1;
          }
        }
        critical_dist[ikO] = 999999.; // Out in the next loop

        unsigned int ivertexO    = critical_index[ikO];
        unsigned int ivertex     = vertices->order(ivertexO);  // This will be deleted
        unsigned int ivertexnext = vertices->order(ivertexO+1);
        vertices->isGood(ivertex) = false; // Delete it!
        //printf("removing vertex: (%d, %d)\n",ivertexO, ivertex);
        
        if (0 == threadIdx.x && 0 == blockIdx.x){ // Really no way of parallelizing this I'm afraid
            double rho =  vertices->rho(ivertex) + vertices->rho(ivertexnext);
            if (rho > 0){ 
              vertices->z(ivertexnext) = (vertices->rho(ivertex) * vertices->z(ivertex) + vertices->rho(ivertexnext) * vertices->z(ivertexnext)) / rho;
            } 
            else{
              vertices->z(ivertexnext) = 0.5 * (vertices->z(ivertex) + vertices->z(ivertexnext));
            } 
            vertices->rho(ivertexnext)  = rho;
            vertices->sw(ivertexnext)  += vertices->sw(ivertexnext);
            for (unsigned int ivertexOO = 0; ivertexOO < nvprev - 1; ++ivertexOO){ // TODO:: Any tricks here?
              if (ivertexOO >= ivertexO){ //As we copy from the next, we go forward in ivertex 
                vertices->order(ivertexOO) =vertices->order(ivertexOO+1);
              }
            }
            vertices->nTrueVertex = vertices->nTrueVertex-1; // Also update nvertex
        }
        for (unsigned int resort = 0; resort < ncritical ; ++resort){
          if (critical_index[resort] > ivertexO) critical_index[resort]--; // critical_index refers to the original vertices->order, so it needs to be updated 
        }
        for (unsigned int itrack = firstElement; itrack < ntracks ; itrack += gridSize){
          if (not tracks->isGood(itrack)) continue;
          if (tracks->kmax(itrack) > ivertexO) tracks->kmax(itrack)--;
          if ((tracks->kmin(itrack) > ivertexO) || ((tracks->kmax(itrack) < (tracks->kmin(itrack) + 1)) && (tracks->kmin(itrack) > 0))) tracks->kmin(itrack)--;
        } 
      }

      
      __syncthreads();
      set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
      __syncthreads(); 
    }
  

  __device__ __forceinline__ void split(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double * osumtkwt, double* beta, double threshold){
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x; // This is going to be the vertex index
    size_t gridSize = blockDim.x * gridDim.x;
    update(ntracks, tracks, vertices, params, osumtkwt, beta, 0.0, true); 

    double epsilon = 1e-3; //Minimum size for split
    unsigned int nvprev = vertices->nTrueVertex;
    
    // First, get Tc per vertex
    for (unsigned int ivertexO = firstElement; ivertexO < nvprev; ivertexO+=gridSize){
      unsigned int ivertex     = vertices->order(ivertexO);
      double Tc = 2 * vertices->swE(ivertex) / vertices->sw(ivertex);
      vertices->aux1(ivertex) = Tc; // Now we save the temperature here
      ////////// printf("For vertex %i, Tc %1.3f\n",ivertexO, Tc); 
    }
    __syncthreads();
    // Now, order them based on Tc, split higher Tc first
    // TODO:: MAJOR OPTIMIZATION POINT //
    // Very annoyingly, we need to sort here per thread. TODO:: Can we do this with something like thrust but better? Seems we would need to atomize the kernel significantly
    // Yes, you can't create an array of non-fixed size in cuda...
    //printf("%i allocate temps\n", threadIdx.x);
    double critical_temp[128]; // 512?
    unsigned int critical_index[128]; // 512?
    //printf("%i allocated\n", threadIdx.x);

    unsigned int ncritical = 0;

//    clock_t start = clock();
    //__shared__ double p1_v[512], z1_v[512], w1_v[512], p2_v[512], z2_v[512], w2_v[512];
    __shared__ double p1, z1, w1, p2, z2, w2;

    for (unsigned int ivertexO = 0; ivertexO < nvprev  ; ivertexO++){
      unsigned int ivertex     = vertices->order(ivertexO);
      ////////// if (threadIdx.x == 2 && blockIdx.x == 0) printf("Vertex %i Tc %1.10f beta %1.10f threshold %1.10f \n",ivertexO, vertices->aux1(ivertex), *beta, threshold);
      if (vertices->aux1(ivertex)*(*beta) > threshold){
        ////////// if (threadIdx.x == 0 && blockIdx.x == 0) printf("I'll split vertex %i\n",ivertexO);
        //printf("%i sees %i needs to split\n", threadIdx.x, ivertexO);
        critical_temp[ncritical]  = std::fabs(vertices->aux1(ivertex));
        critical_index[ncritical] = ivertexO;
        ncritical++;
        if (ncritical == 128) break;
        //printf("%i sees %i splitted\n", threadIdx.x, ivertexO);
      }
    }
    __syncthreads();
 //   clock_t stop = clock();
 //   if (threadIdx.x == 0 && blockIdx.x == 0) printf("Split 1st loop: %i\n\n",(int) (stop-start));
    //printf("%i Verifies splitting\n", threadIdx.x);
    if (ncritical == 0) return;
    //printf("%i Splitting verified\n", threadIdx.x);
    __syncthreads();
    // Yep, this is a very bogus sorting algorithm, not even quicksort, but the size of critical shouldn't be > 10

      for (unsigned int sortO = 0; sortO < ncritical ; ++sortO){//This we might be able to parallelize more. TODO
        //start = clock();
        //stop = clock();
        unsigned int ikO = 0;
        double maxVal = -1.;
        for (unsigned int sort1 = 0; sort1 < ncritical; ++sort1){
          if (critical_temp[sort1] > maxVal){
            maxVal = critical_temp[sort1];
            ikO    = sort1;
          }
        }
        critical_temp[ikO] = -1.; // Out in the next loop
        unsigned int ivertexO    = critical_index[ikO];
        unsigned int ivertex     = vertices->order(ivertexO);  // This will be splitted
        unsigned int ivertexprev = 0;
        unsigned int ivertexnext = nvprev -1; 
        // A little bit of safety here. First is needed to avoid reading the -1 entry of vertices->order. Second is not as far as we don't go over 511 vertices, but better keep it just in case
        if (ivertexO > 0) ivertexprev = vertices->order(ivertexO-1);  // This will be used in a couple of computations
        if (ivertexO < nvprev -1) ivertexnext = vertices->order(ivertexO+1);  // This will be used in a couple of computations

//        if (threadIdx.x == 0 && blockIdx.x == 0) printf("Split loop of 1.-1 event: %i\n\n",(int) (clock()-stop));
//        stop = clock();
        if (threadIdx.x == 0 && blockIdx.x == 0) {
            p1 = 0;
            z1 = 0;
            w1 = 0;
            p2 = 0;
            z2 = 0;
            w2 = 0;
        }
        __syncthreads();
        for (unsigned int itrack = firstElement; itrack < ntracks; itrack+=gridSize) {
          if (not(tracks->isGood(itrack))) continue;
          if (tracks->sum_Z(itrack) > 1.e-100) {
            // winner-takes-all, usually overestimates splitting
            double tl = tracks->z(itrack) < vertices->z(ivertex) ? 1. : 0.;
            double tr = 1. - tl;
            // soften it, especially at low T
            double arg = (tracks->z(itrack) - vertices->z(ivertex)) * sqrt((*beta) * tracks->dz2(itrack));
            if (std::fabs(arg) < 20) {
              double t = exp(-arg);
              tl = t / (t + 1.);
              tr = 1 / (t + 1.);
            }
            double p = vertices->rho(ivertex) * tracks->weight(itrack) * exp(-(*beta) * (tracks->z(itrack)-vertices->z(ivertex))*(tracks->z(itrack)-vertices->z(ivertex))* tracks->dz2(itrack))/ tracks->sum_Z(itrack);
            double w = p * tracks->dz2(itrack);
            atomicAdd(&p1, p*tl);
            atomicAdd(&z1, w*tl*tracks->z(itrack));
            atomicAdd(&w1, w*tl);
            atomicAdd(&p2, p*tr);
            atomicAdd(&z2, w*tr*tracks->z(itrack));
            atomicAdd(&w2, w*tr);
          }
        }
        /*
        for (unsigned int itrack = firstElement; itrack < ntracks; itrack+=gridSize) {
          if (not(tracks->isGood(itrack))) continue;
          if (tracks->sum_Z(itrack) > 1.e-100) {
            // winner-takes-all, usually overestimates splitting
            double tl = tracks->z(itrack) < vertices->z(ivertex) ? 1. : 0.;
            double tr = 1. - tl;
            // soften it, especially at low T
            double arg = (tracks->z(itrack) - vertices->z(ivertex)) * sqrt((*beta) * tracks->dz2(itrack));
            if (std::fabs(arg) < 20) {
              double t = exp(-arg);
              tl = t / (t + 1.);
              tr = 1 / (t + 1.);
            }
            double p = vertices->rho(ivertex) * tracks->weight(itrack) * exp(-(*beta) * (tracks->z(itrack)-vertices->z(ivertex))*(tracks->z(itrack)-vertices->z(ivertex))* tracks->dz2(itrack))/ tracks->sum_Z(itrack);
            double w = p * tracks->dz2(itrack);
            if (itrack == firstElement) {
                p1_v[firstElement] = p*tl;
                z1_v[firstElement] =  w*tl*tracks->z(itrack);
                w1_v[firstElement] =  w*tl;
                p2_v[firstElement] =  p*tr;
                z2_v[firstElement] =  w*tr*tracks->z(itrack);
                w2_v[firstElement] =  w*tr;
            }
            else {
                p1_v[firstElement] += p*tl;
                z1_v[firstElement] +=  w*tl*tracks->z(itrack);
                w1_v[firstElement] +=  w*tl;
                p2_v[firstElement] +=  p*tr;
                z2_v[firstElement] +=  w*tr*tracks->z(itrack);
                w2_v[firstElement] +=  w*tr;
            }
          }
        }
        */
        __syncthreads();
        if (threadIdx.x == 0 && blockIdx.x == 0){
            /*
            p1 = 0;
            z1 = 0;
            w1 = 0;
            p2 = 0;
            z2 = 0;
            w2 = 0;
            
            for (unsigned int i = 0; i<gridSize; i++) {
               p1 += p1_v[i];      
               z1 += z1_v[i];      
               w1 += w1_v[i];      
               p2 += p2_v[i];      
               z2 += z2_v[i];      
               w2 += w2_v[i];      
            }
            */
            
            //printf("Split loop of 1.0 event: %i\n\n",(int) (clock()-stop));
            // stop = clock();
            if (w1 > 0) {
              z1 = z1 / w1;
            }
            else {
              z1 = vertices->z(ivertex) - epsilon;
            }
            if (w2 > 0) {
              z2 = z2 / w2;
            }  
            else {
              z2 = vertices->z(ivertex) + epsilon;
            }
            // reduce split size if there is not enough room
            if ((ivertexO > 0) && (z1 < (0.6 * vertices->z(ivertex) + 0.4 * vertices->z(ivertexprev)))) { // First in the if is ivertexO, as we care on whether the vertex is the leftmost or rightmost
              z1 = 0.6 * vertices->z(ivertex) + 0.4 * vertices->z(ivertexprev);
            }
            if ((ivertexO + 1 < nvprev) && (z2 > (0.6 * vertices->z(ivertex) + 0.4 * vertices->z(ivertexnext)))) {
              z2 = 0.6 * vertices->z(ivertex) + 0.4 * vertices->z(ivertexnext);
            }
        }
        __syncthreads();
        ////////// printf("New Z1, Z2, epsilon: %1.10f %1.10f %1.10f\n", z1, z2, epsilon);
        if (abs(z2-z1) > epsilon){
          //if (z2<=z1) printf("\n\nAttention: z2 <= z1: %f %f\n\n", z2, z1);
          // This we can for sure parallelize, but we need to be careful enclosing stuff on the if threadIdx.x stuff 
          for (unsigned int itrack = firstElement; itrack < ntracks ; itrack+=gridSize){
            if (not(tracks->isGood(itrack))) continue;
            if (tracks->kmin(itrack) > ivertexO) tracks->kmin(itrack)++;
            if ((tracks->kmax(itrack) >= ivertexO) || (tracks->kmax(itrack) == tracks->kmin(itrack))) tracks->kmax(itrack)++;
          }
          if (threadIdx.x == 0 && blockIdx.x == 0){
              //printf("Split loop of 1b event: %i\n\n",(int) (clock()-stop));
              //stop = clock();

              double pk1 = p1 * vertices->rho(ivertex) / (p1 + p2);
              double pk2 = p2 * vertices->rho(ivertex) / (p1 + p2);
              vertices->z(ivertex) = z2;
              vertices->rho(ivertex)  = pk2;
              // Now we need to get the first empty index to save the vertex in
              unsigned int nnew = 999999;
              unsigned int icheck = 0;
              while (nnew == 999999){
                if (not(vertices->isGood(icheck))) nnew = icheck;
                icheck++;
              }
              //printf("found nnew: %d, while ntruevertices: %d\n", nnew, nvprev);
              //nnew = vertices->order(nnew);
              // Insert it into the first available slot           
              vertices->z(nnew)      = z1; 
              vertices->rho(nnew)    = pk1; 
              // And register it as used
              vertices->isGood(nnew) = true;
              // TODO:: this is likely not needed as far as it is reset anytime we call update
              vertices->sw(nnew)     = 0.;
              vertices->se(nnew)     = 0.;
              vertices->swz(nnew)    = 0.;
              vertices->swE(nnew)    = 0.;
              vertices->exp(nnew)    = 0.;
              vertices->exparg(nnew) = 0.;
              //printf("Split loop of 1a event: %i\n\n",(int) (clock()-stop));
              //stop = clock();
              for (unsigned int ivnew = nvprev ; ivnew > ivertexO ; ivnew--){ // As we add a vertex, we update from the back downwards

                ////////// printf("I'm changing order %i for %i at %i\n", vertices->order(ivnew), vertices->order(ivnew-1), ivnew);
                vertices->order(ivnew) = vertices->order(ivnew-1);

              }
              //printf("Split loop of 1c event: %i\n\n",(int) (clock()-stop));
              //stop = clock();
              // And then the new vertex will be at the old split one position
              ////////// printf("I'm adding %i at %i\n", nnew, ivertexO);
              vertices->order(ivertexO) = nnew;
              vertices->nTrueVertex++; // Add one to the count of real vertex
              //stop = clock();
              //printf("Split loop of 1 event: %i\n\n",(int) (stop-start));
          }
          __syncthreads();
          nvprev++; // And to the counter of previous vertices
          // However, the original list of vertices to be split is going to refer to the old vertices_order vector, so we need to add some shenanigans
          for (unsigned int resort = 0; resort < ncritical ; ++resort){
            if (critical_index[resort] > ivertexO) critical_index[resort]++; // critical_index refers to the original vertices->order, so it needs to be updated 
          }
        ////////// printf("Vertex created\n");
      }
    }
   // __syncthreads();
  }


  __device__ __forceinline__ void purge(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double * osumtkwt, double* beta, double rho0){
    // If only one vertex, there is nothing to purge
    if (vertices->nTrueVertex < 2) return;
    constexpr double eps = 1.e-100;
    constexpr int nunique_min = 2; // Hardcoded as in the CPU version. Why?
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x;
    size_t gridSize = blockDim.x * gridDim.x;
    unsigned int nvprev = vertices->nTrueVertex;
    double rhoconst = rho0*exp(-(*beta)*params.dzCutOff*params.dzCutOff); // Damping term in the purging parameter
    // First, set up the track-vertex assignments again
    set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta); // TODO::Probably want to cleanup the input a bit
    __syncthreads();
    for (unsigned int itrack = firstElement; itrack < ntracks; itrack +=gridSize){
      if (not(tracks->isGood(itrack))) continue;
      tracks->aux1(itrack) = ((tracks->sum_Z(itrack) > eps) && (tracks->weight(itrack) > params.uniquetrkminp)) ? 1. / tracks->sum_Z(itrack) : 0.; //invZ
      // WARNING: REUSING  matrices from tracks, but no relation to this quantities whatsoever
      // TODO::Make this more explicit somehow?
      for (unsigned int ivertexO = 0 ; ivertexO < nvprev; ivertexO++){
        unsigned int ivertex = vertices->order(ivertexO);
        if (ivertexO >= tracks->kmin(itrack) && ivertexO < tracks->kmax(itrack)){
          tracks->vert_exp(itrack)(ivertex) = exp(-(*beta)*tracks->dz2(itrack) * ( (tracks->z(itrack)-vertices->z(ivertex))*(tracks->z(itrack)-vertices->z(ivertex)) )); // TODO: either this or std::pow?
        }
        else{
          tracks->vert_exp(itrack)(ivertex) = 0;
        }
      }
    }
    __syncthreads();
    for (unsigned int ivertexO = firstElement ; ivertexO < nvprev;  ivertexO += gridSize){
      unsigned int ivertex = vertices->order(ivertexO);
      double ppcut = params.uniquetrkweight * vertices->rho(ivertex) / (vertices->rho(ivertex)+rhoconst); 
      vertices->aux1(ivertex) = 0; //psump
      vertices->aux2(ivertex) = 0; //pnUnique
      for (unsigned int itrack = 0; itrack < ntracks; itrack++){
        if (not(tracks->isGood(itrack))) continue;
        double p = vertices->rho(ivertex)*tracks->vert_exp(itrack)(ivertex)*tracks->aux1(itrack);
        vertices->aux1(ivertex) += p; //psump
        vertices->aux2(ivertex) += (p > ppcut) ? 1 : 0; //pnUique
      }
    }
    __syncthreads();
    if (0==threadIdx.x && 0==blockIdx.x){
      double sumpmin  = ntracks;
      unsigned int k0 = nvprev;
      for (unsigned int ivertexO = 0; ivertexO < nvprev ; ivertexO++){
        unsigned int ivertex = vertices->order(ivertexO);
        if ((vertices->aux2(ivertex) < nunique_min) && (vertices->aux1(ivertex) < sumpmin)){
          // Will purge the worst one
          sumpmin = vertices->aux1(ivertex);
          k0 = ivertexO;
        }
      }
      if (k0 != nvprev){
        for (unsigned int ivertexOO = 0; ivertexOO < nvprev - 1; ++ivertexOO){ // TODO:: Any tricks here?
          if (ivertexOO >= k0){ //As we copy from the next, we go forward in ivertex 
            vertices->order(ivertexOO) =vertices->order(ivertexOO+1);
          }
        }
        vertices->nTrueVertex = vertices->nTrueVertex-1; // Also update nvertex
        for (unsigned int itrack = firstElement; itrack < ntracks ; itrack += gridSize){
          if (not tracks->isGood(itrack)) continue;
          if (tracks->kmax(itrack) > k0) tracks->kmax(itrack)--;
          if ((tracks->kmin(itrack) > k0) || ((tracks->kmax(itrack) < (tracks->kmin(itrack) + 1)) && (tracks->kmin(itrack) > 0))) tracks->kmin(itrack)--;
        }
      }
    }
    __syncthreads();
    if (nvprev != vertices->nTrueVertex){
      set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
    }
    __syncthreads();
   
  } 
  __device__ __forceinline__ void checkOrder (int n, TrackForPV::VertexForPVSoA* vertices){
  if (threadIdx.x == 0 && blockIdx.x == 0){
    printf("Begin checkOrder %d \n",n);
    for (unsigned int ivertexO =0; ivertexO<vertices->nTrueVertex-1; ivertexO++){
      unsigned int ivertexActual = vertices->order(ivertexO);
      unsigned int ivertexNext = vertices->order(ivertexO+1);
      if (vertices->z(ivertexActual)>=vertices->z(ivertexNext)){
        printf("Error for vertex (%d, %d, %f) and next vertex (%d, %d, %f)\n", ivertexO, ivertexActual, vertices->z(ivertexActual) , ivertexO+1, ivertexNext, vertices->z(ivertexNext));
        break;
      }
      else {  
        printf("Everything ok for vertex (%d, %d, %f) and next vertex (%d, %d, %f)\n", ivertexO, ivertexActual, vertices->z(ivertexActual) , ivertexO+1, ivertexNext, vertices->z(ivertexNext));
      }
      
    }

  }
}
  std::vector<TransientVertex> vertices(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, const std::vector<reco::TransientTrack>& t_tks, double * beta) ;
  std::vector<std::vector<reco::TransientTrack>> clusterize(std::vector<TransientVertex>& pv, clusterParameters params) ;
  void initializeWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
  void getBeta0Wrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
  void thermalizeWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
  void coolingWhileSplittingWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
  void remergeTracksWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
  void resplitTracksWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
  void outlierRejectionWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream);
}
#endif
