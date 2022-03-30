// CUDA include files
#include <cuda_runtime.h>

// CMSSW include files
#include "HeterogeneousCore/CUDAUtilities/interface/cudaCheck.h"
#include "RecoVertex/PrimaryVertexProducer/interface/clusterizerCUDA.h"
#include <stdio.h>
namespace clusterizerCUDA {



__global__ void initializeKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params){
    // TODO:: Study performance of explicitly writing squares vs std::pow. Study performance of atomics for the serial part
    size_t firstElement = threadIdx.x + blockIdx.x * blockDim.x;
    size_t gridSize = blockDim.x * gridDim.x;
    // First the preselection block
    vertices->nTrueVertex = 1; // Will add one here. As it is hard assignment, not really care to separate by threads, a bit lazy it is
    for (unsigned int i = firstElement; i < vertices->stride(); i += gridSize) {
      vertices->sw(i)     = 0.;
      vertices->se(i)     = 0.;
      vertices->swz(i)    = 0.;
      vertices->swE(i)    = 0.;
      vertices->exp(i)    = 0.;
      vertices->exparg(i) = 0.;
      vertices->z(i)      = 0.;
      vertices->rho(i)    = 0.;
      vertices->isGood(i) = false;
      // The order vector is special, -1 means undefined
      vertices->order(i)  = -1;
      if (i == 0){ // TODO:: Cleanup
        vertices->z(i)     = 0.;
        vertices->rho(i)   = 1.; // First vertex has all the swag
        vertices->order(i) = 0 ;  // And it is the only one to loop over
        vertices->isGood(i)= true;
      }
    }
    __syncthreads(); // Synchronize after loop,just in case
}

__global__ void getBeta0Kernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* beta){
    double znew = 0;
    double wnew = 0;
    size_t firstElement = threadIdx.x + blockIdx.x*blockDim.x; // This is going to be the track index
    size_t gridSize = blockDim.x*gridDim.x;
    for (unsigned int itrack = firstElement; itrack < ntracks ; itrack += gridSize){
      if (not(tracks->isGood(itrack))) continue;
      //printf("Track %i is Good, adding stuff\n", itrack);
      tracks->aux1(itrack)  = tracks->weight(itrack)*tracks->dz2(itrack); // Will be sumw
      tracks->aux2(itrack)  = tracks->weight(itrack)*tracks->dz2(itrack)*tracks->z(itrack); // Will be sumwz
    }
    __syncthreads();
    if (0 == threadIdx.x && 0 == blockIdx.x){ // Serialized code. TODO:: Test atomicAdd instead, or better even a syncable sum 
      for (unsigned int itrack = 0; itrack < ntracks ; itrack++){
        if (not(tracks->isGood(itrack))) continue;
        znew += tracks->aux2(itrack); //sumwz
        wnew += tracks->aux1(itrack); //sumw
        //printf("znew, wnew, itrack: %1.16f, %1.16f, %i \n", znew, wnew, itrack);
      }
      vertices->z(0) = znew/wnew; //New z is the quotient 
    }
    __syncthreads();
    for (unsigned int itrack = firstElement; itrack < ntracks ; itrack += gridSize){
      if (not(tracks->isGood(itrack))) continue;
      tracks->aux2(itrack) = tracks->aux1(itrack)*(vertices->z(0) - tracks->z(itrack) )*(vertices->z(0) - tracks->z(itrack))*tracks->dz2(itrack);
    }
    __syncthreads();
    if (0 == threadIdx.x && 0 == blockIdx.x){ // More serialized code. TODO:: Test atomicAdd instead, or better even a syncable sum 
      double a = 0;
      for (unsigned int itrack = 0; itrack < ntracks ; itrack++){
        if (not(tracks->isGood(itrack))) continue;
        a += tracks->aux2(itrack);
      }
      (*beta) = 2 * a/wnew; // Here it is technically 1/beta0, i.e. Tc, but ok to save on memory allocation
      ////////// printf("Beta0 before comparing: %1.3f\n", (*beta));
      ////////// printf("znew, wnew, a: %1.16f, %1.16f, %1.16f \n", znew, wnew, a);

      double betamax_ = 1./params.Tmin; //From the config file, 
      ////////// printf("Betamax, coolingFactor %1.10f, %1.10f\n", betamax_, params.coolingFactor);
      if ((*beta) > 1./betamax_){
        int coolingsteps = 1 - int(std::log((*beta) * betamax_) / std::log(params.coolingFactor)); // A tricky conversion to round the number
        (*beta) = betamax_ * std::pow(params.coolingFactor, coolingsteps);
        ////////// printf("Betamax, coolingFactor, coolingsteps %1.10f, %1.10f %i \n", betamax_, params.coolingFactor, coolingsteps);
      }
      else{
        (*beta) = betamax_ * params.coolingFactor;
      }
      //TODO::Add debugging option
      ////////// printf("1./Beta0 %1.3f\n", 1./(*beta));
      //printf("%1.10f",beta[0]);
    }
    __syncthreads();
}

__global__ void thermalizeKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double * osumtkwt, double* beta, double delta_max0, double rho0){
   thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, delta_max0, rho0);
   //Dif (0 == threadIdx.x && 0 == blockIdx.x){
   //D  for (unsigned int ivertexO = 0; ivertexO < vertices->nTrueVertex ; ivertexO++){
   //D    printf("At T=%1.3f, ivertex=%i, pos=%1.10f, rho=%1.10f", 1./(*beta), ivertexO, vertices->z(vertices->order(ivertexO)), vertices->rho(vertices->order(ivertexO)));
   //D  }
   //D}
   __syncthreads();
}

__global__ void coolingWhileSplittingKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* osumtkwt, double* beta){
  double betafreeze = (1./params.Tmin) * sqrt(params.coolingFactor); //Last T to be updated
  // First the T loop
  ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("Start cooling! \n");
  while ((*beta) < betafreeze) {
    unsigned int nprev = vertices->nTrueVertex;
    ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("New T: %1.5f ; nv = %i \n",1./(*beta), nprev);
    
//    clock_t mstart = clock();
    merge(ntracks, tracks, vertices, params, osumtkwt, beta);
 //   clock_t mstop = clock();
 //   if (threadIdx.x == 0 && 0 == blockIdx.x) printf("Clock for merge: %i\n", (int) (mstop-mstart));
    
    //checkOrder(0, vertices);
    __syncthreads();
    
    //merge(ntracks, tracks, vertices, params, osumtkwt, beta);
    //__syncthreads();
    //checkOrder(1, vertices);
    ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("After merging nv = %i \n", vertices->nTrueVertex);

    while (nprev !=  vertices->nTrueVertex) { // While merge is true
      nprev = vertices->nTrueVertex;
      __syncthreads();
      
 //     clock_t ustart = clock();
      update(ntracks, tracks, vertices, params, osumtkwt, beta, 0.0, false); //Udpdate them 
 //     clock_t ustop = clock();
 //     if (threadIdx.x == 0 && 0 == blockIdx.x) printf("Clock for update: %i\n", (int) (ustop-ustart));
      __syncthreads();
 //     mstart = clock();
      merge(ntracks, tracks, vertices, params, osumtkwt, beta);
 //     mstop = clock();
 //     if (threadIdx.x == 0 && 0 == blockIdx.x) printf("Clock for merge: %i\n", (int) (mstop-mstart));
      
      //update(ntracks, tracks, vertices, params, osumtkwt, beta, 0.0, false); //Udpdate them 

      ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("After merging nv = %i \n", vertices->nTrueVertex);
      //
      __syncthreads();
      //checkOrder(2, vertices);
    }
    //checkOrder(21, vertices);
    ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("After merge loop nv = %i \n", vertices->nTrueVertex);
    
 //   clock_t sstart = clock();
    split(ntracks, tracks, vertices, params, osumtkwt, beta, 1.); // Then split if we need to
 //   clock_t sstop = clock();
 //   if (threadIdx.x == 0 && 0 == blockIdx.x) printf("Clock for split: %i\n", (int) (sstop-sstart));
    
    //split(ntracks, tracks, vertices, params, osumtkwt, beta, 1.); // Then split if we need to
    __syncthreads();
    //checkOrder(3, vertices);
    ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("After splitting nv = %i \n", vertices->nTrueVertex);

    if (0 == threadIdx.x && 0 == blockIdx.x) (*beta) = (*beta) / params.coolingFactor; // Reduce temperature
    ////////// if (0 == threadIdx.x && 0 == blockIdx.x) printf("New T = %1.5f \n", 1./(*beta));

    __syncthreads();
 //   clock_t tstart = clock();
    thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_highT, 0.0); // And recompute everything at the new temperature
 //   clock_t tstop = clock();
 //   if (threadIdx.x == 0 && 0 == blockIdx.x) printf("Clock for thermalize: %i\n", (int) (tstop-tstart));
    __syncthreads();
    //checkOrder(4, vertices);
  }
  // After the T loop, reassign vertices, and update again
 // clock_t srstart = clock();
  set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
 // clock_t srstop = clock();
 // if (threadIdx.x == 0 && 0 == blockIdx.x) printf("Clock for set vtx range: %i\n", (int) (srstop-srstart));
  __syncthreads();
  update(ntracks, tracks, vertices, params, osumtkwt, beta, 0.0, false);
  //checkOrder(5, vertices);
  __syncthreads();
}

__global__ void remergeTracksKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* osumtkwt, double* beta){
  unsigned int nprev = vertices->nTrueVertex;
  merge(ntracks, tracks, vertices, params, osumtkwt, beta);
  __syncthreads();
  while (nprev !=  vertices->nTrueVertex) {
    set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
    __syncthreads();
    update(ntracks, tracks, vertices, params, osumtkwt, beta, 0.0, false);
    __syncthreads();
    nprev = vertices->nTrueVertex;   
    merge(ntracks, tracks, vertices, params, osumtkwt, beta);
    __syncthreads();
  }
}

__global__ void resplitTracksKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* osumtkwt, double* beta){
  unsigned int ntry = 0; 
  double threshold = 1.0;
  unsigned int nprev = vertices->nTrueVertex;
  split(ntracks, tracks, vertices, params, osumtkwt, beta, threshold);
  __syncthreads();
  while (nprev !=  vertices->nTrueVertex && (ntry++<10)) {
    thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_highT, 0.0); // if split, we recompute everything 
    __syncthreads();
    // We might need to merge afterwards!
    nprev = vertices->nTrueVertex;
    __syncthreads();
    merge(ntracks, tracks, vertices, params, osumtkwt, beta);
    __syncthreads();
    while (nprev !=  vertices->nTrueVertex) {
      nprev = vertices->nTrueVertex;
      update(ntracks, tracks, vertices, params, osumtkwt, beta, 0.0, false); //Udpdate them 
      __syncthreads();
      merge(ntracks, tracks, vertices, params, osumtkwt, beta);
      __syncthreads();
    }
    nprev = vertices->nTrueVertex;
    threshold *= 1.1;
    split(ntracks, tracks, vertices, params, osumtkwt, beta, threshold);
    __syncthreads();
  }
}


__global__ void outlierRejectionKernel(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, double* osumtkwt, double* beta){
  double rho0 = 0.0;
  if (params.dzCutOff > 0){
    rho0 = vertices->nTrueVertex > 1 ? 1./vertices->nTrueVertex : 1.;
    for (unsigned int a = 0; a < 5 ; a++){ //Can't be parallelized in any reasonable way
      update(ntracks, tracks, vertices, params, osumtkwt, beta, a*rho0/5., false);
      __syncthreads();
    }
  }
  thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_lowT, rho0); // If we don't do outlier rejection, we just thermalize afther the resplitting step at either rho0=0 or rho0 = 1/nv
  __syncthreads();
  // With outlier rejection we might lose tracks in some vertices, so merge again
  unsigned int nprev = vertices->nTrueVertex;
  __syncthreads();
  merge(ntracks, tracks, vertices, params, osumtkwt, beta);
  __syncthreads();
  while (nprev !=  vertices->nTrueVertex) {
    set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
    __syncthreads();
    update(ntracks, tracks, vertices, params, osumtkwt, beta, rho0, false); // Now at the proper rho0! This changes Z_init and so the whole partition function for tracks
    __syncthreads();
    nprev = vertices->nTrueVertex;
    __syncthreads();
    merge(ntracks, tracks, vertices, params, osumtkwt, beta);
    __syncthreads();
  }
  // Now we go to the purge temperature, without splitting or merging
  double betapurge = 1./params.Tpurge;
  while ((*beta) < betapurge){
    if (0==threadIdx.x && 0 == blockIdx.x)  (*beta) = std::min((*beta)/params.coolingFactor, betapurge);
    __syncthreads();
    thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_lowT, rho0);
    __syncthreads();
  }
  // Now it is time for the purging
  nprev = vertices->nTrueVertex;
  purge(ntracks, tracks, vertices, params, osumtkwt, beta, rho0);
  __syncthreads();
  while (nprev !=  vertices->nTrueVertex) {
    thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_lowT, rho0);
    __syncthreads();
    nprev = vertices->nTrueVertex;
    __syncthreads();
    purge(ntracks, tracks, vertices, params, osumtkwt, beta, rho0);
  }

  // And cool down more to make the assignment harder
  double betastop = 1./params.Tstop;
  while ((*beta) < betastop){
    if (0==threadIdx.x && 0 == blockIdx.x) (*beta) = std::min((*beta)/params.coolingFactor, betastop);
    __syncthreads();
    thermalize(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_lowT, rho0);
    __syncthreads();
  }
  // A final assignment, before sending it to the fitter
  __syncthreads();
  set_vtx_range(ntracks, tracks, vertices, params, osumtkwt, beta);
}

#ifdef __CUDACC__
std::vector<TransientVertex> vertices(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, clusterParameters params, const std::vector<reco::TransientTrack>& t_tks, double* beta) {
  const unsigned int nv = vertices->nTrueVertex;
  for (unsigned int ivertexO = 0;ivertexO < nv; ivertexO++) {
    unsigned int ivertex = vertices->order(ivertexO);
    //std::cout << "V:" << ivertexO << " " << ivertex << " z: " << vertices->z(ivertex) << std::endl;
  }
  for (unsigned int ivertexO = 0;ivertexO < nv; ivertexO++) {
    unsigned int ivertex = vertices->order(ivertexO);
    if ( vertices->rho(ivertex) > 10000 || abs(vertices->z(ivertex)) > 10000) {
      vertices->rho(ivertex) = 0;
      vertices->z(ivertex) = 0;
    }
  }
  std::vector<TransientVertex> clusters;
  const double mintrkweight_ = 0.5;
  double rho0 = vertices->nTrueVertex > 1 ? 1./vertices->nTrueVertex : 1.;
  const auto z_sum_init = rho0*exp(-(*beta)*params.dzCutOff*params.dzCutOff);
  std::vector<std::vector<unsigned int> > vtx_track_indices(vertices->nTrueVertex);

  for (unsigned int i = 0; i < ntracks; i++) {
    if (not(tracks->isGood(i))) continue;
    const auto kmin = tracks->kmin(i);
    const auto kmax = tracks->kmax(i);
    for (auto k = kmin; k < kmax; k++) {
      unsigned int ivertex = vertices->order(k);
      vertices->exp(ivertex) = exp(-(*beta) * std::pow( tracks->z(i) - vertices->z(ivertex), 2) * tracks->dz2(i)) ;
    }

    //local_exp_list_range(y.exp_arg, y.exp, kmin, kmax);


    tracks->sum_Z(i) = z_sum_init;
    for (auto k = kmin; k < kmax; k++) {
      unsigned int ivertex = vertices->order(k);
      tracks->sum_Z(i) += vertices->rho(ivertex) * vertices->exp(ivertex);
    }
    const double invZ = tracks->sum_Z(i) > 1e-100 ? 1. / tracks->sum_Z(i) : 0.0;

    for (auto k = kmin; k < kmax; k++) {
      unsigned int ivertex = vertices->order(k);
      double p = vertices->rho(ivertex) * vertices->exp(ivertex) *invZ; 
      if (p > mintrkweight_) {
        // assign  track i -> vertex k (hard, mintrkweight_ should be >= 0.5 here
        vtx_track_indices[k].push_back(i);
        break;
      }
    }

  }  // track loop
  /*
  for (unsigned int itrack=0; itrack<ntracks; itrack++){
    if (not(tracks->isGood(itrack))) continue;
    for (auto k = tracks->kmin(itrack); k<tracks->kmax(itrack); k++){
      vtx_track_indices[k].push_back(itrack);
    }
  }
  */
  GlobalError dummyError(0.01, 0, 0.01, 0., 0., 0.01);
  for (unsigned int k = 0; k < vertices->nTrueVertex; k++) {
    if (!vtx_track_indices[k].empty()) {
      unsigned int ivertex = vertices->order(k);
      GlobalPoint pos(0, 0, vertices->z(ivertex));
      std::vector<reco::TransientTrack> vertexTracks;
      for (auto i : vtx_track_indices[k]) {
        //std::cout << "Adding track " << i << " to vertex " << k << std::endl;
        vertexTracks.push_back(t_tks[i]);
      }
      TransientVertex v(pos, dummyError, vertexTracks, 0);
      clusters.push_back(v);
    }
  }
  return clusters;
}

std::vector<std::vector<reco::TransientTrack>> clusterize(std::vector<TransientVertex>& pv, clusterParameters params) {
  std::vector<std::vector<reco::TransientTrack> > clusters;
  // vector<TransientVertex>&& pv = vertices(tracks);


  if (pv.empty()) {
    return clusters;
  }

  // fill into clusters and merge
  std::vector<reco::TransientTrack> aCluster = pv.begin()->originalTracks();

  for (auto k = pv.begin() + 1; k != pv.end(); k++) {
    if (std::abs(k->position().z() - (k - 1)->position().z()) > (2 * params.vertexSize)) {
      // close a cluster
      if (aCluster.size() > 1) {
        clusters.push_back(aCluster);
      }
#ifdef DEBUG
      else {
        std::cout << " one track cluster at " << k->position().z() << "  suppressed" << std::endl;
      }
#endif
      aCluster.clear();
    }
    for (unsigned int i = 0; i < k->originalTracks().size(); i++) {
      aCluster.push_back(k->originalTracks()[i]);
    }
  }
  clusters.emplace_back(std::move(aCluster));

  return clusters;
}
  // Only on GPUs, of course...
void initializeWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Initialize is vectorized across tracks
  unsigned int blockSize = 512;
  unsigned int gridSize  = 1; // 512 vertices always. TODO::Need to get a way of setting this up based on TrackForPV::VertexForPVSoA::stride(), but the object here is on device, not on host
  initializeKernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params);  
  cudaCheck(cudaGetLastError());
}

void getBeta0Wrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Beta0 is vectorized across tracks
  unsigned int blockSize = 512; //Run 512 threads in a single block
  unsigned int gridSize  = 1;
  getBeta0Kernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params, beta);
  cudaCheck(cudaGetLastError());
}

void thermalizeWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Thermalize is vectorized -mostly- across tracks
  unsigned int blockSize = 512;
  unsigned int gridSize  = 1;
  thermalizeKernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params, osumtkwt, beta, params.delta_highT, 0.0); // At the first iteration, rho0=0, no purging of any kind
  cudaCheck(cudaGetLastError());
}


void coolingWhileSplittingWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Vectorized across tracks, as it calls thermalize and update copiusly
  unsigned int blockSize = 512;
  unsigned int gridSize  = 1;
  coolingWhileSplittingKernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params, osumtkwt, beta); // At the first iteration, rho0=0, no purging of any kind 
  cudaCheck(cudaGetLastError());
}

void remergeTracksWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Vectorized across tracks, as it calls update
  unsigned int blockSize = 512;
  unsigned int gridSize  = 1;
  remergeTracksKernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params, osumtkwt, beta); 
  cudaCheck(cudaGetLastError());
}

void resplitTracksWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Vectorized across tracks, as it calls update
  unsigned int blockSize = 512;
  unsigned int gridSize  = 1;
  resplitTracksKernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params, osumtkwt, beta);
  cudaCheck(cudaGetLastError());
}

void outlierRejectionWrapper(unsigned int ntracks, TrackForPV::TrackForPVSoA* tracks, TrackForPV::VertexForPVSoA* vertices, double* beta, double* osumtkwt, clusterParameters params, cudaStream_t stream){
  // Vectorized across tracks, as it calls update
  unsigned int blockSize = 512;
  unsigned int gridSize  = 1;
  outlierRejectionKernel<<<gridSize, blockSize,0,stream>>>(ntracks, tracks, vertices, params, osumtkwt, beta);
  cudaCheck(cudaGetLastError());
}

#endif
} // namespace clusterizerCUDA
