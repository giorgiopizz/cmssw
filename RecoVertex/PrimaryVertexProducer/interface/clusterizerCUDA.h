#ifndef clusterizerCUDA_h
#define clusterizerCUDA_h
#include "CUDADataFormats/Track/interface/TrackForPVHeterogeneous.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include <cstddef>
#include <cstdint>
#include "CUDADataFormats/Vertex/interface/ZVertexHeterogeneous.h"

struct track_t {
    std::vector<double> zpca_vec;                  // z-coordinate at point of closest approach to the beamline
    std::vector<double> dz2_vec;                   // square of the error of z(pca)
    std::vector<double> sum_Z_vec;                 // track contribution to the partition function, Z
    std::vector<double> tkwt_vec;                  // track weight, close to 1.0 for most tracks
    std::vector<unsigned int> kmin;                // index of the first cluster within zrange
    std::vector<unsigned int> kmax;                // 1 + index of the last cluster within zrange
    std::vector<const reco::TransientTrack *> tt;  // a pointer to the Transient Track

    double osumtkwt;  // 1. / (sum of all track weights)

    void addItemSorted(double new_zpca, double new_dz2, const reco::TransientTrack *new_tt, double new_tkwt) {
      // sort tracks with decreasing resolution (note that dz2 = 1/sigma^2)
      unsigned int i = 0;
      for (i = 0; i < zpca_vec.size(); i++) {
        if (new_dz2 > dz2_vec[i])
          break;
      }
      insertItem(i, new_zpca, new_dz2, new_tt, new_tkwt);
    }

    void insertItem(
        unsigned int i, double new_zpca, double new_dz2, const reco::TransientTrack *new_tt, double new_tkwt) {
      zpca_vec.insert(zpca_vec.begin() + i, new_zpca);
      dz2_vec.insert(dz2_vec.begin() + i, new_dz2);
      tt.insert(tt.begin() + i, new_tt);
      tkwt_vec.insert(tkwt_vec.begin() + i, new_tkwt);
      sum_Z_vec.insert(sum_Z_vec.begin() + i, 1.0);
      kmin.insert(kmin.begin() + i, 0);
      kmax.insert(kmax.begin() + i, 0);
    }

    unsigned int getSize() const { return zpca_vec.size(); }

    // has to be called everytime the items are modified
    void extractRaw() {
      zpca = &zpca_vec.front();
      dz2 = &dz2_vec.front();
      tkwt = &tkwt_vec.front();
      sum_Z = &sum_Z_vec.front();
      kmin_in = &kmin.front();
      kmax_in = &kmax.front();
    }

    // pointers to the first element of vectors, needed for vectorized code
    double *__restrict__ zpca;
    double *__restrict__ dz2;
    double *__restrict__ tkwt;
    double *__restrict__ sum_Z;
    unsigned int *__restrict__ kmin_in;
    unsigned int *__restrict__ kmax_in;
};

struct vertex_t {
    std::vector<double> zvtx_vec;  // z coordinate
    std::vector<double> rho_vec;   // vertex "mass" for mass-constrained clustering
    // --- temporary numbers, used during update
    std::vector<double> exp_arg_vec;
    std::vector<double> exp_vec;
    std::vector<double> sw_vec;
    std::vector<double> swz_vec;
    std::vector<double> se_vec;
    std::vector<double> swE_vec;

    unsigned int getSize() const { return zvtx_vec.size(); }

    void addItem(double new_zvtx, double new_rho) {
      zvtx_vec.push_back(new_zvtx);
      rho_vec.push_back(new_rho);
      exp_arg_vec.push_back(0.0);
      exp_vec.push_back(0.0);
      sw_vec.push_back(0.0);
      swz_vec.push_back(0.0);
      se_vec.push_back(0.0);
      swE_vec.push_back(0.0);

      extractRaw();
    }

    void insertItem(unsigned int k, double new_zvtx, double new_rho, track_t &tks) {
      zvtx_vec.insert(zvtx_vec.begin() + k, new_zvtx);
      rho_vec.insert(rho_vec.begin() + k, new_rho);

      exp_arg_vec.insert(exp_arg_vec.begin() + k, 0.0);
      exp_vec.insert(exp_vec.begin() + k, 0.0);
      sw_vec.insert(sw_vec.begin() + k, 0.0);
      swz_vec.insert(swz_vec.begin() + k, 0.0);
      se_vec.insert(se_vec.begin() + k, 0.0);
      swE_vec.insert(swE_vec.begin() + k, 0.0);

      // adjust vertex lists of tracks
      for (unsigned int i = 0; i < tks.getSize(); i++) {
        if (tks.kmin[i] > k) {
          tks.kmin[i]++;
        }
        if ((tks.kmax[i] >= k) || (tks.kmax[i] == tks.kmin[i])) {
          tks.kmax[i]++;
        }
      }

      extractRaw();
    }

    void removeItem(unsigned int k, track_t &tks) {
      zvtx_vec.erase(zvtx_vec.begin() + k);
      rho_vec.erase(rho_vec.begin() + k);

      exp_arg_vec.erase(exp_arg_vec.begin() + k);
      exp_vec.erase(exp_vec.begin() + k);
      sw_vec.erase(sw_vec.begin() + k);
      swz_vec.erase(swz_vec.begin() + k);
      se_vec.erase(se_vec.begin() + k);
      swE_vec.erase(swE_vec.begin() + k);

      // adjust vertex lists of tracks
      for (unsigned int i = 0; i < tks.getSize(); i++) {
        if (tks.kmax[i] > k) {
          tks.kmax[i]--;
        }
        if ((tks.kmin[i] > k) || (((tks.kmax[i] < (tks.kmin[i] + 1)) && (tks.kmin[i] > 0)))) {
          tks.kmin[i]--;
        }
      }

      extractRaw();
    }

    // pointers to the first element of vectors, needed for vectorized code
    double *__restrict__ zvtx;
    double *__restrict__ rho;
    double *__restrict__ exp_arg;
    double *__restrict__ exp;
    double *__restrict__ sw;
    double *__restrict__ swz;
    double *__restrict__ se;
    double *__restrict__ swE;

    // has to be called everytime the items are modified
    void extractRaw() {
      zvtx = &zvtx_vec.front();
      rho = &rho_vec.front();
      exp = &exp_vec.front();
      sw = &sw_vec.front();
      swz = &swz_vec.front();
      se = &se_vec.front();
      swE = &swE_vec.front();
      exp_arg = &exp_arg_vec.front();
    }
};

namespace clusterizerCUDA {
  void kernel_update_wrapper(unsigned int nt, unsigned int nv, double beta, bool updateTc, double rho0, double osumtkwt, double Z_init, const unsigned int * tkmin, const unsigned int * tkmax, double * tz, double * tdz2, double * twgt, double * vrho, double * vz, double * tsumz, double * vsw, double * vswe, double * delta, cudaStream_t stream);
}

#endif
