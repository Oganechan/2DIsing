#ifndef WANGLANDAU_HPP
#define WANGLANDAU_HPP

#include "calculator.hpp"
#include <quadmath.h>
#include <unordered_map>

namespace algorithm {

class WangLandau {
  public:
    explicit WangLandau(const Config &config);

    void run();

    __float128 compute_partition_function(double T) const;
    __float128 compute_average_energy(double T) const;
    __float128 compute_free_energy(double T) const;
    __float128 compute_entropy(double T) const;
    __float128 compute_heat_capacity(double T) const;
    __float128 compute_energy_cumulant(double T) const;

    double find_Tc() const;

    auto get_ln_g() const {
        std::unordered_map<int, double> ln_g_norm = ln_g_;
        for (auto &[key, ln_g] : ln_g_norm)
            ln_g -= S0_;

        return ln_g_norm;
    }

  private:
    Random rng_;

    lattice::Atoms atoms_;
    const physics::Calculator calc_;

    const double final_ln_f_;
    const double shake_probability_;

    std::unordered_map<int, double> ln_g_;
    std::unordered_map<int, int> histogram_;

    double ln_f_ = 1.0;
    double current_energy_;
    __float128 S0_;

    void sweep();
    void step();

    bool try_flip(int old_key, int new_key, double new_energy);

    int sweeps_since_check_ = 0;
    void check_flatness();
};

} // namespace algorithm

#endif // WANGLANDAU_HPP
