#include "wanglandau.hpp"
#include <cmath>
#include <limits>
#include <quadmath.h>
#include <vector>

namespace {

int energy_to_key(double energy) {
    return static_cast<int>(std::round(energy * 100.0));
}

double key_to_energy(int key) { return key / 100.0; }

} // namespace

algorithm::WangLandau::WangLandau(const Config &config)
    : atoms_(config, rng_), calc_(atoms_, config),
      final_ln_f_(config.get<double>("simulation.final_ln_f")),
      shake_probability_(config.get<double>("simulation.shake_probability")) {}

void algorithm::WangLandau::run() {
    ln_f_ = 1.0;
    ln_g_.clear();
    histogram_.clear();

    current_energy_ = calc_.calculate_total_energy();

    while (ln_f_ > final_ln_f_) {
        sweep();
    }

    int min_energy_key = std::numeric_limits<int>::max();
    for (const auto &[key, ln_g] : ln_g_)
        if (key < min_energy_key)
            min_energy_key = key;

    S0_ = ln_g_[min_energy_key];
}

void algorithm::WangLandau::sweep() {
    int step_count = atoms_.get_magnetic_count();
    while (step_count-- > 0)
        step();

    if (sweeps_since_check_++ >= 1'000) {
        sweeps_since_check_ = 0;
        check_flatness();
    }
}

void algorithm::WangLandau::step() {
    const auto &magnetic_atoms = atoms_.get_magnetic_atoms();

    if (rng_.bernoulli(shake_probability_)) {
        const int magnetic_count = atoms_.get_magnetic_count();
        std::vector<int> flipped_indices;
        std::vector<int> flipped_ids;

        for (int i = 0; i < 0.25 * magnetic_count; ++i) {
            int idx = rng_.uniform_int(0, magnetic_count - 1);
            int atom_id = magnetic_atoms[idx];
            atoms_.flip_spin(atom_id);
            flipped_indices.push_back(idx);
            flipped_ids.push_back(atom_id);
        }

        double old_E = current_energy_;
        double new_E = calc_.calculate_total_energy();

        int old_key = energy_to_key(old_E);
        int new_key = energy_to_key(new_E);

        if (!try_flip(old_key, new_key, new_E))
            for (int atom_id : flipped_ids)
                atoms_.flip_spin(atom_id);

        return;
    }

    int random_atom_id =
        magnetic_atoms[rng_.uniform_int<int>(0, magnetic_atoms.size() - 1)];

    double old_E = current_energy_;
    double dE = calc_.calculate_flip_energy(random_atom_id);
    double new_E = old_E + dE;

    int old_key = energy_to_key(old_E);
    int new_key = energy_to_key(new_E);

    if (try_flip(old_key, new_key, new_E))
        atoms_.flip_spin(random_atom_id);
}

bool algorithm::WangLandau::try_flip(int old_key, int new_key,
                                     double new_energy) {
    auto [it_ln_old, _] = ln_g_.try_emplace(old_key, 0.0);
    auto [it_ln_new, __] = ln_g_.try_emplace(new_key, 0.0);
    auto [it_hist_old, ___] = histogram_.try_emplace(old_key, 0);
    auto [it_hist_new, ____] = histogram_.try_emplace(new_key, 0);

    double acceptance =
        std::min(1.0, std::exp(it_ln_old->second - it_ln_new->second));
    bool accepted = rng_.bernoulli(acceptance);

    if (accepted) {
        current_energy_ = new_energy;
        it_ln_new->second += ln_f_;
        it_hist_new->second++;
    } else {
        it_ln_old->second += ln_f_;
        it_hist_old->second++;
    }

    return accepted;
}

void algorithm::WangLandau::check_flatness() {
    if (histogram_.size() == 0)
        return;

    int total = 0;
    int min_h = std::numeric_limits<int>::max();

    for (const auto &[energy, count] : histogram_) {
        total += count;
        if (count < min_h)
            min_h = count;
    }

    double avg = static_cast<double>(total) / histogram_.size();
    bool is_flat = (min_h > 0.8 * avg);

    if (is_flat) {
        ln_f_ *= 0.5;
        histogram_.clear();
    }
}

__float128 algorithm::WangLandau::compute_partition_function(double T) const {
    __float128 Z = 0.0q;
    for (const auto &[key, ln_g] : ln_g_) {
        __float128 E = key_to_energy(key);
        Z += expq(static_cast<__float128>(ln_g) - S0_ -
                  E / static_cast<__float128>(T));
    }

    return Z;
}

__float128 algorithm::WangLandau::compute_average_energy(double T) const {
    __float128 Z = 0.0q;
    __float128 E_avg = 0.0q;

    for (const auto &[key, ln_g] : ln_g_) {
        __float128 E = key_to_energy(key);
        __float128 weight = expq(static_cast<__float128>(ln_g) - S0_ -
                                 E / static_cast<__float128>(T));
        Z += weight;
        E_avg += E * weight;
    }

    return E_avg / Z;
}

__float128 algorithm::WangLandau::compute_free_energy(double T) const {
    __float128 Z = compute_partition_function(T);
    __float128 F = -static_cast<__float128>(T) * logq(Z);

    return F;
}

__float128 algorithm::WangLandau::compute_entropy(double T) const {
    int N = atoms_.get_magnetic_count();
    __float128 E = compute_average_energy(T);
    __float128 F = compute_free_energy(T);
    __float128 S = (E - F) / static_cast<__float128>(T);

    return S / N;
}

__float128 algorithm::WangLandau::compute_heat_capacity(double T) const {
    __float128 Z = 0.0q;
    __float128 E_avg = 0.0q;
    __float128 E2_avg = 0.0q;
    int N = atoms_.get_magnetic_count();

    for (const auto &[key, ln_g] : ln_g_) {
        __float128 E = key_to_energy(key);
        __float128 weight = expq(static_cast<__float128>(ln_g) - S0_ -
                                 E / static_cast<__float128>(T));
        Z += weight;
        E_avg += E * weight;
        E2_avg += E * E * weight;
    }

    E_avg /= Z;
    E2_avg /= Z;

    return (E2_avg - E_avg * E_avg) / (N * T * T);
}

__float128 algorithm::WangLandau::compute_energy_cumulant(double T) const {
    __float128 Z = 0.0q;
    __float128 E2_avg = 0.0q;
    __float128 E4_avg = 0.0q;

    for (const auto &[key, ln_g] : ln_g_) {
        __float128 E = key_to_energy(key);
        __float128 weight = expq(static_cast<__float128>(ln_g) - S0_ -
                                 E / static_cast<__float128>(T));
        Z += weight;
        E2_avg += E * E * weight;
        E4_avg += E * E * E * E * weight;
    }
    E2_avg /= Z;
    E4_avg /= Z;

    return 1.0q - E4_avg / (3.0q * E2_avg * E2_avg);
}

double algorithm::WangLandau::find_Tc() const {
    __float128 max_C = -1.0;
    double T_c = 0.0;

    for (double T = 0.0; T < 100.0; T += 0.01) {
        __float128 C = compute_heat_capacity(T);
        if (C > max_C && !isnanq(C) && !isinfq(C)) {
            max_C = C;
            T_c = T;
        }
    }

    return T_c;
}