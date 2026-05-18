#include "atoms.hpp"
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

lattice::Atoms::Atoms(const Config &config, Random &rng)
    : rng_(rng), geometry_(config),
      concentration_(config.get<double>("physical.concentration")) {

    int atom_count = geometry_.get_atom_count();
    spins_.resize(atom_count);
    magnetic_mask_.resize(atom_count, true);
    magnetic_atoms_.reserve(atom_count);

    initialize_random();
    set_random_defects();
}

void lattice::Atoms::update_magnetic_list() {
    magnetic_atoms_.clear();
    magnetic_count_ = 0;
    for (size_t i = 0; i < magnetic_mask_.size(); ++i) {
        if (magnetic_mask_[i]) {
            magnetic_atoms_.push_back(i);
            ++magnetic_count_;
        }
    }
}

void lattice::Atoms::set_random_defects() {
    const int atom_count = geometry_.get_atom_count();
    std::fill(magnetic_mask_.begin(), magnetic_mask_.end(), true);

    const int defect_count =
        static_cast<int>((1.0 - concentration_) * atom_count);

    std::vector<int> atom_indices(atom_count);
    std::iota(atom_indices.begin(), atom_indices.end(), 0);

    for (size_t i = atom_count - 1; i > 0; --i) {
        size_t j = rng_.uniform_int<size_t>(0, i);
        std::swap(atom_indices[i], atom_indices[j]);
    }

    for (size_t i = 0; i < defect_count; ++i)
        magnetic_mask_[atom_indices[i]] = false;

    update_magnetic_list();
}

void lattice::Atoms::initialize_random() {
    for (size_t i = 0; i < spins_.size(); ++i)
        spins_[i].randomize(rng_);
}

void lattice::Atoms::initialize_ferromagnetic() {
    for (size_t i = 0; i < spins_.size(); ++i)
        spins_[i].set(1);
}

void lattice::Atoms::initialize_antiferromagnetic() {
    const auto &indices = geometry_.get_atom_indices();
    for (size_t i = 0; i < spins_.size(); ++i) {
        auto [x, y, incel] = indices[i];
        int spin = ((x + y) % 2 == 0) ? 1 : -1;
        spins_[i].set(spin);
    }
}

void lattice::Atoms::initialize_stripped() {
    const auto &indices = geometry_.get_atom_indices();
    for (size_t i = 0; i < spins_.size(); ++i) {
        auto [x, y, incel] = indices[i];
        int spin = (x % 2 == 0) ? 1 : -1;
        spins_[i].set(spin);
    }
}

void lattice::Atoms::flip_spin(int id) { spins_[id].flip(); }