#ifndef ATOMS_HPP
#define ATOMS_HPP

#include "geometry.hpp"
#include "spin.hpp"
#include <vector>

namespace lattice {

class Atoms {
  public:
    explicit Atoms(const Config &config, Random &rng);

    void set_random_defects();
    void initialize_random();
    void initialize_ferromagnetic();
    void initialize_antiferromagnetic();
    void initialize_stripped();

    void flip_spin(int id);

    double get_concentration() const { return concentration_; }
    int get_magnetic_count() const { return magnetic_count_; }
    bool get_magnetic_state(int id) const { return magnetic_mask_[id]; }

    auto &get_geometry() const { return geometry_; }
    auto &get_spins() const { return spins_; }
    auto &get_magnetic_atoms() const { return magnetic_atoms_; }

  private:
    Random &rng_;

    const Geometry geometry_;
    const double concentration_;

    std::vector<IsingSpin> spins_;
    std::vector<bool> magnetic_mask_;
    std::vector<int> magnetic_atoms_;
    int magnetic_count_ = 0;

    void update_magnetic_list();
};

} // namespace lattice

#endif // ATOMS_HPP