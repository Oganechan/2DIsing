#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "config.hpp"

namespace lattice
{

  enum class CrystalType
  {
    RECTANGULAR,
    RECTANGULAR_ANISOTROPIC,
    TRIANGULAR
  };

  class Topology
  {
  public:
    explicit Topology(const Config &config);

    int get_system_size() const { return system_size_; }
    int get_shell_count() const { return shell_count_; }
    int get_basis_count() const { return basis_count_; }
    int get_atom_count() const { return atom_count_; }

    CrystalType get_crystal_type() const { return crystal_type_; }

    auto &get_lattice_vectors() const { return lattice_vectors_; }
    auto &get_directions() const { return directions_; }

    auto &get_atom_indices() const { return atom_indices_; }
    auto &get_atom_positions() const { return atom_positions_; }

  private:
    const int system_size_, shell_count_;
    int basis_count_, atom_count_;

    const CrystalType crystal_type_;

    std::vector<std::array<double, 2>> basis_vectors_;
    std::vector<std::array<double, 2>> lattice_vectors_;
    std::vector<std::pair<int, int>> directions_;

    std::vector<std::array<int, 3>> atom_indices_;
    std::vector<std::array<double, 2>> atom_positions_;

    void initialize_lattice();
    void compute_indices();
    void compute_positions();
  };

} // namespace lattice

#endif // TOPOLOGY_HPP