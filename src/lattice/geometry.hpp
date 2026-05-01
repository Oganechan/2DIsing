#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include "../../external/config/config.hpp"
#include <array>
#include <vector>

namespace lattice {

enum class CrystalType { RECTANGULAR };
enum class BoundaryType { PERIODIC, HARD };

class Geometry {
  public:
    explicit Geometry(const Config &config);

    int get_system_size() const { return system_size_; }
    int get_shell_count() const { return shell_count_; }
    int get_basis_count() const { return basis_count_; }
    int get_atom_count() const { return atom_count_; }

    const auto &get_atom_indices() const { return atom_indices_; }
    const auto &get_atom_positions() const { return atom_positions_; }
    const auto &get_neighbor_table() const { return neighbor_table_; }

  private:
    const int system_size_, shell_count_;
    int basis_count_, atom_count_;

    const CrystalType crystal_type_;
    const BoundaryType boundary_type_;

    std::vector<std::array<double, 2>> basis_vectors_;
    std::vector<std::array<double, 2>> lattice_vectors_;

    std::vector<std::array<int, 3>> atom_indices_;
    std::vector<std::array<double, 2>> atom_positions_;
    std::vector<std::vector<std::vector<int>>> neighbor_table_;

    void initialize_lattice(const CrystalType type);
    void compute_indices();
    void compute_positions();
    void compute_neighbors();
};

} // namespace lattice

#endif // GEOMETRY_HPP