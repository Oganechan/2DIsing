#include "geometry.hpp"
#include <array>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>

namespace {

std::array<int, 3> get_cell_index(int atom_id, int basis_count,
                                  int system_size) {
    int atom_in_cell_id = atom_id % basis_count;
    int temp = atom_id / basis_count;
    int cell_i = temp % system_size;
    int cell_j = temp / system_size;

    return {cell_i, cell_j, atom_in_cell_id};
}

int get_atom_id(int cell_i, int cell_j, int incel, int basis_count,
                int system_size) {
    int taco_i = (cell_i % system_size + system_size) % system_size;
    int taco_j = (cell_j % system_size + system_size) % system_size;

    return (taco_i + taco_j * system_size) * basis_count + incel;
}

std::array<double, 2>
compute_atom_position(int cell_i, int cell_j, int incel,
                      const std::vector<std::array<double, 2>> &lattice_vectors,
                      const std::vector<std::array<double, 2>> &basis_vectors,
                      int system_size) {

    double x = cell_i * lattice_vectors[0][0] + cell_j * lattice_vectors[1][0] +
               basis_vectors[incel][0];
    double y = cell_i * lattice_vectors[0][1] + cell_j * lattice_vectors[1][1] +
               basis_vectors[incel][1];

    x = std::fmod(x, system_size * lattice_vectors[0][0]);
    y = std::fmod(y, system_size * lattice_vectors[1][1]);

    return {x, y};
}

double
compute_min_distance(double delta_x, double delta_y,
                     const std::vector<std::array<double, 2>> &lattice_vectors,
                     int system_size, bool periodic) {

    double min_distance_sq = delta_x * delta_x + delta_y * delta_y;

    if (periodic) {
        for (int n = -1; n <= 1; ++n) {
            for (int m = -1; m <= 1; ++m) {
                double taco_dx =
                    (n * lattice_vectors[0][0] + m * lattice_vectors[1][0]) *
                        system_size +
                    delta_x;
                double taco_dy =
                    (n * lattice_vectors[0][1] + m * lattice_vectors[1][1]) *
                        system_size +
                    delta_y;
                double distance_sq = taco_dx * taco_dx + taco_dy * taco_dy;
                if (distance_sq < min_distance_sq) {
                    min_distance_sq = distance_sq;
                }
            }
        }
    }

    return std::sqrt(min_distance_sq);
}

lattice::CrystalType parse_crystal_type(const std::string &type) {
    const std::unordered_map<std::string, lattice::CrystalType> mapping = {
        {"rectangular", lattice::CrystalType::RECTANGULAR}};

    if (auto it = mapping.find(type); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknown crystal type: " + type);
}

lattice::BoundaryType parse_boundary_type(const std::string &type) {
    const std::unordered_map<std::string, lattice::BoundaryType> mapping = {
        {"periodic", lattice::BoundaryType::PERIODIC},
        {"hard", lattice::BoundaryType::HARD}};

    if (auto it = mapping.find(type); it != mapping.end())
        return it->second;

    throw std::invalid_argument("Unknown boundary type: " + type);
}

} // namespace

lattice::Geometry::Geometry(const Config &config)
    : system_size_(config.get<int>("lattice.system_size")),
      shell_count_(config.get<int>("lattice.shell_count")),
      crystal_type_(
          parse_crystal_type(config.get<std::string>("lattice.crystal_type"))),
      boundary_type_(parse_boundary_type(
          config.get<std::string>("lattice.boundary_conditions"))) {
    initialize_lattice(crystal_type_);
    compute_indices();
    compute_positions();
    compute_neighbors();
}

void lattice::Geometry::initialize_lattice(const lattice::CrystalType type) {
    switch (type) {
    case lattice::CrystalType::RECTANGULAR:
        basis_vectors_ = {{0.0, 0.0}};
        lattice_vectors_ = {{1.0, 0.0}, {0.0, 1.0}};
        break;
    }

    basis_count_ = basis_vectors_.size();
    atom_count_ = system_size_ * system_size_ * basis_count_;
}

void lattice::Geometry::compute_indices() {
    atom_indices_.resize(atom_count_);
    for (size_t i = 0; i < atom_count_; ++i)
        atom_indices_[i] = get_cell_index(i, basis_count_, system_size_);
}

void lattice::Geometry::compute_positions() {
    atom_positions_.resize(atom_count_);
    for (size_t i = 0; i < atom_count_; ++i) {
        auto [cell_i, cell_j, incel] = atom_indices_[i];
        atom_positions_[i] =
            compute_atom_position(cell_i, cell_j, incel, lattice_vectors_,
                                  basis_vectors_, system_size_);
    }
}

void lattice::Geometry::compute_neighbors() {
    neighbor_table_.resize(atom_count_);
    for (auto &atom_neighbors : neighbor_table_)
        atom_neighbors.resize(shell_count_);

    int32_t max_cell_search = shell_count_ * 2;

    for (int32_t atom_id = 0; atom_id < atom_count_; ++atom_id) {
        std::vector<std::pair<int32_t, double>> atom_neighbors;

        auto [cell_i, cell_j, atom_in_cell_id] = atom_indices_[atom_id];

        for (int32_t di = -max_cell_search; di <= max_cell_search; ++di)
            for (int32_t dj = -max_cell_search; dj <= max_cell_search; ++dj)
                for (int32_t new_incel = 0; new_incel < basis_count_;
                     ++new_incel) {

                    int32_t neighbor_id =
                        get_atom_id(cell_i + di, cell_j + dj, new_incel,
                                    basis_count_, system_size_);

                    if (neighbor_id < 0 || neighbor_id >= atom_count_ ||
                        atom_id == neighbor_id)
                        continue;

                    if (boundary_type_ == BoundaryType::HARD) {
                        auto [ni, nj, pos] = atom_indices_[neighbor_id];
                        if (ni < 0 || ni >= system_size_ || nj < 0 ||
                            nj >= system_size_)
                            continue;
                    }

                    double delta_x = atom_positions_[neighbor_id][0] -
                                     atom_positions_[atom_id][0];
                    double delta_y = atom_positions_[neighbor_id][1] -
                                     atom_positions_[atom_id][1];

                    double distance = compute_min_distance(
                        delta_x, delta_y, lattice_vectors_, system_size_,
                        boundary_type_ == lattice::BoundaryType::PERIODIC);

                    atom_neighbors.emplace_back(neighbor_id, distance);
                }

        std::sort(
            atom_neighbors.begin(), atom_neighbors.end(),
            [](const auto &a, const auto &b) { return a.second < b.second; });

        double current_distance = atom_neighbors[0].second;
        int current_shell = 0;

        for (const auto &[neighbor_id, distance] : atom_neighbors) {
            if (std::abs(distance - current_distance) > 1.0e-5) {
                current_distance = distance;
                ++current_shell;
                if (current_shell >= shell_count_)
                    break;
            }
            if (current_shell < shell_count_)
                neighbor_table_[atom_id][current_shell].push_back(neighbor_id);
        }

        for (int32_t shell = 0; shell < shell_count_; ++shell)
            std::sort(neighbor_table_[atom_id][shell].begin(),
                      neighbor_table_[atom_id][shell].end());
    }
}