#include "calculator.hpp"
#include <cmath>
#include <cstdlib>
#include <omp.h>
#include <stdexcept>
#include <vector>

physics::Calculator::Calculator(const lattice::Atoms &atoms,
                                const Config &config)
    : atoms_(atoms), geometry_(atoms.get_geometry()),
      exchange_constants_(
          config.get<std::vector<double>>("physical.exchange_constants")) {
    if (exchange_constants_.size() < geometry_.get_shell_count())
        throw std::runtime_error("Not enough exchange constants");

    E_min = precalc_minimal_energy();
    build_weighted_neighbors();
}

double physics::Calculator::calculate_total_energy() const {
    double total_energy = 0.0;
    const auto &spins = atoms_.get_spins();

    // #pragma omp parallel for reduction(+ : total_energy)
    for (int atom_id : atoms_.get_magnetic_atoms()) {
        const auto &spin1 = spins[atom_id];

        for (int shell = 0; shell < geometry_.get_shell_count(); ++shell) {
            double J = exchange_constants_[shell];

            for (int neighbor_id :
                 geometry_.get_neighbor_table()[atom_id][shell])
                if (atoms_.get_magnetic_state(neighbor_id) &&
                    neighbor_id > atom_id) {
                    const auto &spin2 = spins[neighbor_id];
                    total_energy += -J * static_cast<double>(spin1.dot(spin2));
                }
        }
    }

    return total_energy;
}

double physics::Calculator::calculate_flip_energy(int atom_id) const {
    double delta = 0.0;
    const auto &spins = atoms_.get_spins();
    const auto &spin1 = spins[atom_id];

    for (const auto &[neighbor_id, J] : weighted_neighbors_[atom_id]) {
        const auto &spin2 = spins[neighbor_id];
        delta += 2.0 * J * static_cast<double>(spin1.dot(spin2));
    }

    return delta;
}

double physics::Calculator::calculate_total_magnetization() const {
    double total_magnetization = 0.0;
    const auto &spins = atoms_.get_spins();

    for (int atom_id : atoms_.get_magnetic_atoms())
        total_magnetization += static_cast<double>(spins[atom_id].value());

    return total_magnetization;
}

double physics::Calculator::calculate_total_order() const {
    const auto &spins = atoms_.get_spins();
    const auto &indices = geometry_.get_atom_indices();

    std::array<double, 4> sub_mag = {0.0, 0.0, 0.0, 0.0};
    std::array<int, 4> sub_count = {0, 0, 0, 0};

    for (int atom_id : atoms_.get_magnetic_atoms()) {
        auto [x, y, incel] = indices[atom_id];
        int sub = (x % 2) * 2 + (y % 2);

        sub_mag[sub] += static_cast<double>(spins[atom_id].value());
        sub_count[sub]++;
    }

    for (int i = 0; i < 4; ++i)
        if (sub_count[i] > 0)
            sub_mag[i] /= sub_count[i];

    double m_avg = (sub_mag[0] + sub_mag[1] + sub_mag[2] + sub_mag[3]) / 4.0;
    double order = 0.0;
    for (int i = 0; i < 4; ++i) {
        double diff = sub_mag[i] - m_avg;
        order += diff * diff;
    }

    return std::sqrt(order / 4.0);
}

double physics::Calculator::precalc_minimal_energy() const {
    double E_min = 0.0;
    const auto &neighbor_table = geometry_.get_neighbor_table();

    for (int atom_id : atoms_.get_magnetic_atoms())
        for (int shell = 0; shell < geometry_.get_shell_count(); ++shell) {
            const double J = exchange_constants_[shell];
            for (int neighbor_id : neighbor_table[atom_id][shell])
                if (atoms_.get_magnetic_state(neighbor_id) &&
                    neighbor_id > atom_id)
                    E_min -= std::abs(J);
        }

    return E_min;
}

void physics::Calculator::build_weighted_neighbors() {
    const auto &neighbor_table = geometry_.get_neighbor_table();
    weighted_neighbors_.resize(geometry_.get_atom_count());

    for (int atom_id : atoms_.get_magnetic_atoms()) {
        auto &neighbors = weighted_neighbors_[atom_id];
        for (int shell = 0; shell < geometry_.get_shell_count(); ++shell) {
            double J = exchange_constants_[shell];
            for (int neighbor_id : neighbor_table[atom_id][shell]) {
                double effective_J =
                    atoms_.get_magnetic_state(neighbor_id) ? J : 0.0;
                neighbors.emplace_back(neighbor_id, effective_J);
            }
        }
    }
}