#include "metropolis.hpp"
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

algorithm::Metropolis::Metropolis(const Config &config)
    : atoms_(config, rng_), calc_(atoms_, config),
      num_sweeps_(config.get<int>("simulation.number_sweeps")),
      num_term_sweeps_(config.get<int>("simulation.number_term_sweeps")) {
    int num_steps = 41;
    double max_T = 2.0, min_T = 0.0;
    double step = (max_T - min_T) / (num_steps - 1);

    temperatures_.reserve(num_steps);
    for (int i = 0; i < num_steps; ++i) {
        double T = max_T - i * step;
        temperatures_.push_back(T);
    }
}

void algorithm::Metropolis::run() {
    std::filesystem::create_directories(output_dir_);
    std::cout << "Metropolis simulation started. Snapshots will be saved to: "
              << output_dir_ << "\n";

    current_energy_ = calc_.calculate_total_energy();

    for (double T : temperatures_) {
        current_temperature_ = T;
        for (int i = 0; i < num_term_sweeps_; ++i)
            sweep();

        int snapshots_per_temp = 10;
        int sweep_interval = num_sweeps_ / snapshots_per_temp;
        for (int i = 0; i < num_sweeps_; ++i) {
            sweep();

            if (i % sweep_interval == 0 &&
                snapshot_counter_ < snapshots_per_temp) {
                make_snap();
                snapshot_counter_++;
            }
        }
        snapshot_counter_ = 0;
    }

    std::cout << "Metropolis simulation completed!\n";
}

void algorithm::Metropolis::sweep() {
    int num_steps = atoms_.get_magnetic_count();
    while (num_steps-- > 0) {
        step();
    }
}

void algorithm::Metropolis::step() {
    const auto &magnetic_atoms = atoms_.get_magnetic_atoms();
    int random_idx = rng_.uniform_int<int>(0, magnetic_atoms.size() - 1);
    int atom_id = magnetic_atoms[random_idx];

    double dE = calc_.calculate_flip_energy(atom_id);

    if (dE < 0.0) {
        atoms_.flip_spin(atom_id);
        current_energy_ += dE;
    } else {
        double acceptance = std::exp(-dE / current_temperature_);
        if (rng_.bernoulli(acceptance)) {
            atoms_.flip_spin(atom_id);
            current_energy_ += dE;
        }
    }
}

void algorithm::Metropolis::make_snap() {
    std::string filename = "snapshot_T" +
                           std::to_string(current_temperature_).substr(0, 6) +
                           "_" + std::to_string(snapshot_counter_) + ".txt";
    std::filesystem::path filepath = output_dir_ / filename;

    std::ofstream f(filepath);

    const auto &spins = atoms_.get_spins();
    const auto &positions = atoms_.get_geometry().get_atom_positions();

    f << "# Metropolis snapshot\n";
    f << "# Temperature: " << current_temperature_ << "\n";
    f << "# Frustration parameter: "
      << calc_.get_exchange_constants().at(1) /
             calc_.get_exchange_constants().at(0)
      << "\n";
    f << "# Concentration: " << atoms_.get_concentration() << "\n";
    f << "# Energy: " << current_energy_ << "\n";
    f << "# Magnetization: " << calc_.calculate_total_magnetization() << "\n";
    f << "# Order parameter: " << calc_.calculate_total_order() << "\n";
    f << "# Format: x y spin magnetic(1/0)\n";

    for (size_t i = 0; i < spins.size(); ++i) {
        auto [x, y] = positions[i];
        int spin_value = static_cast<int>(spins[i].value());
        int is_magnetic = atoms_.get_magnetic_state(i) ? 1 : 0;

        f << x << " " << y << " " << spin_value << " " << is_magnetic << "\n";
    }
    f.close();
}