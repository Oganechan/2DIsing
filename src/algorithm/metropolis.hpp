#ifndef METROPOLIS_HPP
#define METROPOLIS_HPP

#include "../physics/calculator.hpp"
#include <filesystem>
#include <vector>

namespace algorithm {

class Metropolis {
  public:
    explicit Metropolis(const Config &config);

    void run();

    void set_output_dir(const std::filesystem::path &dir) { output_dir_ = dir; }

  private:
    Random rng_;

    std::filesystem::path output_dir_;

    lattice::Atoms atoms_;
    const physics::Calculator calc_;

    const int num_sweeps_;
    const int num_term_sweeps_;

    std::vector<double> temperatures_;
    double current_temperature_;
    double current_energy_;
    int snapshot_counter_ = 0;

    void sweep();
    void step();

    void make_snap();
};

} // namespace algorithm

#endif // METROPOLIS_HPP