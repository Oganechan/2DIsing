#ifndef SIMULATION_HPP
#define SIMULATION_HPP

#include "config.hpp"
#include <filesystem>
#include <map>
#include <vector>

namespace fs = std::filesystem;

class Simulation {
  public:
    Simulation(Config &base_cfg);

    void run();
    void save_results() const;

  private:
    Config base_cfg_;
    fs::path output_dir_;

    const double concentration_;
    const double final_ln_f_;
    const int n_realizations_;
    const std::vector<int> lattice_sizes_;
    const std::vector<double> r_values_;

    std::map<double, std::map<int, double>> all_Tc_;

    std::string generate_simulation_name() const;
    void config_backup();
};

#endif // SIMULATION_HPP