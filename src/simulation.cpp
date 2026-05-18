#include "simulation.hpp"
#include "algorithm/wanglandau.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <omp.h>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace {
template <typename T> T compute_average(const std::vector<T> &values) {
    T sum = std::accumulate(values.begin(), values.end(), 0.0);
    return sum / values.size();
}

double compute_stddev(const std::vector<double> &values) {
    double mean = compute_average(values);

    double sum_sq = 0.0;
    for (double v : values)
        sum_sq += (v - mean) * (v - mean);

    return std::sqrt(sum_sq / (values.size() - 1));
}

double extrapolate_Tc(const std::vector<double> &inv_L,
                      const std::vector<double> &Tc) {
    if (inv_L.size() < 2)
        return Tc[0];

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_x2 = 0.0;
    int n = inv_L.size();

    for (int i = 0; i < n; ++i) {
        sum_x += inv_L[i];
        sum_y += Tc[i];
        sum_xy += inv_L[i] * Tc[i];
        sum_x2 += inv_L[i] * inv_L[i];
    }

    double slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - sum_x * sum_x);
    double intercept = (sum_y - slope * sum_x) / n;

    return intercept;
}

std::string float128_to_string(__float128 value) {
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.4Qf", value);
    return std::string(buf);
}

} // namespace

Simulation::Simulation(Config &base_cfg)
    : base_cfg_(base_cfg),
      concentration_(base_cfg.get<double>("physical.concentration")),
      final_ln_f_(base_cfg.get<double>("simulation.final_ln_f")),
      n_realizations_(base_cfg.get<int>("simulation.number_realizations")),
      lattice_sizes_(
          base_cfg.get<std::vector<int>>("simulation.lattice_sizes")),
      r_values_(base_cfg.get<std::vector<double>>("simulation.r_values")) {}

void Simulation::run() {
    const int num_steps = 1'000;
    const double min_T = 0.01, max_T = 100.0;
    const double log10_min_T = std::log10(min_T),
                 log10_max_T = std::log10(max_T);
    const double step = (log10_max_T - log10_min_T) / (num_steps - 1);

    std::string sim_name = generate_simulation_name();
    output_dir_ = fs::path("results") / sim_name;
    fs::create_directories(output_dir_.string());

    config_backup();

    std::cout << "Results will be saved to: " << output_dir_ << "\n\n";

    for (double r : r_values_) {
        std::cout << "Progressing r = " << r << "\n";

        fs::path r_path = output_dir_ / ("r" + std::to_string(r).substr(0, 4));
        fs::create_directories(r_path);

        for (int L : lattice_sizes_) {
            std::cout << " L = " << L << std::endl;

            std::vector<double> Tc_realize(n_realizations_);
            std::vector<std::vector<__float128>> E_realize(
                num_steps, std::vector<__float128>(n_realizations_));
            std::vector<std::vector<__float128>> C_realize(
                num_steps, std::vector<__float128>(n_realizations_));
            std::vector<std::vector<__float128>> S_realize(
                num_steps, std::vector<__float128>(n_realizations_));

#pragma omp parallel for
            for (int rep = 0; rep < n_realizations_; ++rep) {
                Config config = base_cfg_;
                config.set<int>("lattice.system_size", L);
                config.set<int>("lattice.shell_count", 2);
                config.set<std::vector<double>>("physical.exchange_constants",
                                                {-1.0, -r});

                algorithm::WangLandau wl(config);
                wl.run();

                for (int i = 0; i < num_steps; ++i) {
                    double exponent = log10_min_T + i * step;
                    double T = std::pow(10, exponent);
                    __float128 E = wl.compute_average_energy(T);
                    __float128 C = wl.compute_heat_capacity(T);
                    __float128 S = wl.compute_entropy(T);

                    E_realize[i][rep] = E;
                    C_realize[i][rep] = C;
                    S_realize[i][rep] = S;
                }

                double Tc = wl.find_Tc();
                Tc_realize[rep] = Tc;
            }

            std::ofstream f(r_path / ("L" + std::to_string(L) + "results.csv"));
            f << "T,E,C,S\n";
            for (int i = 0; i < num_steps; ++i) {
                double exponent = log10_min_T + i * step;
                double T = std::pow(10, exponent);

                __float128 average_E = compute_average(E_realize[i]);
                __float128 average_C = compute_average(C_realize[i]);
                __float128 average_S = compute_average(S_realize[i]);

                f << T << ",";
                f << float128_to_string(average_E) << ",";
                f << float128_to_string(average_C) << ",";
                f << float128_to_string(average_S) << "\n";
            }
            f.close();

            all_Tc_[r][L] = compute_average(Tc_realize);
        }
    }
}

void Simulation::save_results() const {
    std::ofstream raw(output_dir_ / "raw_results.csv");
    raw << "L,r,Tc\n";
    for (int L : lattice_sizes_) {
        for (double r : r_values_) {
            raw << L << "," << r << "," << all_Tc_.at(r).at(L) << "\n";
        }
    }
    raw.close();

    std::ofstream results(output_dir_ / "results.csv");
    results << "c,r,Tc\n";
    for (double r : r_values_) {
        std::vector<double> inv_L;
        std::vector<double> Tc_means;

        for (int L : lattice_sizes_) {
            auto &Tc_mean = all_Tc_.at(r).at(L);
            inv_L.push_back(1.0 / L);
            Tc_means.push_back(Tc_mean);
        }

        double Tc_extrap = extrapolate_Tc(inv_L, Tc_means);
        results << concentration_ << "," << r << "," << Tc_extrap << "\n";
    }
    results.close();
}

std::string Simulation::generate_simulation_name() const {
    std::ostringstream name;
    name << "R" << n_realizations_ << "_";
    name << "f" << final_ln_f_ << "_";
    name << "c" << std::fixed << std::setprecision(2) << concentration_ << "_";

    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm *timeinfo = std::localtime(&t);
    std::ostringstream timestamp;
    timestamp << std::put_time(timeinfo, "%Y%m%d_%H%M%S");

    name << timestamp.str();

    return name.str();
}

void Simulation::config_backup() {
    fs::path config_file_path = output_dir_ / "config_backup.json";
    std::ofstream config_file(config_file_path);
    config_file << base_cfg_.get().dump(4);
    config_file.close();
}