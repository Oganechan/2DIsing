#include "../src/algorithm/wanglandau.hpp"
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <quadmath.h>

std::string float128_to_string(__float128 value) {
    char buf[128];
    quadmath_snprintf(buf, sizeof(buf), "%.4Qg", value);
    return std::string(buf);
}

int main() {
    try {
        double J1 = -1.0, J2 = 0.0;
        int L = 12;

        Config config;
        config.set<int>("lattice.system_size", L);
        config.set<int>("lattice.shell_count", 2);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_conditions", "periodic");
        config.set<double>("physical.concentration", 1.0);
        config.set<std::vector<double>>("physical.exchange_constants",
                                        {J1, J2});
        config.set<double>("simulation.final_ln_f", 1.0e-6);
        config.set<double>("simulation.shake_probability", 0.001);

        algorithm::WangLandau wl(config);
        std::cout << "Start WL Simulation...\n";

        std::chrono::time_point start =
            std::chrono::high_resolution_clock::now();

        wl.run();

        std::chrono::time_point end = std::chrono::high_resolution_clock::now();
        long duration =
            std::chrono::duration_cast<std::chrono::seconds>(end - start)
                .count();

        double T_c = wl.find_Tc();
        __float128 C = wl.compute_heat_capacity(T_c);

        std::cout << "\n=== Simulation Results ===\n";
        std::cout << "Energy states visited: " << wl.get_ln_g().size() << "\n";
        std::cout << "Critical Temperature: " << T_c << "\n";
        std::cout << "Heat capacity at Tc: " << float128_to_string(C) << "\n";
        std::cout << "Total time: " << duration << " seconds\n";

        std::ofstream result_file("test_wl_results.csv");
        result_file << "T,Z,E,C,U,S\n";

        result_file << std::scientific << std::setprecision(8);

        const double min_T = 0.01, max_T = 100.0;
        const int num_steps = 500;
        for (int i = 0; i < num_steps; ++i) {
            double exponent =
                std::log10(min_T) +
                i * (std::log10(max_T) - std::log10(min_T)) / (num_steps - 1);
            double T = std::pow(10, exponent);
            __float128 Z = wl.compute_partition_function(T);
            __float128 E = wl.compute_average_energy(T);
            __float128 C = wl.compute_heat_capacity(T);
            __float128 U = wl.compute_energy_cumulant(T);
            __float128 S = wl.compute_entropy(T);

            result_file << T << ",";
            result_file << float128_to_string(Z) << ",";
            result_file << float128_to_string(E) << ",";
            result_file << float128_to_string(C) << ",";
            result_file << float128_to_string(U) << ",";
            result_file << float128_to_string(S) << "\n";
        }
        result_file.close();

        auto ln_g_ = wl.get_ln_g();
        std::vector<int64_t> keys;
        for (const auto &[key, ln_g] : ln_g_) {
            keys.push_back(key);
        }
        std::sort(keys.begin(), keys.end());

        std::ofstream ln_g_file("test_wl_ln_g.csv");
        ln_g_file << "key,E,ln_g,g\n";
        for (int key : keys) {
            double ln_g = ln_g_.at(key);
            double E = key / 100.0;

            ln_g_file << key << "," << std::fixed << std::setprecision(2) << E
                      << "," << std::scientific << std::setprecision(10) << ln_g
                      << "\n";
        }
        ln_g_file.close();

        system("python3 ../scripts/visualize_plots.py test_wl_results.csv");
        system("python3 ../scripts/visualize_ln_g.py");

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}