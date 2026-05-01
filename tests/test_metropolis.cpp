#include "../src/algorithm/metropolis.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>

int main() {
    try {
        Config config;
        config.set<int>("lattice.system_size", 32);
        config.set<int>("lattice.shell_count", 2);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_conditions", "periodic");
        config.set<double>("physical.concentration", 0.6);
        config.set<std::vector<double>>("physical.exchange_constants",
                                        {-1.0, -1.0});
        config.set<int>("simulation.number_sweeps", 1'000);
        config.set<int>("simulation.number_term_sweeps", 10'000);

        algorithm::Metropolis metropolis(config);
        metropolis.set_output_dir("./metropolis_snaps");

        metropolis.run();

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}