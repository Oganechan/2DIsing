#include "../external/config/config.hpp"
#include "simulation.hpp"
#include <cstdlib>
#include <exception>
#include <iostream>
#include <omp.h>

int main() {
    try {
        Config config;
        config.load("config.json");

        omp_set_num_threads(config.get<int>("simulation.number_threads"));
        int max_threads = omp_get_max_threads();
        std::cout << "Max threads: " << max_threads << "\n";

        auto start_time = std::chrono::high_resolution_clock::now();

        Simulation sim(config);
        std::cout << "Simulation started...\n\n";
        sim.run();
        sim.save_results();

        std::cout << "\nSimulation completed!\n";

        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(
            end_time - start_time);

        std::cout << "Total time: " << total_duration.count() << " seconds\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}