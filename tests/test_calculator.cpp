#include "../src/physics/calculator.hpp"
#include <exception>
#include <iostream>

int main() {
    try {
        Random rng_;

        Config config;
        config.set<int>("lattice.system_size", 16);
        config.set<int>("lattice.shell_count", 2);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_conditions", "periodic");
        config.set<double>("physical.concentration", 1.0);
        config.set<std::vector<double>>("physical.exchange_constants",
                                        {-1.0, -0.5});

        lattice::Atoms atoms(config, rng_);
        physics::Calculator calculator(atoms, config);

        std::cout << "Ferromagnetic\n";
        atoms.initialize_ferromagnetic();
        double E_fm = calculator.calculate_total_energy();
        double E_fm_min = calculator.calculate_minimal_energy();
        double M_fm = calculator.calculate_total_magnetization();
        double O_fm = calculator.calculate_total_order();
        std::cout << "Energy: " << E_fm << "\n";
        std::cout << "Minimal Energy: " << E_fm_min << "\n";
        std::cout << "Magnetization: " << M_fm << "\n";
        std::cout << "Order: " << O_fm << "\n\n";

        std::cout << "Antiferromagnetic\n";
        atoms.initialize_antiferromagnetic();
        double E_neel = calculator.calculate_total_energy();
        double E_neel_min = calculator.calculate_minimal_energy();
        double M_neel = calculator.calculate_total_magnetization();
        double O_neel = calculator.calculate_total_order();
        std::cout << "Energy: " << E_neel << "\n";
        std::cout << "Minimal Energy: " << E_neel_min << "\n";
        std::cout << "Magnetization: " << M_neel << "\n";
        std::cout << "Order: " << O_neel << "\n\n";

        std::cout << "Striped\n";
        atoms.initialize_stripped();
        double E_stripe = calculator.calculate_total_energy();
        double E_stripe_min = calculator.calculate_minimal_energy();
        double M_stripe = calculator.calculate_total_magnetization();
        double O_stripe = calculator.calculate_total_order();
        std::cout << "Energy: " << E_stripe << "\n";
        std::cout << "Minimal Energy: " << E_stripe_min << "\n";
        std::cout << "Magnetization: " << M_stripe << "\n";
        std::cout << "Order: " << O_stripe << "\n\n";

        std::cout << "Random\n";
        atoms.initialize_random();
        double E_rand = calculator.calculate_total_energy();
        double E_rand_min = calculator.calculate_minimal_energy();
        double M_rand = calculator.calculate_total_magnetization();
        double O_rand = calculator.calculate_total_order();
        std::cout << "Energy: " << E_rand << "\n";
        std::cout << "Minimal Energy: " << E_rand_min << "\n";
        std::cout << "Magnetization: " << M_rand << "\n";
        std::cout << "Order: " << O_rand << "\n";
        atoms.flip_spin(0);
        double E_rand_flip = calculator.calculate_flip_energy(0);
        std::cout << "Shift energy 1 spin: " << E_rand_flip << "\n";
        std::cout << "Energy after flip 1 spin: " << E_rand + E_rand_flip
                  << "\n\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}