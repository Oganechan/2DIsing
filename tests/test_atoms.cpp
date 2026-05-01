#include "../src/lattice/atoms.hpp"
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>

int main() {
    try {
        Random rng;

        Config config;
        config.set<int>("lattice.system_size", 100);
        config.set<int>("lattice.shell_count", 1);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_conditions", "periodic");
        config.set<double>("physical.concentration", 1.0);

        lattice::Atoms atoms0(config, rng);

        int total_atoms = atoms0.get_geometry().get_atom_count();
        int magnetic_count = atoms0.get_magnetic_count();

        std::cout << "Total atoms: " << total_atoms << "\n";
        std::cout << "Magnetic atoms: " << magnetic_count << "\n";

        atoms0.initialize_random();

        std::cout << "First 10 spins: ";
        for (size_t i = 0; i < 10; ++i) {
            std::cout << (int)atoms0.get_spins()[i].value() << " ";
        }
        std::cout << "\n";

        config.set<double>("physical.concentration", 0.9);
        lattice::Atoms atoms1(config, rng);
        magnetic_count = atoms0.get_magnetic_count();
        std::cout << "Magnetic atoms (c=0.9): " << magnetic_count << "\n";

        config.set<double>("physical.concentration", 0.5);
        lattice::Atoms atoms2(config, rng);
        magnetic_count = atoms2.get_magnetic_count();
        std::cout << "Magnetic atoms (c=0.5): " << magnetic_count << "\n";

        config.set<double>("physical.concentration", 0.1);
        lattice::Atoms atoms3(config, rng);
        magnetic_count = atoms3.get_magnetic_count();
        std::cout << "Magnetic atoms (c=0.1): " << magnetic_count << "\n";

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}