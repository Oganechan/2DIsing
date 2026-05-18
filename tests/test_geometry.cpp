#include "../external/config/config.hpp"
#include "../src/lattice/geometry.hpp"
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <string>

int main() {
    try {
        Config config;
        config.set<int>("lattice.system_size", 8);
        config.set<int>("lattice.shell_count", 2);
        config.set<std::string>("lattice.crystal_type", "rectangular");
        config.set<std::string>("lattice.boundary_conditions", "periodic");

        lattice::Geometry geometry(config);
        auto positions = geometry.get_atom_positions();
        auto neighbors = geometry.get_neighbor_table();

        std::ofstream pos_file("atom_positions.csv");
        pos_file << "atom_id,x,y\n";
        for (size_t i = 0; i < positions.size(); ++i) {
            auto [x, y] = positions[i];
            pos_file << i << "," << x << "," << y << "\n";
        }
        pos_file.close();

        std::ofstream neigh_file("neighbors.csv");
        neigh_file << "atom_id,shell,neighbor_id\n";
        for (size_t atom = 0; atom < neighbors.size(); ++atom)
            for (size_t shell = 0; shell < geometry.get_shell_count(); ++shell)
                for (int neighbor : neighbors[atom][shell])
                    neigh_file << atom << "," << shell << "," << neighbor
                               << "\n";
        neigh_file.close();

        std::cout << "System size: " << geometry.get_system_size() << "\n";
        std::cout << "Shell count: " << geometry.get_shell_count() << "\n";
        std::cout << "Basis count: " << geometry.get_basis_count() << "\n";
        std::cout << "Total atoms: " << geometry.get_atom_count() << "\n";

        system("python3 ../scripts/visualize_lattice.py");
    } catch (std::exception &e) {
        std::cerr << "Error: " << e.what() << "\n";
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}