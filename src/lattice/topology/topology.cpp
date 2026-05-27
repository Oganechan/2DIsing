#include "topology.hpp"

namespace
{

    lattice::CrystalType parse_crystal_type(const std::string &type)
    {
        const std::unordered_map<std::string, lattice::CrystalType> mapping = {
            {"rectangular", lattice::CrystalType::RECTANGULAR},
            {"rectangular_anisotropic", lattice::CrystalType::RECTANGULAR_ANISOTROPIC},
            {"triangular", lattice::CrystalType::TRIANGULAR}};

        if (auto it = mapping.find(type); it != mapping.end())
            return it->second;

        throw std::invalid_argument("Unknown crystal type: " + type);
    }

    std::array<int, 3> get_cell_index(int atom_id, int basis_count,
                                      int system_size)
    {
        int atom_in_cell_id = atom_id % basis_count;
        int temp = atom_id / basis_count;
        int cell_i = temp % system_size;
        int cell_j = temp / system_size;

        return {cell_i, cell_j, atom_in_cell_id};
    }

    std::array<double, 2>
    compute_atom_position(int cell_i, int cell_j, int incel,
                          const std::vector<std::array<double, 2>> &lattice_vectors,
                          const std::vector<std::array<double, 2>> &basis_vectors,
                          int system_size)
    {

        double x = cell_i * lattice_vectors[0][0] + cell_j * lattice_vectors[1][0] +
                   basis_vectors[incel][0];
        double y = cell_i * lattice_vectors[0][1] + cell_j * lattice_vectors[1][1] +
                   basis_vectors[incel][1];

        x = std::fmod(x, system_size * lattice_vectors[0][0]);
        y = std::fmod(y, system_size * lattice_vectors[1][1]);

        return {x, y};
    }

}

lattice::Topology::Topology(const Config &config)
    : system_size_(config.get<int>("lattice.system_size")),
      shell_count_(config.get<int>("lattice.shell_count",
                                   config.get<std::vector<double>>("physical.exchange_constants").size())),
      crystal_type_(
          parse_crystal_type(config.get<std::string>("lattice.crystal_type")))
{
    initialize_lattice();
    compute_indices();
    compute_positions();
}

void lattice::Topology::initialize_lattice()
{
    constexpr double SQRT_3 = 1.7320508075;

    switch (crystal_type_)
    {
    case lattice::CrystalType::RECTANGULAR:
        basis_vectors_ = {{0.0, 0.0}};
        lattice_vectors_ = {{1.0, 0.0}, {0.0, 1.0}};
        directions_ = {{1, 0}, {0, 1}};
        break;

    case lattice::CrystalType::RECTANGULAR_ANISOTROPIC:
        basis_vectors_ = {{0.0, 0.0}};
        lattice_vectors_ = {{1.0, 0.0}, {0.0, 1.1}};
        directions_ = {{1, 0}, {0, 1}};
        break;

    case lattice::CrystalType::TRIANGULAR:
        basis_vectors_ = {{0.0, 0.0}};
        lattice_vectors_ = {{1.0, 0.0}, {0.5, SQRT_3 / 2}};
        directions_ = {{1, 0}, {0, 1}, {-1, 1}};
        break;
    }

    basis_count_ = basis_vectors_.size();
    atom_count_ = system_size_ * system_size_ * basis_count_;
}

void lattice::Topology::compute_indices()
{
    atom_indices_.resize(atom_count_);
    for (size_t i = 0; i < atom_count_; ++i)
        atom_indices_[i] = get_cell_index(i, basis_count_, system_size_);
}

void lattice::Topology::compute_positions()
{

    atom_positions_.resize(atom_count_);

    for (size_t i = 0; i < atom_count_; ++i)
    {
        auto [cell_i, cell_j, incel] = atom_indices_[i];
        atom_positions_[i] =
            compute_atom_position(cell_i, cell_j, incel, lattice_vectors_,
                                  basis_vectors_, system_size_);
    }
}
