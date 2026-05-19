#include "interaction.hpp"

namespace
{

    lattice::BoundaryType parse_boundary_type(const std::string &type)
    {
        const std::unordered_map<std::string, lattice::BoundaryType> mapping = {
            {"hard", lattice::BoundaryType::HARD},
            {"periodic", lattice::BoundaryType::PERIODIC}};

        if (auto it = mapping.find(type); it != mapping.end())
            return it->second;

        throw std::invalid_argument("Unknown boundary type: " + type);
    }

    lattice::InteractionType parse_interaction_type(const std::string &type)
    {
        const std::unordered_map<std::string, lattice::InteractionType> mapping = {
            {"direction", lattice::InteractionType::DIRECTION},
            {"shell", lattice::InteractionType::SHELL}};

        if (auto it = mapping.find(type); it != mapping.end())
            return it->second;

        throw std::invalid_argument("Unknown interaction type: " + type);
    }

    int get_atom_id(int cell_i, int cell_j, int incel, int basis_count,
                    int system_size)
    {
        int taco_i = (cell_i % system_size + system_size) % system_size;
        int taco_j = (cell_j % system_size + system_size) % system_size;

        return (taco_i + taco_j * system_size) * basis_count + incel;
    }

    double
    compute_min_distance(double delta_x, double delta_y,
                         const std::vector<std::array<double, 2>> &lattice_vectors,
                         int system_size, bool periodic)
    {

        double min_distance_sq = delta_x * delta_x + delta_y * delta_y;

        if (periodic)
        {
            for (int n = -1; n <= 1; ++n)
            {
                for (int m = -1; m <= 1; ++m)
                {
                    double taco_dx =
                        (n * lattice_vectors[0][0] + m * lattice_vectors[1][0]) *
                            system_size +
                        delta_x;
                    double taco_dy =
                        (n * lattice_vectors[0][1] + m * lattice_vectors[1][1]) *
                            system_size +
                        delta_y;
                    double distance_sq = taco_dx * taco_dx + taco_dy * taco_dy;
                    if (distance_sq < min_distance_sq)
                    {
                        min_distance_sq = distance_sq;
                    }
                }
            }
        }

        return std::sqrt(min_distance_sq);
    }

}

lattice::Interaction::Interaction(const Topology &topology, const Config &config)
    : topology_(topology),
      boundary_type_(parse_boundary_type(config.get<std::string>("lattice.boundary_conditions"))),
      interaction_type_(parse_interaction_type(config.get<std::string>("lattice.interaction_type")))
{
    compute();
}

void lattice::Interaction::compute()
{
    neighbor_table_.resize(topology_.get_atom_count());

    switch (interaction_type_)
    {
    case lattice::InteractionType::DIRECTION:
        compute_direction();
        break;

    case lattice::InteractionType::SHELL:
        compute_shell();
        break;
    }
}

void lattice::Interaction::compute_shell()
{
    const int system_size = topology_.get_system_size();
    const int atom_count = topology_.get_atom_count();
    const int basis_count = topology_.get_basis_count();
    const int shell_count = topology_.get_shell_count();

    const auto &atom_indices = topology_.get_atom_indices();
    const auto &atom_positions = topology_.get_atom_positions();

    int max_cell_search = shell_count * 2;

    for (int atom_id = 0; atom_id < atom_count; ++atom_id)
    {
        neighbor_table_[atom_id].resize(shell_count);

        std::vector<std::pair<int, double>> atom_neighbors;
        auto [cell_i, cell_j, incel] = atom_indices[atom_id];

        for (int di = -max_cell_search; di <= max_cell_search; ++di)
            for (int dj = -max_cell_search; dj <= max_cell_search; ++dj)
                for (int new_incel = 0; new_incel < basis_count; ++new_incel)
                {
                    int neighbor_id = get_atom_id(cell_i + di, cell_j + dj, new_incel, basis_count, system_size);

                    if (neighbor_id < 0 || neighbor_id >= atom_count ||
                        atom_id == neighbor_id)
                        continue;

                    if (boundary_type_ == BoundaryType::HARD)
                    {
                        auto [ni, nj, pos] = atom_indices[neighbor_id];
                        if (ni < 0 || ni >= system_size || nj < 0 ||
                            nj >= system_size)
                            continue;
                    }

                    double delta_x = atom_positions[neighbor_id][0] -
                                     atom_positions[atom_id][0];
                    double delta_y = atom_positions[neighbor_id][1] -
                                     atom_positions[atom_id][1];

                    double distance = compute_min_distance(
                        delta_x, delta_y, topology_.get_lattice_vectors(), system_size,
                        boundary_type_ == lattice::BoundaryType::PERIODIC);

                    atom_neighbors.emplace_back(neighbor_id, distance);
                }

        std::sort(
            atom_neighbors.begin(), atom_neighbors.end(),
            [](const auto &a, const auto &b)
            { return a.second < b.second; });

        double current_distance = atom_neighbors[0].second;
        int current_shell = 0;

        for (const auto &[neighbor_id, distance] : atom_neighbors)
        {
            if (std::abs(distance - current_distance) > 1.0e-5)
            {
                current_distance = distance;
                ++current_shell;
                if (current_shell >= shell_count)
                    break;
            }
            if (current_shell < shell_count)
                neighbor_table_[atom_id][current_shell].push_back(neighbor_id);
        }

        for (int shell = 0; shell < shell_count; ++shell)
            std::sort(neighbor_table_[atom_id][shell].begin(),
                      neighbor_table_[atom_id][shell].end());
    }
}

void lattice::Interaction::compute_direction()
{
    const int system_size = topology_.get_system_size();
    const int atom_count = topology_.get_atom_count();
    const int basis_count = topology_.get_basis_count();
    const int shell_count = topology_.get_shell_count();

    const auto &atom_indices = topology_.get_atom_indices();
    const auto &directions = topology_.get_directions();

    for (int atom_id = 0; atom_id < atom_count; ++atom_id)
    {
        neighbor_table_[atom_id].resize(shell_count);

        auto [cell_i, cell_j, incel] = atom_indices[atom_id];

        for (int dir = 0; dir < shell_count; ++dir)
        {
            auto [di, dj] = directions[dir];

            for (int sign : {-1, 1})
            {
                int ni = cell_i + sign * di;
                int nj = cell_j + sign * dj;

                if (boundary_type_ == BoundaryType::PERIODIC)
                {
                    ni = (ni + system_size) % system_size;
                    nj = (nj + system_size) % system_size;
                }

                if (boundary_type_ == BoundaryType::HARD)
                    if (ni < 0 || ni >= system_size || nj < 0 || nj >= system_size)
                        continue;

                int neighbor_id = get_atom_id(ni, nj, incel, basis_count, system_size);

                if (neighbor_id < 0 || neighbor_id >= atom_count)
                    continue;

                if (neighbor_id == atom_id)
                    continue;

                neighbor_table_[atom_id][dir].push_back(neighbor_id);
            }
        }
    }
}
