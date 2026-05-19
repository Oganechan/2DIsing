#ifndef INTERACTION_HPP
#define INTERACTION_HPP

#include "topology.hpp"

namespace lattice
{

  enum class BoundaryType
  {
    PERIODIC,
    HARD
  };

  enum class InteractionType
  {
    SHELL,
    DIRECTION
  };

  class Interaction
  {
  public:
    explicit Interaction(const Topology &topology, const Config &config);

    auto &get_neighbor_table() const { return neighbor_table_; }

  private:
    const Topology &topology_;

    const BoundaryType boundary_type_;
    const InteractionType interaction_type_;

    std::vector<std::vector<std::vector<int>>> neighbor_table_;

    void compute();

    void compute_shell();
    void compute_direction();
  };

} // namespace lattice

#endif // INTERACTION_HPP