#ifndef CALCULATOR_HPP
#define CALCULATOR_HPP

#include "../lattice/atoms.hpp"
#include <utility>
#include <vector>

namespace physics {

class Calculator {
  public:
    explicit Calculator(const lattice::Atoms &atoms, const Config &config);

    double calculate_total_energy() const;
    double calculate_flip_energy(int atom_id) const;
    double calculate_total_magnetization() const;
    double calculate_total_order() const;

    double calculate_minimal_energy() const { return E_min; }

    auto &get_exchange_constants() const { return exchange_constants_; }

  private:
    const lattice::Atoms &atoms_;
    const lattice::Geometry &geometry_;

    const std::vector<double> exchange_constants_;

    double E_min;
    double precalc_minimal_energy() const;

    std::vector<std::vector<std::pair<int, double>>> weighted_neighbors_;
    void build_weighted_neighbors();
};

} // namespace physics

#endif // CALCULATOR_HPP