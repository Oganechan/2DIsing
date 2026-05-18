#ifndef SPIN_HPP
#define SPIN_HPP

#include "../../external/random/random.hpp"
#include <cstdint>

namespace lattice {

class IsingSpin {
  public:
    IsingSpin() : value_(1) {}
    explicit IsingSpin(int8_t value) : value_(value) {}

    int8_t value() const { return value_; }

    void randomize() { value_ = thread_local_random().bernoulli(0.5) ? 1 : -1; }
    void flip() { value_ = -value_; }
    void set(int8_t new_value) { value_ = new_value; }
    int8_t dot(const IsingSpin &other) const { return value_ * other.value_; }
    int8_t dot(const int8_t other) const { return value_ * other; }

  private:
    int8_t value_;
};

} // namespace lattice

#endif // SPIN_HPP