#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstddef>
#include <random>
#include <stdexcept>

class Random {
  public:
    explicit Random() : engine_(std::random_device{}()) {}
    explicit Random(size_t sd) : engine_(sd) {}

    template <typename T> T uniform_int(T min, T max) {
        static_assert(std::is_integral<T>::value);
        std::uniform_int_distribution<T> dist(min, max);
        return dist(engine_);
    }

    template <typename T> T uniform_real(T min, T max) {
        static_assert(std::is_floating_point<T>::value);
        std::uniform_real_distribution<T> dist(min, max);
        return dist(engine_);
    }

    bool bernoulli(double p) {
        if (p < 0.0 || p > 1.0)
            throw std::out_of_range("");
        std::bernoulli_distribution dist(p);
        return dist(engine_);
    }

    auto &engine() {return engine_;}

  private:
    std::mt19937 engine_;

    Random(const Random &) = delete;
    Random &operator=(const Random &) = delete;
    Random(Random &&) = default;
    Random &operator=(Random &&) = default;
};

inline Random &thread_local_random() {
    thread_local Random rng;
    return rng;
}

inline Random &thread_local_random(size_t sd) {
    thread_local Random rng(sd);
    return rng;
}

#endif // RANDOM_HPP