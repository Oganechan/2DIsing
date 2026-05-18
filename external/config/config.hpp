#ifndef CONFIG_HPP
#define CONFIG_HPP

#include "nlohmann_json.hpp"
#include <fstream>
#include <string>

using json = nlohmann::json;

class Config {
  public:
    Config() : config_(json::object()) {}

    void load(const std::string &filepath) {
        json new_config;
        std::ifstream file(filepath);

        if (!file.is_open())
            throw std::runtime_error("Config file not found: " + filepath);

        file >> new_config;

        config_ = std::move(new_config);
    }

    template <typename T>
    T get(const std::string &key, const T &default_value = T()) const {
        json current = config_;
        size_t pos = 0;
        std::string k = key;

        while ((pos = k.find('.')) != std::string::npos) {
            std::string part = k.substr(0, pos);

            if (!current.contains(part))
                return default_value;

            current = current[part];
            k = k.substr(pos + 1);
        }

        return current.value(k, default_value);
    }

    json get() { return config_; }

    template <typename T> void set(const std::string &key, const T &value) {
        json *current = &config_;
        size_t pos = 0;
        std::string k = key;

        while ((pos = k.find('.')) != std::string::npos) {
            std::string part = k.substr(0, pos);

            if (!current->contains(part))
                (*current)[part] = json::object();

            current = &(*current)[part];
            k = k.substr(pos + 1);
        }

        (*current)[k] = value;
    }

  private:
    json config_;
};

#endif // CONFIG_HPP