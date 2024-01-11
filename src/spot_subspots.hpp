#ifndef BINXENIUM_SPOT_SUBSPOTS_HPP
#define BINXENIUM_SPOT_SUBSPOTS_HPP

#include <string>
#include <map>
#include <iostream>

#include "subspots.hpp"

template<typename T, typename K>
class SpotSubspots {
private:
    T spot_id;
    std::map<K, Subspot<K>> subspots; // use `id` as key

public:
    // constructors
    SpotSubspots() = delete; // disable the default constructor without any arguments

    SpotSubspots(const SpotSubspots &) = default; // default copy constructor
    SpotSubspots(SpotSubspots &&) = default; // default move constructor
    SpotSubspots & operator=(const SpotSubspots &) = default; // default copy-assignment constructor
    SpotSubspots & operator=(SpotSubspots &&) = default; // default move-assignment constructor
    
    SpotSubspots(T);
    SpotSubspots(T, const Subspot<K> &);
    SpotSubspots(T, Subspot<K> &&);

    // getters
    T get_spot_id() const;
    const std::map<K, Subspot<K>> & get_subspots() const;
    const Subspot<K> & get_subspot(K) const;

    // others
    void add_subspot(const Subspot<K> &, bool = false);
    void add_subspot(Subspot<K> &&, bool = false);
};

template<typename T, typename K>
SpotSubspots<T, K>::SpotSubspots(T spot_id): spot_id(spot_id) {}

template<typename T, typename K>
SpotSubspots<T, K>::SpotSubspots(T spot_id, const Subspot<K> &subspot): spot_id(spot_id) {
    this->add_subspot(subspot);
}

template<typename T, typename K>
SpotSubspots<T, K>::SpotSubspots(T spot_id, Subspot<K> &&subspot): spot_id(spot_id) {
    this->add_subspot(std::move(subspot));
}

template<typename T, typename K>
T
SpotSubspots<T, K>::get_spot_id() const {
    return this->spot_id;
}

template<typename T, typename K>
const std::map<K, Subspot<K>> &
SpotSubspots<T, K>::get_subspots() const {
    return this->subspots;
}

template<typename T, typename K>
const Subspot<K> &
SpotSubspots<T, K>::get_subspot(K id) const {
    if (this->subspots.size() == 0) {
        throw std::invalid_argument("The map is empty.");
    }

    const auto it = this->subspots.find(id);
    if (it == this->subspots.end()) {
        throw std::invalid_argument("Cannot find key in map.");
    }

    return it->second;
}

template<typename T, typename K>
void
SpotSubspots<T, K>::add_subspot(const Subspot<K> &val, bool force) {
    const auto it = this->subspots.find(val.get_id());
    if (it != this->subspots.end()) {
        if (!force) {
            std::cerr << "Key is already in the map." << std::endl;
            return;
        }

        this->subspots.at(val.get_id()) = val;
    }

    this->subspots.emplace(std::make_pair(val.get_id(), val));
}

template<typename T, typename K>
void
SpotSubspots<T, K>::add_subspot(Subspot<K> &&val, bool force) {
    const auto it = this->subspots.find(val.get_id());
    if (it != this->subspots.end()) {
        if (!force) {
            std::cerr << "Key is already in the map." << std::endl;
            return;
        }

        this->subspots.at(val.get_id()) = std::move(val);
    }

    this->subspots.emplace(std::make_pair(val.get_id(), std::move(val)));
}

#endif // BINXENIUM_SPOT_SUBSPOTS_HPP
