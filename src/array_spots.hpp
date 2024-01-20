#ifndef BINXENIUM_ARRAY_SPOTS_HPP
#define BINXENIUM_ARRAY_SPOTS_HPP

#include <vector>
#include <algorithm>


template <typename T1, typename T2>
class ArraySpots1D {
private:
    T1 array_idx;
    std::vector<T2> spot_idx;
    double min_loc = static_cast<double>(INT32_MAX);
    double max_loc = -static_cast<double>(INT32_MAX);

public:
    // constructors
    ArraySpots1D() = delete; // disable the default constructor without any arguments

    ArraySpots1D(const ArraySpots1D &) = default; // default copy constructor
    ArraySpots1D(ArraySpots1D &&) = default; // default move constructor
    ArraySpots1D & operator=(const ArraySpots1D &) = default; // default copy-assignment constructor
    ArraySpots1D & operator=(ArraySpots1D &&) = default; // default move-assignment constructor

    ArraySpots1D(T1);
    ArraySpots1D(T1, T2);
    ArraySpots1D(T1, const std::vector<T2> &);

    // getters
    T1 get_array_idx() const;
    const std::vector<T2> & get_spot_idx() const;
    double get_min_loc() const;
    double get_max_loc() const;

    // setters
    void set_min_loc(double);
    void set_max_loc(double);

    // operator overloading
    friend bool operator==(
        const ArraySpots1D &l,
        const ArraySpots1D &r
    ) {
        if (l.array_idx == r.array_idx) return true;

        return false;
    }

    friend bool operator!=(
        const ArraySpots1D &l,
        const ArraySpots1D &r
    ) {
        return !(l == r);
    }

    friend bool operator<(
        const ArraySpots1D &l,
        const ArraySpots1D &r
    ) {
        if (l.min_loc < r.min_loc) return true;

        return false;
    }

    // others
    void add_spot_idx(T2);
    void add_spot_idx(const std::vector<T2> &);
};

template <typename T1, typename T2>
ArraySpots1D<T1, T2>::ArraySpots1D(T1 array_idx): array_idx(array_idx) {}

template <typename T1, typename T2>
ArraySpots1D<T1, T2>::ArraySpots1D(T1 array_idx, T2 idx): array_idx(array_idx) {
    this->spot_idx.emplace_back(idx);
}

template <typename T1, typename T2>
ArraySpots1D<T1, T2>::ArraySpots1D(T1 array_idx, const std::vector<T2> &idx): array_idx(array_idx) {
    for (auto it = idx.begin(); it != idx.end(); it++) {
        this->add_spot_idx(*it);
    }
}

template <typename T1, typename T2>
T1
ArraySpots1D<T1, T2>::get_array_idx() const {
    return this->array_idx;
}

template <typename T1, typename T2>
const std::vector<T2> &
ArraySpots1D<T1, T2>::get_spot_idx() const {
    return this->spot_idx;
}

template <typename T1, typename T2>
double
ArraySpots1D<T1, T2>::get_min_loc() const {
    return this->min_loc;
}

template <typename T1, typename T2>
double
ArraySpots1D<T1, T2>::get_max_loc() const {
    return this->max_loc;
}

template <typename T1, typename T2>
void
ArraySpots1D<T1, T2>::set_min_loc(double loc) {
    if (loc < this->min_loc) {
        this->min_loc = loc;
    }
}

template <typename T1, typename T2>
void
ArraySpots1D<T1, T2>::set_max_loc(double loc) {
    if (loc > this->max_loc) {
        this->max_loc = loc;
    }
}

template <typename T1, typename T2>
void
ArraySpots1D<T1, T2>::add_spot_idx(T2 idx) {
    const auto it = std::find(
        this->spot_idx.begin(),
        this->spot_idx.end(),
        idx
    );

    if (it == this->spot_idx.end()) {
        this->spot_idx.emplace_back(idx);
    }
}

template <typename T1, typename T2>
void
ArraySpots1D<T1, T2>::add_spot_idx(const std::vector<T2> &idx) {
    for (auto it = idx.begin(); it != idx.end(); it++) {
        this->add_spot_idx(*it);
    }
}


#endif // BINXENIUM_ARRAY_SPOTS_HPP
