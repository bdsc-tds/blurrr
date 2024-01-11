#ifndef BINXENIUM_ARRAY_SPOTS_HPP
#define BINXENIUM_ARRAY_SPOTS_HPP

#include <set>
#include <vector>

class ArraySpots {
private:
    int array_idx;
    std::set<int> spot_idx;
    double min_loc = static_cast<double>(INT32_MAX);
    double max_loc = 0;

public:
    // constructors
    ArraySpots() = delete; // disable the default constructor without any arguments

    ArraySpots(const ArraySpots &) = default; // default copy constructor
    ArraySpots(ArraySpots &&) = default; // default move constructor
    ArraySpots & operator=(const ArraySpots &) = default; // default copy-assignment constructor
    ArraySpots & operator=(ArraySpots &&) = default; // default move-assignment constructor

    ArraySpots(int);
    ArraySpots(int, int);
    ArraySpots(int, const std::set<int> &);
    ArraySpots(int, const std::vector<int> &);

    // getters
    int get_array_idx() const;
    const std::set<int> & get_spot_idx() const;
    double get_min_loc() const;
    double get_max_loc() const;

    // setters
    void set_min_loc(double);
    void set_max_loc(double);

    // operator overloading
    friend bool operator==(const ArraySpots &, const ArraySpots &);
    friend bool operator<(const ArraySpots &, const ArraySpots &);

    // others
    void add_spot_idx(int);
};

ArraySpots::ArraySpots(int array_idx): array_idx(array_idx) {}

ArraySpots::ArraySpots(int array_idx, int idx): array_idx(array_idx) {
    this->spot_idx.emplace(idx);
}

ArraySpots::ArraySpots(int array_idx, const std::set<int> &idx): array_idx(array_idx) {
    this->spot_idx.insert(idx.begin(), idx.end());
}

ArraySpots::ArraySpots(int array_idx, const std::vector<int> &idx): array_idx(array_idx) {
    this->spot_idx.insert(idx.begin(), idx.end());
}

int
ArraySpots::get_array_idx() const {
    return this->array_idx;
}

const std::set<int> &
ArraySpots::get_spot_idx() const {
    return this->spot_idx;
}

double
ArraySpots::get_min_loc() const {
    return this->min_loc;
}

double
ArraySpots::get_max_loc() const {
    return this->max_loc;
}

void
ArraySpots::set_min_loc(double loc) {
    if (loc < this->min_loc) {
        this->min_loc = loc;
    }
}

void
ArraySpots::set_max_loc(double loc) {
    if (loc > this->max_loc) {
        this->max_loc = loc;
    }
}

bool
operator==(const ArraySpots &l, const ArraySpots &r) {
    if (l.array_idx == r.array_idx) return true;

    return false;
}

bool
operator<(const ArraySpots &l, const ArraySpots &r) {
    if (l.array_idx < r.array_idx) return true;

    return false;
}

void
ArraySpots::add_spot_idx(int idx) {
    this->spot_idx.emplace(idx);
}

#endif // BINXENIUM_ARRAY_SPOTS_HPP
