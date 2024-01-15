#ifndef BINXENIUM_UTILS_HPP
#define BINXENIUM_UTILS_HPP

#include <RcppArmadillo.h>
#include <csignal>
#include <cmath>
#include <vector>
#include <iostream>

#include "assignment.hpp"

static volatile sig_atomic_t early_stop = 0;


template <typename T1, typename T2>
bool
is_in_spot(
    const T1 &mole_coords,
    const T2 &spot_coords,
    const double radius
) {
    if (mole_coords.n_elem != 2 || mole_coords.n_elem != spot_coords.n_elem) {
        throw std::invalid_argument("Invalid arguments.");
    }

    return pow(mole_coords(0) - spot_coords(0), 2) + pow(mole_coords(1) - spot_coords(1), 2) <= pow(radius, 2);
}

template <typename T>
bool
is_in_subspot(
    const T &mole_coords,
    const T &center_coords,
    const T &fan_coords_1,
    const T &fan_coords_2,
    const double radius
) {
    // The first point of the subspot (center) is left.
    const auto diff_m_1 = mole_coords - fan_coords_1;
    const auto diff_c_1 = center_coords - fan_coords_1;
    const auto diff_1_2 = fan_coords_1 - fan_coords_2;

    const auto __m_12 = diff_m_1[1] * diff_1_2[0] - diff_m_1[0] * diff_1_2[1];
    const auto __c_12 = diff_c_1[1] * diff_1_2[0] - diff_c_1[0] * diff_1_2[1];

    if (__m_12 * __c_12 < 0 && !is_in_spot(
        mole_coords,
        center_coords,
        radius
    )) {
        return false;
    }

    // The second point of the subspot (fan_point_1) is left.
    const auto diff_m_c = mole_coords - center_coords;
    const auto diff_1_c = fan_coords_1 - center_coords;
    const auto diff_2_c = fan_coords_2 - center_coords;

    const auto __m_2c = diff_m_c[1] * diff_2_c[0] - diff_m_c[0] * diff_2_c[1];
    const auto __1_2c = diff_1_c[1] * diff_2_c[0] - diff_1_c[0] * diff_2_c[1];

    if (__m_2c * __1_2c < 0) {
        return false;
    }

    // The third point of the subspot (fan_point_2) is left.
    const auto __m_1c = diff_m_c[1] * diff_1_c[0] - diff_m_c[0] * diff_1_c[1];
    const auto __2_1c = diff_2_c[1] * diff_1_c[0] - diff_2_c[0] * diff_1_c[1];

    return __m_1c * __2_1c >= 0;
}

template <typename T>
void
print_thread_hits(const std::vector<T> &arr) {
    if (arr.size() > 0) {
        for (size_t i = 0; i < arr.size(); i++)
        std::cout << "[DEBUG] Thread " << i << " is hit " << arr[i]
                    << " times.\n";
        std::cout << std::endl;
    }
}

static void
sig_handler(int _) {
    (void) _;
    std::cout << "\nStopping..." << std::endl;

    early_stop = 1;
}

template <typename T, typename C>
arma::Mat<T>
convert2arma_assigned(const Assignment<T, T, C> &val, const bool base_1 = true) {
    arma::Mat<T> ret(val.get_assigned_num(), 2);

    if (val.get_assigned_num() > 0) {
        ret.col(0) = arma::Col<T>(val.get_assigned_mole()) + static_cast<T>(base_1);
        ret.col(1) = arma::Col<T>(val.get_assigned_to()) + static_cast<T>(base_1);
    }

    return ret;
}

template <typename T, typename C>
arma::Mat<T>
convert2arma_ambi_assigned(const Assignment<T, T, C> &val, const bool base_1 = true) {
    arma::Mat<T> ret(val.get_ambi_assigned_num(), 2);

    if (val.get_ambi_assigned_num() > 0) {
        ret.col(0) = arma::Col<T>(val.get_ambi_assigned_mole()) + static_cast<T>(base_1);
        ret.col(1) = arma::Col<T>(val.get_ambi_assigned_to()) + static_cast<T>(base_1);
    }

    return ret;
}

template <typename T, typename S>
arma::Mat<T>
convert2arma_assign_to_count(const Assignment<S, T, T> &val, const bool base_1 = true) {
    arma::Mat<T> ret(val.get_assign_to_count_num(), 2);
    const auto assign_to_count = val.get_assign_to_count();

    auto it = assign_to_count.begin();
    T row_idx = 0;
    while (it != assign_to_count.end()) {
        const arma::Row<T> __count(
            std::vector<T>{it->first + static_cast<T>(base_1), it->second}
        );
        ret.row(row_idx) = std::move(__count);

        it++;
        row_idx++;
    }

    return ret;
}

#endif // BINXENIUM_UTILS_HPP
