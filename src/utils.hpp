#ifndef BINXENIUM_UTILS_HPP
#define BINXENIUM_UTILS_HPP

#include <RcppArmadillo.h>
#include <cmath>
#include <vector>
#include <algorithm>

#include "assignment.hpp"


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

template <typename T1, typename T2, typename T3>
void build_index(
    const T1 &vec,
    std::map<T2, T3> &idx
) {
    if (idx.size() > 0) idx.clear();

    for (T3 i = 0; i < static_cast<T3>(vec.size()); i++) {
        const auto it = idx.find(static_cast<T2>(vec[i]));

        if (it == idx.end()) {
            idx.emplace(std::make_pair(
                static_cast<T2>(vec[i]),
                i
            ));
        } else {
            throw std::invalid_argument("Found duplicate elements.");
        }
    }
}

template <typename T>
void
print_thread_hits(const std::vector<T> &arr) {
    if (arr.size() > 0) {
        for (size_t i = 0; i < arr.size(); i++) {
            std::cout << "[DEBUG] Thread " << i << " is hit " << arr[i] << " times.\n";
        }
        
        std::cout << std::endl;
    }
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

template <typename T1, typename T2>
std::map<T1, T2> &
marge_maps(
    std::map<T1, T2> &m1,
    const std::map<T1, T2> &m2
) {
    for (auto it = m2.begin(); it != m2.end(); it++) {
        const auto __it = m1.find(it->first);

        if (__it != m1.end()) {
            throw std::invalid_argument("Found duplicate elements.");
        }

        m1.insert(std::make_pair(
            it->first,
            it->second
        ));
    }

    return m1;
}

template <typename T>
std::vector<T> &
merge_vectors(
    std::vector<T> &v1,
    const std::vector<T> &v2
) {
    v1.insert(v1.end(), v2.begin(), v2.end());

    return v1;
}

template <typename T>
std::vector<T> &
merge_vectors(
    std::vector<T> &v1,
    const std::vector<T> &v2,
    const bool remove_dups
) {
    if (!remove_dups) {
        return merge_vectors(v1, v2);
    }

    for (T i: v2) {
        const auto it = std::find(v1.begin(), v1.end(), i);

        if (it == v1.end()) {
            v1.emplace_back(i);
        }
    }

    return v1;
    
}

template <typename T>
arma::Mat<T> &
merge_arma_mats(
    arma::Mat<T> &m1,
    const arma::Mat<T> &m2
) {
    if (m1.n_rows != m2.n_rows || m1.n_cols != m2.n_cols) {
        throw std::invalid_argument("Unmatched dimmensions of matrices.");
    }

    m1 = m1 + m2;

    return m1;
}


#endif // BINXENIUM_UTILS_HPP
