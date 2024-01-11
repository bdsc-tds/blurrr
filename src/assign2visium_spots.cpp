// [[Rcpp::plugins("cpp11" openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <map>
// #include <indicators/cursor_control.hpp>
// #include <indicators/progress_bar.hpp>

#include "array_spots.hpp"
#include "utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


void
create_array_spots(
    const arma::umat &array_coords,
    const arma::mat &img_coords,
    std::vector<ArraySpots> &array_cols,
    std::vector<ArraySpots> &array_rows
) {
    if (array_coords.n_rows != img_coords.n_rows || array_coords.n_cols != 2 || array_coords.n_cols != img_coords.n_cols) {
        throw std::invalid_argument("Invalid arguments.");
    }

    for (int i = 0; i < static_cast<int>(array_coords.n_rows); i++) {
        const arma::urowvec __array_coords_i = array_coords.row(i);
        const arma::rowvec __img_coords_i = img_coords.row(i);

        ArraySpots __col = ArraySpots(__array_coords_i(0), i);
        ArraySpots __row = ArraySpots(__array_coords_i(1), i);

        const auto __array_cols_loc = std::find(
            array_cols.begin(),
            array_cols.end(),
            __col
        );
        const auto __array_rows_loc = std::find(
            array_rows.begin(),
            array_rows.end(),
            __row
        );

        if (__array_cols_loc == array_cols.end()) {
            __col.set_min_loc(__img_coords_i(0));
            __col.set_max_loc(__img_coords_i(0));

            array_cols.emplace_back(__col);
        } else {
            __array_cols_loc->set_min_loc(__img_coords_i(0));
            __array_cols_loc->set_max_loc(__img_coords_i(0));

            __array_cols_loc->add_spot_idx(i);
        }

        if (__array_rows_loc == array_rows.end()) {
            __row.set_min_loc(__img_coords_i(1));
            __row.set_max_loc(__img_coords_i(1));

            array_rows.emplace_back(__row);
        } else {
            __array_rows_loc->set_min_loc(__img_coords_i(1));
            __array_rows_loc->set_max_loc(__img_coords_i(1));

            __array_rows_loc->add_spot_idx(i);
        }
    }

    std::sort(array_cols.begin(), array_cols.end());
    std::sort(array_rows.begin(), array_rows.end());
}


template <typename T>
void
get_candidate_spot_idx(
    const double coord,
    const std::vector<ArraySpots> &array,
    std::set<T> &candidate_spots
) {
    if (!candidate_spots.empty()) candidate_spots.clear();

    auto left_it = array.begin();
    auto right_it = array.begin() + 1;

    while(right_it != array.end()) {
        if ((coord >= left_it->get_min_loc() && coord <= right_it->get_max_loc()) || (coord >= right_it->get_min_loc() && coord <= left_it->get_max_loc())) {
            auto __it = left_it->get_spot_idx();
            candidate_spots.insert(
                __it.begin(),
                __it.end()
            );

            __it = right_it->get_spot_idx();
            candidate_spots.insert(
                __it.begin(),
                __it.end()
            );
        }

        left_it++;
        right_it++;
    }
}


// [[Rcpp::export]]
Rcpp::List
assign2visium_spots(
    const arma::mat &mole_coords,
    const arma::umat &array_coords,
    const arma::mat &img_coords,
    const double spot_radius
) {
    std::vector<arma::uword> assigned_idx_moles, assigned_idx_spots, ambi_assigned_idx_moles, ambi_assigned_idx_spots;
    std::map<arma::uword, arma::uword> assigned_spot_counts;

    // build index of spots
    std::vector<ArraySpots> array_cols;
    std::vector<ArraySpots> array_rows;
    create_array_spots(array_coords, img_coords, array_cols, array_rows);

    for (arma::uword mole_idx = 0; mole_idx < mole_coords.n_rows; mole_idx++) {
        std::set<arma::uword> candidate_spots_col;
        std::set<arma::uword> candidate_spots_row;

        get_candidate_spot_idx(
            mole_coords(mole_idx, 0),
            array_cols,
            candidate_spots_col
        );

        get_candidate_spot_idx(
            mole_coords(mole_idx, 1),
            array_rows,
            candidate_spots_row
        );

        if (!candidate_spots_col.empty() && !candidate_spots_row.empty()) {
            std::set<arma::uword> candidate_spots;
            std::set_intersection(
                candidate_spots_col.begin(),
                candidate_spots_col.end(),
                candidate_spots_row.begin(),
                candidate_spots_row.end(),
                std::inserter(candidate_spots, candidate_spots.begin())
            );

            std::vector<arma::uword> __assigned_idx_moles, __assigned_idx_spots;
            for (auto it = candidate_spots.begin(); it != candidate_spots.end(); it++) {
                if (is_in_spot(
                    mole_coords.row(mole_idx),
                    img_coords.row(*it),
                    spot_radius
                )) {
                    __assigned_idx_moles.emplace_back(mole_idx);
                    __assigned_idx_spots.emplace_back(*it);
                }
            }

            if (__assigned_idx_moles.size() == 1) {
                const auto it = assigned_spot_counts.find(__assigned_idx_spots[0]);
                if (it != assigned_spot_counts.end()) {
                    (it->second)++;
                } else {
                    assigned_spot_counts.emplace(std::make_pair(__assigned_idx_spots[0], 1));
                }

                assigned_idx_moles.emplace_back(std::move(__assigned_idx_moles[0]));
                assigned_idx_spots.emplace_back(std::move(__assigned_idx_spots[0]));
            } else if (__assigned_idx_moles.size() > 1) {
                ambi_assigned_idx_moles.insert(
                    ambi_assigned_idx_moles.end(),
                    __assigned_idx_moles.begin(),
                    __assigned_idx_moles.end()
                );
                ambi_assigned_idx_spots.insert(
                    ambi_assigned_idx_spots.end(),
                    __assigned_idx_spots.begin(),
                    __assigned_idx_spots.end()
                );
            }
        }
    }

    // adjust to the R index and convert to matrix for return
    if (assigned_idx_moles.size() != assigned_idx_spots.size()) {
        throw std::runtime_error("Internal error.");
    }
    arma::umat assignment2spots(assigned_idx_moles.size(), 2);

    if (assigned_idx_moles.size() > 0) {
        assignment2spots.col(0) = arma::uvec(assigned_idx_moles) + 1;
        assignment2spots.col(1) = arma::uvec(assigned_idx_spots) + 1;
    }

    if (ambi_assigned_idx_moles.size() != ambi_assigned_idx_spots.size()) {
        throw std::runtime_error("Internal error.");
    }
    arma::umat ambi_assignment2spots(ambi_assigned_idx_moles.size(), 2);

    if (ambi_assigned_idx_moles.size() > 0) {
        ambi_assignment2spots.col(0) = arma::uvec(ambi_assigned_idx_moles) + 1;
        ambi_assignment2spots.col(1) = arma::uvec(ambi_assigned_idx_spots) + 1;
    }

    arma::umat count_of_spots(assigned_spot_counts.size(), 2);
    auto it = assigned_spot_counts.begin();
    arma::uword row_idx = 0;
    while (it != assigned_spot_counts.end()) {
        const arma::urowvec __count_of_spots(
            std::vector<arma::uword>{it->first + 1, it->second}
        );
        count_of_spots.row(row_idx) = std::move(__count_of_spots);

        it++;
        row_idx++;
    }

    return Rcpp::List::create(
        Rcpp::_["assignment2Spots"] = assignment2spots,
        Rcpp::_["ambiAssignment2Spots"] = ambi_assignment2spots,
        Rcpp::_["countOfSpots"] = count_of_spots,
        Rcpp::_["propAssigned"] = static_cast<double>(assigned_idx_moles.size()) / static_cast<double>(mole_coords.n_rows)
    );
}
