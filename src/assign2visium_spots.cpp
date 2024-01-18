// [[Rcpp::plugins("cpp11" openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <map>
#include <indicators/cursor_control.hpp>
#include <indicators/progress_bar.hpp>

#include "array_spots.hpp"
#include "assignment.hpp"
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
    const double spot_radius,
    const int thread_num = 1,
    const bool verbose = true
) {
    std::vector<int> thread_hits;

#ifdef _OPENMP
    omp_set_max_active_levels(2);
    omp_set_num_threads(thread_num);

    for (int i = 0; i < thread_num; i++)
        thread_hits.emplace_back(0);

    if (verbose) {
        std::cout << "[DEBUG] The number of threads is " << thread_num << std::endl;
    }
#endif

    Assignment<arma::uword, arma::uword, arma::uword> assignment;

    // build index of spots
    std::vector<ArraySpots> array_cols;
    std::vector<ArraySpots> array_rows;
    create_array_spots(array_coords, img_coords, array_cols, array_rows);

    // progress bar
    indicators::show_console_cursor(false);
    indicators::ProgressBar pb{
        indicators::option::MaxProgress{mole_coords.n_rows - 1},
        indicators::option::BarWidth{50},
        indicators::option::Start{" ["},
        indicators::option::Fill{"█"},
        indicators::option::Lead{"█"},
        indicators::option::Remainder{"-"},
        indicators::option::End{"]"},
        indicators::option::PrefixText{"Binning to spots"},
        indicators::option::ForegroundColor{indicators::Color::blue},
        indicators::option::ShowElapsedTime{true},
        indicators::option::ShowRemainingTime{true},
        indicators::option::FontStyles{
            std::vector<indicators::FontStyle>{indicators::FontStyle::bold}
        }
    };

#pragma omp declare reduction(red_assign:Assignment<arma::uword, arma::uword, arma::uword>:omp_out += omp_in) initializer(omp_priv = omp_orig)

#pragma omp parallel for shared(pb, mole_coords, array_cols, array_rows, img_coords, spot_radius) reduction(red_assign:assignment)
    for (arma::uword mole_idx = 0; mole_idx < mole_coords.n_rows; mole_idx++) {

#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        pb.tick();

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
                assignment.add_assigned(
                    __assigned_idx_moles[0],
                    __assigned_idx_spots[0]
                );

                assignment.add_assign_to_count(
                    __assigned_idx_spots[0],
                    1
                );
            } else if (__assigned_idx_moles.size() > 1) {
                assignment.add_ambi_assigned(
                    __assigned_idx_moles,
                    __assigned_idx_spots
                );
            }
        }
    }

    indicators::show_console_cursor(true);

#ifdef _OPENMP
    if (verbose) {
        print_thread_hits(thread_hits);
    }
#endif

    return Rcpp::List::create(
        Rcpp::_["assignment2Spots"] = convert2arma_assigned(assignment),
        Rcpp::_["ambiAssignment2Spots"] = convert2arma_ambi_assigned(assignment),
        Rcpp::_["countOfSpots"] = convert2arma_assign_to_count(assignment),
        Rcpp::_["propAssigned"] = static_cast<double>(assignment.get_assigned_num()) / static_cast<double>(mole_coords.n_rows)
    );
}
