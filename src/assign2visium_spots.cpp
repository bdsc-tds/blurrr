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


template <typename T1, typename T2>
void
create_array_spots_1d(
    const arma::umat &array_coords,
    const arma::mat &img_coords,
    std::vector<ArraySpots1D<T1, T2>> &array_cols,
    std::vector<ArraySpots1D<T1, T2>> &array_rows,
    std::vector<int> &thread_hits
) {
#pragma omp single
{
    if (array_coords.n_rows != img_coords.n_rows || array_coords.n_cols != 2 || array_coords.n_cols != img_coords.n_cols) {
        throw std::invalid_argument("Invalid arguments.");
    }
}

#pragma omp declare reduction(red_vec_array_spots:std::vector<ArraySpots1D<T1, T2>>:omp_out = merge_vectors<ArraySpots1D<T1, T2>>(omp_out, omp_in)) initializer(omp_priv = omp_orig)

#pragma omp for reduction(red_vec_array_spots:array_cols, array_rows)
    for (T2 i = 0; i < static_cast<T2>(array_coords.n_rows); i++) {

#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        const arma::urowvec __array_coords_i = array_coords.row(i);
        const arma::rowvec __img_coords_i = img_coords.row(i);

        ArraySpots1D<T1, T2> __col(
            static_cast<T1>(__array_coords_i(0)),
            i
        );
        ArraySpots1D<T1, T2> __row(
            static_cast<T1>(__array_coords_i(1)),
            i
        );

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

#pragma omp single
{
    std::sort(array_cols.begin(), array_cols.end());
    std::sort(array_rows.begin(), array_rows.end());
}
}

template <typename T1, typename T2>
bool
comp_min_loc(
    const ArraySpots1D<T1, T2> &l,
    const double r
) {
    if (static_cast<double>(l.get_min_loc()) < r) return true;

    return false;
}

template <typename T1, typename T2>
bool
comp_max_loc(
    const ArraySpots1D<T1, T2> &l,
    const double r
) {
    if (static_cast<double>(l.get_max_loc()) < r) return true;

    return false;
}

template <typename T1, typename T2>
void
get_candidate_spot_idx(
    const double coord,
    const std::vector<ArraySpots1D<T1, T2>> &array,
    std::vector<T2> &candidate_spots
) {
    if (!candidate_spots.empty()) candidate_spots.clear();

    const auto min_loc_lower_bound = std::lower_bound(
        array.begin(),
        array.end(),
        coord,
        comp_min_loc<T1, T2>
    );

    const auto max_loc_lower_bound = std::lower_bound(
        array.begin(),
        array.end(),
        coord,
        comp_max_loc<T1, T2>
    );

    if (max_loc_lower_bound == array.end()) {
        candidate_spots = merge_vectors(
            candidate_spots,
            (max_loc_lower_bound - 1)->get_spot_idx()
        );

        return;
    }

    candidate_spots = merge_vectors(
        candidate_spots,
        max_loc_lower_bound->get_spot_idx()
    );

    if (*max_loc_lower_bound != *min_loc_lower_bound || max_loc_lower_bound == array.begin()) return;

    candidate_spots = merge_vectors(
        candidate_spots,
        (max_loc_lower_bound - 1)->get_spot_idx(),
        true
    );
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
    omp_set_max_active_levels(1);
    omp_set_num_threads(thread_num);

    for (int i = 0; i < thread_num; i++)
        thread_hits.emplace_back(0);

    if (verbose) {
        std::cout << "[DEBUG] The number of threads is " << thread_num << std::endl;
    }
#endif

    // variables assigned in parallel
    std::vector<ArraySpots1D<arma::uword, arma::uword>> array_cols;
    std::vector<ArraySpots1D<arma::uword, arma::uword>> array_rows;
    Assignment<arma::uword, arma::uword, arma::uword> assignment;

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
    

    // run time measurement
    double time_start, time_end_array_spots, time_end_assignment;

#pragma omp declare reduction(red_assign:Assignment<arma::uword, arma::uword, arma::uword>:omp_out += omp_in) initializer(omp_priv = omp_orig)

#pragma omp parallel shared(mole_coords, array_coords, img_coords, spot_radius)
{

#ifdef _OPENMP
#pragma omp single
{
    time_start = omp_get_wtime();
}
#endif

    create_array_spots_1d(
        array_coords,
        img_coords,
        array_cols,
        array_rows,
        thread_hits
    );

#ifdef _OPENMP
#pragma omp single
{
    time_end_array_spots = omp_get_wtime();
}
#endif

#pragma omp for reduction(red_assign:assignment)
    for (arma::uword mole_idx = 0; mole_idx < mole_coords.n_rows; mole_idx++) {

#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        pb.tick();

        std::vector<arma::uword> candidate_spots_col;
        std::vector<arma::uword> candidate_spots_row;

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
            candidate_spots_col = merge_vectors(
                candidate_spots_col,
                candidate_spots_row,
                true
            );

            std::vector<arma::uword> __assigned_idx_moles, __assigned_idx_spots;
            for (auto it = candidate_spots_col.begin(); it != candidate_spots_col.end(); it++) {
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
            } else if (__assigned_idx_moles.size() > 1) {
                assignment.add_ambi_assigned(
                    __assigned_idx_moles,
                    __assigned_idx_spots
                );
            }
        }
    }

#ifdef _OPENMP
#pragma omp single
{
    time_end_assignment = omp_get_wtime();
}
#endif

}

    indicators::show_console_cursor(true);

#ifdef _OPENMP
    if (verbose) {
        print_thread_hits(thread_hits);
    }
#endif

    if (verbose) {
        std::cout << "[DEBUG] Mapping spots to arrays takes " << time_end_array_spots - time_start << " seconds.\n";
        std::cout << "[DEBUG] Assigning molecules to spots takes " << time_end_assignment - time_end_array_spots << " seconds.\n";
        std::cout << "[DEBUG] In total it takes " << time_end_assignment - time_start << " seconds.\n";
    }

    return Rcpp::List::create(
        Rcpp::_["assignment2Spots"] = convert2arma_assigned(assignment),
        Rcpp::_["ambiAssignment2Spots"] = convert2arma_ambi_assigned(assignment),
        Rcpp::_["countOfSpots"] = convert2arma_assign_to_count(assignment),
        Rcpp::_["propAssigned"] = static_cast<double>(assignment.get_assigned_num()) / static_cast<double>(mole_coords.n_rows)
    );
}
