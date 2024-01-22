// [[Rcpp::plugins("cpp11" openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <map>
#include <string>

#include "spot_subspots.hpp"
#include "subspots.hpp"
#include "assignment.hpp"
#include "utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


template<typename T>
void
create_spot_subspots(
    const arma::umat &subspot_ids,
    const Rcpp::CharacterVector &subspot_assigned_barcodes,
    const std::map<std::string, T> &spot_map,
    std::map<std::string, SpotSubspots<T, T>> &spot_subspot_map
) {
    if (subspot_ids.n_rows != subspot_assigned_barcodes.length() || subspot_ids.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    if (spot_subspot_map.size() > 0) spot_subspot_map.clear();

    for (T i = 0; i < static_cast<T>(subspot_ids.n_rows); i++) {
        const auto spot_it = spot_map.find(static_cast<std::string>(subspot_assigned_barcodes[i]));

        if (spot_it == spot_map.end()) {
            throw std::invalid_argument("Unknown spot.");
        }

        const auto subspot_it = spot_subspot_map.find(static_cast<std::string>(subspot_assigned_barcodes[i]));

        if (subspot_it == spot_subspot_map.end()) {
            spot_subspot_map.emplace(std::make_pair(
                static_cast<std::string>(subspot_assigned_barcodes[i]),
                SpotSubspots<T, T>(
                    spot_it->second,
                    Subspot<T>(
                        i,
                        subspot_ids(i, 0),
                        subspot_ids(i, 1)
                    )
                )
            ));
        } else {
            (subspot_it->second).add_subspot(
                Subspot<T>(
                    i,
                    subspot_ids(i, 0),
                    subspot_ids(i, 1)
                )
            );
        }
    }
}


// [[Rcpp::export]]
Rcpp::List
assign2visium_subspots(
    const arma::mat &mole_coords,
    const Rcpp::CharacterVector &mole_assigned_barcodes,
    const arma::umat &subspot_ids,
    const arma::mat &img_subspot_coords,
    const Rcpp::CharacterVector &subspot_assigned_barcodes,
    const arma::mat &img_spot_coords,
    const Rcpp::CharacterVector &spot_barcodes,
    const double spot_radius,
    const int thread_num,
    const bool verbose
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

    if (mole_coords.n_rows != mole_assigned_barcodes.length() || mole_coords.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    if (subspot_ids.n_rows != img_subspot_coords.n_rows || subspot_ids.n_rows != subspot_assigned_barcodes.length() || subspot_ids.n_cols != 2 || img_subspot_coords.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    if (img_spot_coords.n_rows != spot_barcodes.length() || img_spot_coords.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    Assignment<arma::uword, arma::uword, arma::uword> assignment;

    // build the index of spots to extract the coordinates in `img_spot_coords`
    std::map<std::string, arma::uword> spot_map;
    build_index(spot_barcodes, spot_map);

    // build a map from spot to subspots
    std::map<std::string, SpotSubspots<arma::uword, arma::uword>> spot_subspot_map;
    create_spot_subspots(
        subspot_ids,
        subspot_assigned_barcodes,
        spot_map,
        spot_subspot_map
    );

    // run time measurement
    double time_start, time_end;

#pragma omp declare reduction(red_assign:Assignment<arma::uword, arma::uword, arma::uword>:omp_out += omp_in) initializer(omp_priv = omp_orig)

#ifdef _OPENMP
    time_start = omp_get_wtime();
#endif

    // assign molecules to subspots
#pragma omp parallel for shared(mole_coords, spot_map, spot_subspot_map, mole_assigned_barcodes, img_spot_coords, spot_radius) reduction(red_assign:assignment)
    for (arma::uword idx = 0; idx < mole_coords.n_rows; idx++) {

#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        const auto it_spot = spot_map.find(
            static_cast<std::string>(mole_assigned_barcodes[idx])
        );
        const auto it_spot_subspot = spot_subspot_map.find(
            static_cast<std::string>(mole_assigned_barcodes[idx])
        );

        if (it_spot == spot_map.end() || it_spot_subspot == spot_subspot_map.end()) {
            throw std::invalid_argument("Unknown spot!");
        }

        // coords of spot center
        const auto __spot_coords = img_spot_coords.row(it_spot->second);

        // coords of molecule
        const auto __mole_coords = mole_coords.row(idx);

        // subspots
        const auto __subspots = (it_spot_subspot->second).get_subspots();

        // iterate over subspots and find the assignment
        std::vector<arma::uword> __assigned_idx_moles, __assigned_idx_subspots;
        for (auto it_subspot = __subspots.begin(); it_subspot != __subspots.end(); it_subspot++) {
            const auto __subspot_idx = (it_subspot->second).get_idx();
            const auto __it_third_pt = __subspots.find(
                (it_subspot->second).get_third_pt_id()
            );

            if (__it_third_pt == __subspots.end()) {
                throw std::runtime_error("Cannot find a subspot!");
            }
            const auto __third_pt_idx = (__it_third_pt->second).get_idx();

            if (is_in_subspot(
                __mole_coords,
                __spot_coords,
                img_subspot_coords.row(__subspot_idx),
                img_subspot_coords.row(__third_pt_idx),
                spot_radius
            )) {
                __assigned_idx_moles.emplace_back(idx);
                __assigned_idx_subspots.emplace_back(__subspot_idx);
            }
        }

        if (__assigned_idx_subspots.size() == 0) {
            throw std::runtime_error("Unable to assign to any subspots!");
        } else if (__assigned_idx_subspots.size() == 1) {
            assignment.add_assigned(
                    __assigned_idx_moles[0],
                    __assigned_idx_subspots[0]
                );

                assignment.add_assign_to_count(
                    __assigned_idx_subspots[0],
                    1
                );
        } else if (__assigned_idx_subspots.size() > 1) {
            assignment.add_ambi_assigned(
                    __assigned_idx_moles,
                    __assigned_idx_subspots
                );
        }
    }

#ifdef _OPENMP
    time_end = omp_get_wtime();

    if (verbose) {
        print_thread_hits(thread_hits);
    }
#endif

    if (verbose) {
        std::cout << "[DEBUG] Assigning molecules to subspots takes " << time_end - time_start << " seconds.\n";
        std::cout << std::endl;
    }

    return Rcpp::List::create(
        Rcpp::_["assignment2Subspots"] = convert2arma_assigned(assignment),
        Rcpp::_["ambiAssignment2Subspots"] = convert2arma_ambi_assigned(assignment),
        Rcpp::_["countOfSubspots"] = convert2arma_assign_to_count(assignment)
    );
}
