// [[Rcpp::plugins("cpp11" openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
// #include <indicators/cursor_control.hpp>
// #include <indicators/progress_bar.hpp>

#include "spot_subspots.hpp"
#include "subspots.hpp"
#include "utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


template<typename T>
void
create_spot_map(
    const Rcpp::CharacterVector &barcodes,
    std::map<std::string, T> &spot_map
) {
    if (spot_map.size() > 0) spot_map.clear();

    for (T i = 0; i < static_cast<T>(barcodes.length()); i++) {
        spot_map.emplace(std::make_pair(barcodes[i], i));
    }
}


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
    const double spot_radius
) {
    if (mole_coords.n_rows != mole_assigned_barcodes.length() || mole_coords.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    if (subspot_ids.n_rows != img_subspot_coords.n_rows || subspot_ids.n_rows != subspot_assigned_barcodes.length() || subspot_ids.n_cols != 2 || img_subspot_coords.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    if (img_spot_coords.n_rows != spot_barcodes.length() || img_spot_coords.n_cols != 2) {
        throw std::invalid_argument("Invalid arguments.");
    }

    std::vector<arma::uword> assigned_idx_moles, assigned_idx_subspots, ambi_assigned_idx_moles, ambi_assigned_idx_subspots;
    std::map<arma::uword, arma::uword> assigned_subspot_counts;

    // build the index of spots to extract the coordinates in `img_spot_coords`
    std::map<std::string, arma::uword> spot_map;
    create_spot_map(spot_barcodes, spot_map);

    // build a map from spot to subspots
    std::map<std::string, SpotSubspots<arma::uword, arma::uword>> spot_subspot_map;
    create_spot_subspots(
        subspot_ids,
        subspot_assigned_barcodes,
        spot_map,
        spot_subspot_map
    );

    // assign molecules to subspots
    for (arma::uword idx = 0; idx < mole_coords.n_rows; idx++) {
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
            const auto it = assigned_subspot_counts.find(__assigned_idx_subspots[0]);
            if (it != assigned_subspot_counts.end()) {
                (it->second)++;
            } else {
                assigned_subspot_counts.emplace(
                    std::make_pair(__assigned_idx_subspots[0], 1)
                );
            }

            assigned_idx_moles.emplace_back(std::move(__assigned_idx_moles[0]));
            assigned_idx_subspots.emplace_back(std::move(__assigned_idx_subspots[0]));
        } else if (__assigned_idx_subspots.size() > 1) {
            ambi_assigned_idx_moles.insert(
                ambi_assigned_idx_moles.end(),
                __assigned_idx_moles.begin(),
                __assigned_idx_moles.end()
            );
            ambi_assigned_idx_subspots.insert(
                ambi_assigned_idx_subspots.end(),
                __assigned_idx_subspots.begin(),
                __assigned_idx_subspots.end()
            );
        }
    }

    // adjust to the R index and convert to matrix for return
    if (assigned_idx_moles.size() != assigned_idx_subspots.size()) {
        throw std::runtime_error("Internal error.");
    }
    arma::umat assignment2subspots(assigned_idx_moles.size(), 2);

    if (assigned_idx_moles.size() > 0) {
        assignment2subspots.col(0) = arma::uvec(assigned_idx_moles) + 1;
        assignment2subspots.col(1) = arma::uvec(assigned_idx_subspots) + 1;
    }

    if (ambi_assigned_idx_moles.size() != ambi_assigned_idx_subspots.size()) {
        throw std::invalid_argument("Internal error.");
    }
    arma::umat ambi_assignment2subspots(ambi_assigned_idx_moles.size(), 2);

    if (ambi_assigned_idx_moles.size() > 0) {
        ambi_assignment2subspots.col(0) = arma::uvec(ambi_assigned_idx_moles) + 1;
        ambi_assignment2subspots.col(1) = arma::uvec(ambi_assigned_idx_subspots) + 1;
    }

    arma::umat count_of_subspots(assigned_subspot_counts.size(), 2);
    auto it = assigned_subspot_counts.begin();
    arma::uword row_idx = 0;
    while (it != assigned_subspot_counts.end()) {
        const arma::urowvec __count_of_subspots(
            std::vector<arma::uword>{it->first + 1, it->second}
        );
        count_of_subspots.row(row_idx) = std::move(__count_of_subspots);

        it++;
        row_idx++;
    }

    return Rcpp::List::create(
        Rcpp::_["assignment2Subspots"] = assignment2subspots,
        Rcpp::_["ambiAssignment2Subspots"] = ambi_assignment2subspots,
        Rcpp::_["countOfSubspots"] = count_of_subspots
    );
}
