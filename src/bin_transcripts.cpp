// [[Rcpp::plugins("cpp11" openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <tuple>
#include <map>
#include <string>

#include "utils.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif


template <typename T1, typename T2, typename T3, typename T4, typename T5>
void
build_bridged_idx(
    const T1 &keys,
    const T2 &vals,
    const std::map<T3, T4> &val_idx,
    std::map<T5, T4> &idx
) {
    if (keys.size() != vals.size()) {
        throw std::invalid_argument("Inputs are of different lengths.");
    }

    if (!idx.empty()) idx.clear();

    for (size_t i = 0; i < static_cast<size_t>(keys.size()); i++) {
        const auto it_key = idx.find(static_cast<T5>(keys[i]));
        const auto it_val = val_idx.find(static_cast<T3>(vals[i]));

        if (it_key != idx.end()) {
            throw std::invalid_argument("Found duplicate elements.");
        }

        if (it_val == val_idx.end()) {
            throw std::invalid_argument("Cannot find key in map.");
        }

        idx.emplace(std::make_pair(
            static_cast<T5>(keys[i]),
            it_val->second
        ));
    }
}

template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
void
build_bridged_idx(
    const T1 &keys,
    const T2 &vals,
    const std::map<T3, T4> &val_idx,
    const std::map<T6, T4> &filter_map,
    std::map<T5, T4> &idx
) {
    if (keys.size() != vals.size()) {
        throw std::invalid_argument("Inputs are of different lengths.");
    }

    if (!idx.empty()) idx.clear();

    for (size_t i = 0; i < static_cast<size_t>(keys.size()); i++) {
        const auto it_filter = filter_map.find(static_cast<T6>(keys[i]));
        if (it_filter == filter_map.end()) continue;

        const auto it_key = idx.find(static_cast<T5>(keys[i]));
        const auto it_val = val_idx.find(static_cast<T3>(vals[i]));

        if (it_key != idx.end()) {
            throw std::invalid_argument("Found duplicate elements.");
        }

        if (it_val == val_idx.end()) {
            throw std::invalid_argument("Cannot find key in map.");
        }

        idx.emplace(std::make_pair(
            static_cast<T5>(keys[i]),
            it_val->second
        ));
    }
}


// [[Rcpp::export]]
arma::umat
bin_transcripts(
    const Rcpp::CharacterVector &unique_spot_names,
    const Rcpp::CharacterVector &unique_gene_names,
    const Rcpp::CharacterVector &assignment_mole_names,
    const Rcpp::CharacterVector &assignment_spot_names,
    const Rcpp::CharacterVector &feature_mole_names,
    const Rcpp::CharacterVector &feature_gene_names,
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

    // create a map from spot name to index
    std::map<std::string, arma::uword> unique_spot_map;
    build_index(unique_spot_names, unique_spot_map);
    
    // create a map from gene name to index
    std::map<std::string, arma::uword> unique_gene_map;
    build_index(unique_gene_names, unique_gene_map);

    // variables assigned in parallel
    std::vector<std::tuple<std::string, arma::uword, arma::uword>> mole_assignment_feature;
    arma::umat gene_by_spot_count(
        unique_gene_map.size(),
        unique_spot_map.size(),
        arma::fill::zeros
    );

    // run time measurement
    double time_start_map, time_end_map, time_start_map2spot_feature, time_end_map2spot_feature, time_end_bin, time_end;

    // build maps for assignments and features
    std::map<std::string, arma::uword> assignment_map;
    std::map<std::string, arma::uword> feature_map;
#ifdef _OPENMP
    time_start_map = omp_get_wtime();
#endif
    build_bridged_idx(
        assignment_mole_names,
        assignment_spot_names,
        unique_spot_map,
        assignment_map
    );

    build_bridged_idx(
        feature_mole_names,
        feature_gene_names,
        unique_gene_map,
        assignment_map,
        feature_map
    );
#ifdef _OPENMP
    time_end_map = omp_get_wtime();
#endif

// #pragma omp declare reduction(red_map:std::map<std::string, arma::uword>:omp_out=marge_maps<std::string, arma::uword>(omp_out, omp_in)) initializer(omp_priv = omp_orig)

#pragma omp declare reduction(red_m_a_f:std::vector<std::tuple<std::string, arma::uword, arma::uword>>:omp_out = merge_vectors<std::tuple<std::string, arma::uword, arma::uword>>(omp_out, omp_in)) initializer(omp_priv = omp_orig)

#pragma omp declare reduction(red_mat:arma::Mat<arma::uword>:omp_out = merge_arma_mats<arma::uword>(omp_out, omp_in)) initializer(omp_priv = omp_orig)

#pragma omp parallel shared(assignment_mole_names, assignment_spot_names, feature_mole_names, feature_gene_names, unique_spot_map, unique_gene_map, assignment_map, feature_map)
{

#ifdef _OPENMP
#pragma omp single
{
    time_start_map2spot_feature = omp_get_wtime();
}

#endif
#pragma omp for reduction(red_m_a_f:mole_assignment_feature)
    for (size_t i = 0; i < static_cast<size_t>(assignment_mole_names.size()); i++) {

#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        const auto it_assignment = assignment_map.find(
            static_cast<std::string>(assignment_mole_names[i])
        );
        const auto it_feature = feature_map.find(
            static_cast<std::string>(assignment_mole_names[i])
        );

        if (it_assignment == assignment_map.end() || it_feature == feature_map.end()) {
            throw std::runtime_error("Missing elements.");
        }

        mole_assignment_feature.emplace_back(std::make_tuple(
            it_assignment->first,
            it_assignment->second,
            it_feature->second
        ));
    }

#ifdef _OPENMP
#pragma omp single
{
    time_end_map2spot_feature = omp_get_wtime();
}
#endif

#pragma omp for reduction(red_mat:gene_by_spot_count)
    for (size_t i = 0; i < mole_assignment_feature.size(); i++) {

#ifdef _OPENMP
#pragma omp atomic update
        thread_hits[omp_get_thread_num()]++;
#endif

        gene_by_spot_count(
            std::get<2>(mole_assignment_feature[i]),
            std::get<1>(mole_assignment_feature[i])
        )++;
    }

#ifdef _OPENMP
#pragma omp single
{
    time_end_bin = omp_get_wtime();
}
#endif

}

#ifdef _OPENMP
    time_end = omp_get_wtime();

    if (verbose) {
        print_thread_hits(thread_hits);

        std::cout << "[DEBUG] [UNPARALLELIZED] Mapping molecules to (sub)spots and features takes " << time_end_map - time_start_map << " seconds.\n";
        std::cout << "[DEBUG] [PARALLELIZED] Combining maps of molecules to (sub)spots and features takes " << time_end_map2spot_feature - time_start_map2spot_feature << " seconds.\n";
        std::cout << "[DEBUG] Binning moleculae per (sub)spot takes " << time_end_bin - time_end_map2spot_feature << " seconds.\n";
        std::cout << "[DEBUG] In total it takes " << time_end - time_start_map << " seconds.\n";
        std::cout << std::endl;
    }
#endif
    
    return gene_by_spot_count;
}
