#pragma once

#include <cstddef>
#include <string>
#include <vector>

// One analysis of Table 1 of the paper (complete cases, no padding yet).
//   "pearson"  : x, y continuous, scaled to [0, 1]
//   "wilcoxon" : x continuous scaled to [0.2, 1], y = group (1 = higher toxicity)
//   "chi2"     : x = technique (1 = VMAT/HT), y = label (1 = grade 4 lymphopenia)
struct Analysis
{
    std::string name;
    std::string kind;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<int> center;   // emulated center of each row
    size_t nEff = 0;
};

struct Cohort
{
    std::vector<Analysis> analyses;
    size_t numCenters = 1;
};

// Reads the cohort CSV (header line, empty field = missing value):
//   patient_index,center,pelvic_dmean_gy,cervical_dmean_gy,alc_week1,anc_week1,
//   lymphopenia_grade_max,neutropenia_grade_max,technique_vmat_ht
// and builds the five analyses of Table 1. If numCenters > 0, row i goes to
// center i mod numCenters, otherwise the center column is used.
Cohort loadCohortCsv(const std::string& path, size_t numCenters = 0);

// Synthetic cohort of size n (fixed seed) with one analysis of each kind, for the benchmarks.
Cohort syntheticCohort(size_t n, size_t numCenters, unsigned seed = 20260925u);

// Helpers
size_t nextPow2(size_t n);
std::vector<double> padTo(const std::vector<double>& v, size_t slots, double fill = 0.0);
std::vector<double> indicatorVector(size_t nEff, size_t slots);   // 1 on the real slots, 0 on the padding
std::vector<std::vector<double>> splitByCenter(const std::vector<double>& v,
                                               const std::vector<int>& center,
                                               size_t numCenters, size_t slots);

// Plaintext versions, to compare with the encrypted results
double pearsonPlain(const std::vector<double>& x, const std::vector<double>& y);
// Wilcoxon rank-sum z (normal approximation, average ranks, no tie or continuity correction)
double wilcoxonPlain(const std::vector<double>& x, const std::vector<double>& g);
// chi-squared for the 2x2 table of t against l (no Yates correction)
double chi2Plain(const std::vector<double>& t, const std::vector<double>& l);
