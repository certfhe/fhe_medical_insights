#include "cohort-io.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>

namespace {

const double NA = std::numeric_limits<double>::quiet_NaN();

std::string trim(const std::string& s)
{
    size_t a = s.find_first_not_of(" \t\r\n\"");
    if (a == std::string::npos) return "";
    size_t b = s.find_last_not_of(" \t\r\n\"");
    return s.substr(a, b - a + 1);
}

std::vector<std::string> splitCsvLine(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    for (char ch : line) {
        if (ch == ',') { out.push_back(trim(cur)); cur.clear(); }
        else cur += ch;
    }
    out.push_back(trim(cur));
    return out;
}

double toDouble(const std::string& s)
{
    if (s.empty() || s == "NA" || s == "nan" || s == "NaN") return NA;
    std::string t = s;
    std::replace(t.begin(), t.end(), ',', '.');   // accept a decimal comma
    return std::stod(t);
}

// min-max scaling to [lo, hi] (r and the ranks do not change)
std::vector<double> normalise(const std::vector<double>& x, double lo, double hi)
{
    double mn = *std::min_element(x.begin(), x.end());
    double mx = *std::max_element(x.begin(), x.end());
    std::vector<double> out(x.size());
    for (size_t i = 0; i < x.size(); ++i)
        out[i] = (mx > mn) ? lo + (hi - lo) * (x[i] - mn) / (mx - mn) : 0.5 * (lo + hi);
    return out;
}

Analysis makeAnalysis(const std::string& name, const std::string& kind,
                      const std::vector<double>& xRaw, const std::vector<double>& yRaw,
                      const std::vector<int>& centerRaw)
{
    // complete cases only
    Analysis a;
    a.name = name;
    a.kind = kind;
    std::vector<double> x, y;
    for (size_t i = 0; i < xRaw.size(); ++i) {
        if (std::isnan(xRaw[i]) || std::isnan(yRaw[i])) continue;
        x.push_back(xRaw[i]);
        y.push_back(yRaw[i]);
        a.center.push_back(centerRaw[i]);
    }
    a.nEff = x.size();
    if (kind == "pearson") {
        a.x = normalise(x, 0.0, 1.0);
        a.y = normalise(y, 0.0, 1.0);
    } else if (kind == "wilcoxon") {
        a.x = normalise(x, 0.2, 1.0);   // above the padding value 0
        a.y = y;                        // 0/1 group label
    } else {
        a.x = x;                        // 0/1 technique
        a.y = y;                        // 0/1 label
    }
    return a;
}

} // namespace

size_t nextPow2(size_t n)
{
    size_t p = 1;
    while (p < n) p <<= 1;
    return std::max<size_t>(p, 4);
}

std::vector<double> padTo(const std::vector<double>& v, size_t slots, double fill)
{
    if (v.size() > slots) throw std::runtime_error("padTo: vector longer than slots");
    std::vector<double> out(slots, fill);
    std::copy(v.begin(), v.end(), out.begin());
    return out;
}

std::vector<double> indicatorVector(size_t nEff, size_t slots)
{
    std::vector<double> ind(slots, 0.0);
    for (size_t i = 0; i < nEff && i < slots; ++i) ind[i] = 1.0;
    return ind;
}

std::vector<std::vector<double>> splitByCenter(const std::vector<double>& v,
                                               const std::vector<int>& center,
                                               size_t numCenters, size_t slots)
{
    // center c gets its own rows, zeros elsewhere
    std::vector<std::vector<double>> parts(numCenters, std::vector<double>(slots, 0.0));
    for (size_t i = 0; i < v.size(); ++i) {
        size_t c = static_cast<size_t>(center[i]) % numCenters;
        parts[c][i] = v[i];
    }
    return parts;
}

Cohort loadCohortCsv(const std::string& path, size_t numCenters)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open CSV file: " + path);

    std::string line;
    if (!std::getline(in, line)) throw std::runtime_error("empty CSV file: " + path);
    if (!line.empty() && line[0] == '\xEF') line = line.substr(3);   // UTF-8 BOM
    std::vector<std::string> header = splitCsvLine(line);
    std::map<std::string, size_t> col;
    for (size_t i = 0; i < header.size(); ++i) col[header[i]] = i;

    const char* required[] = {"center", "pelvic_dmean_gy", "cervical_dmean_gy", "alc_week1", "anc_week1",
                              "lymphopenia_grade_max", "neutropenia_grade_max", "technique_vmat_ht"};
    for (const char* c : required)
        if (!col.count(c)) throw std::runtime_error(std::string("CSV column missing: ") + c);

    std::vector<double> pelvic, cervical, alc, anc, lyGrade, neGrade, tech;
    std::vector<int> center;
    size_t row = 0;
    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        std::vector<std::string> f = splitCsvLine(line);
        if (f.size() < header.size()) f.resize(header.size());
        pelvic.push_back(toDouble(f[col["pelvic_dmean_gy"]]));
        cervical.push_back(toDouble(f[col["cervical_dmean_gy"]]));
        alc.push_back(toDouble(f[col["alc_week1"]]));
        anc.push_back(toDouble(f[col["anc_week1"]]));
        lyGrade.push_back(toDouble(f[col["lymphopenia_grade_max"]]));
        neGrade.push_back(toDouble(f[col["neutropenia_grade_max"]]));
        tech.push_back(toDouble(f[col["technique_vmat_ht"]]));
        int c = numCenters > 0 ? static_cast<int>(row % numCenters)
                               : static_cast<int>(std::lround(toDouble(f[col["center"]])));
        center.push_back(c);
        ++row;
    }

    Cohort cohort;
    cohort.numCenters = numCenters > 0 ? numCenters
                                       : static_cast<size_t>(*std::max_element(center.begin(), center.end())) + 1;

    auto ge3 = [](const std::vector<double>& g) {
        std::vector<double> out(g.size());
        for (size_t i = 0; i < g.size(); ++i) out[i] = std::isnan(g[i]) ? NA : (g[i] >= 3 ? 1.0 : 0.0);
        return out;
    };
    auto eq4 = [](const std::vector<double>& g) {
        std::vector<double> out(g.size());
        for (size_t i = 0; i < g.size(); ++i) out[i] = std::isnan(g[i]) ? NA : (g[i] == 4 ? 1.0 : 0.0);
        return out;
    };

    cohort.analyses.push_back(makeAnalysis("Pearson: pelvic Dmean vs ALC week 1", "pearson", pelvic, alc, center));
    cohort.analyses.push_back(makeAnalysis("Pearson: cervical Dmean vs ANC week 1", "pearson", cervical, anc, center));
    cohort.analyses.push_back(makeAnalysis("Wilcoxon: pelvic Dmean, grade>=3 lymphopenia vs lower", "wilcoxon", pelvic, ge3(lyGrade), center));
    cohort.analyses.push_back(makeAnalysis("Wilcoxon: cervical Dmean, grade>=3 neutropenia vs lower", "wilcoxon", cervical, ge3(neGrade), center));
    cohort.analyses.push_back(makeAnalysis("Chi-squared: technique (VMAT/HT) vs grade 4 lymphopenia", "chi2", tech, eq4(lyGrade), center));
    return cohort;
}

Cohort syntheticCohort(size_t n, size_t numCenters, unsigned seed)
{
    // Synthetic data with the same structure as the real cohort: dose, blood count
    // (lower for higher doses), toxicity group, technique and a grade 4 label.
    // Only used for the benchmarks.
    std::mt19937_64 rng(seed + static_cast<unsigned>(n));
    std::uniform_real_distribution<double> U(0.0, 1.0);
    std::normal_distribution<double> N(0.0, 0.18);
    std::normal_distribution<double> G(0.0, 0.20);

    // the effects are strong enough that no statistic is close to zero
    // (a relative error of a value close to zero says nothing)
    std::vector<double> dose(n), count(n), grp(n), tech(n), lab(n);
    std::vector<int> center(n);
    for (size_t i = 0; i < n; ++i) {
        dose[i] = U(rng);
        count[i] = std::min(1.0, std::max(0.0, 0.7 - 0.5 * dose[i] + N(rng)));
        grp[i] = (dose[i] + G(rng) > 0.55) ? 1.0 : 0.0;
        tech[i] = (U(rng) < 0.4) ? 1.0 : 0.0;
        lab[i] = (U(rng) < (tech[i] > 0.5 ? 0.65 : 0.3)) ? 1.0 : 0.0;
        center[i] = static_cast<int>(i % std::max<size_t>(numCenters, 1));
    }
    // at n = 4 a random draw can leave a group or a cell empty, so use fixed values
    if (n <= 4) {
        double d4[4] = {0.15, 0.35, 0.65, 0.90}, c4[4] = {0.80, 0.55, 0.40, 0.20};
        double g4[4] = {0, 0, 1, 1}, t4[4] = {0, 1, 0, 1}, l4[4] = {0, 1, 0, 1};
        for (size_t i = 0; i < n; ++i) { dose[i] = d4[i]; count[i] = c4[i]; grp[i] = g4[i]; tech[i] = t4[i]; lab[i] = l4[i]; }
    }
    Cohort cohort;
    cohort.numCenters = std::max<size_t>(numCenters, 1);
    cohort.analyses.push_back(makeAnalysis("Pearson (synthetic)", "pearson", dose, count, center));
    cohort.analyses.push_back(makeAnalysis("Wilcoxon (synthetic)", "wilcoxon", dose, grp, center));
    cohort.analyses.push_back(makeAnalysis("Chi-squared (synthetic)", "chi2", tech, lab, center));
    return cohort;
}

double pearsonPlain(const std::vector<double>& x, const std::vector<double>& y)
{
    const size_t n = x.size();
    double meanX = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double meanY = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double num = 0.0, denX = 0.0, denY = 0.0;
    for (size_t i = 0; i < n; ++i) {
        double dx = x[i] - meanX, dy = y[i] - meanY;
        num += dx * dy; denX += dx * dx; denY += dy * dy;
    }
    return num / (std::sqrt(denX) * std::sqrt(denY));
}

double wilcoxonPlain(const std::vector<double>& x, const std::vector<double>& g)
{
    const size_t n = x.size();
    // average ranks for ties
    std::vector<double> R(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        double r = 0.5;
        for (size_t j = 0; j < n; ++j) {
            if (x[i] == x[j]) r += 0.5;        // includes j == i
            else if (x[i] > x[j]) r += 1.0;
        }
        R[i] = r;
    }
    double S1 = 0.0, n1 = 0.0;
    for (size_t i = 0; i < n; ++i)
        if (g[i] > 0.5) { S1 += R[i]; n1 += 1.0; }
    double n0 = static_cast<double>(n) - n1;
    double D1 = S1 - n1 * (n + 1.0) / 2.0;                 // = U1 - n0 n1 / 2
    double D2 = n0 * n1 * (n + 1.0) / 12.0;
    return D1 / std::sqrt(D2);
}

double chi2Plain(const std::vector<double>& t, const std::vector<double>& l)
{
    double a = 0, b = 0, c = 0, d = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        a += (1 - t[i]) * (1 - l[i]);
        b += (1 - t[i]) * l[i];
        c += t[i] * (1 - l[i]);
        d += t[i] * l[i];
    }
    double N = a + b + c + d;
    return N * (a * d - b * c) * (a * d - b * c) / ((a + b) * (c + d) * (a + c) * (b + d));
}
