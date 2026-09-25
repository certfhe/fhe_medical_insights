// Single-party FHE version of the three tests (no threshold decryption), with the
// same masked release as mpc-statistics. This is the "1-party (FHE)" row of the
// tables in the paper.
//
// Usage:
//   statistics-test [--n N] [--runs R] [--csv FILE] [--no-wilcoxon] [--bootstrap] [--verbose]
//                   [--depth D] [--intbits I] [--decbits F]
//                   [--chi2-depth E] [--chi2-intbits I] [--chi2-decbits F]
//                   [--compare-depth C] [--rank-intbits I] [--rank-decbits F]
//
// The contexts are the same as in mpc-statistics, so that the FHE and THE rows only
// differ by the threshold part. The ranking context has depth C + 4 here (no
// interactive bootstrapping); --bootstrap uses normal CKKS bootstrapping instead.

#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "statistics.h"
#include "ranking.h"
#include "cohort-io.h"
#include "mask.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <sstream>
#include <string>

using namespace lbcrypto;

template <class F>
auto timed(F&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    auto res = f();
    auto end = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration<double>(end - start).count();
    return std::make_pair(t, res);
}

static bool g_verbose = false;

struct Party {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> keys;
    size_t slots = 0;
    bool bootstrap = false;
};

static Party makeParty(usint intBits, usint decBits, usint depth, size_t matrixSize, size_t slots,
                       bool needMatrixKeys, bool bootstrap) {
    Party p;
    p.slots = slots;
    p.bootstrap = bootstrap;
    std::cout << "Generating single-party context (slots = " << slots << ", depth = " << depth
              << (bootstrap ? " + bootstrapping" : "") << ")..." << std::endl;
    p.cc = generateCryptoContext(intBits, decBits, depth, slots, bootstrap, 0, g_verbose);

    std::set<int32_t> indexSet;
    if (needMatrixKeys)
        for (int32_t idx : getRotationIndices(matrixSize)) indexSet.insert(idx);
    for (size_t i = 1; i < slots; i *= 2) {
        indexSet.insert(static_cast<int32_t>(i));
        indexSet.insert(-static_cast<int32_t>(i));
    }
    std::vector<int32_t> indices(indexSet.begin(), indexSet.end());

    if (g_verbose) {
        // sizes for Table 7: a fresh ciphertext has 2*N*L words,
        // a key-switching key 2*dnum*(L+|P|)*N words
        const auto rnsParams = std::dynamic_pointer_cast<CryptoParametersRNS>(p.cc->GetCryptoParameters());
        const size_t numQ = p.cc->GetCryptoParameters()->GetElementParams()->GetParams().size();
        const size_t numP = (rnsParams->GetParamsP() ? rnsParams->GetParamsP()->GetParams().size() : 0);
        const uint32_t dnum = rnsParams->GetNumPartQ();
        const double mib = 1024.0 * 1024.0;
        const double N = static_cast<double>(p.cc->GetRingDimension());
        const size_t numKeys = 1 + indices.size();   // relinearisation + rotation keys
        std::cout << "Ring dimension / CRT primes in Q / in P / dnum : " << p.cc->GetRingDimension() << " / " << numQ
                  << " / " << numP << " / " << dnum << std::endl
                  << "Fresh ciphertext (2 N L words)                 : " << 2.0 * N * numQ * 8.0 / mib << " MiB" << std::endl
                  << "Evaluation keys (1 relin + " << indices.size() << " rotation)  : " << numKeys << " keys, about "
                  << numKeys * 2.0 * dnum * (numQ + numP) * N * 8.0 / mib << " MiB of public material" << std::endl;
    }

    p.keys = keyGeneration(p.cc, indices, slots, bootstrap, g_verbose);
    return p;
}

static Ciphertext<DCRTPoly> encryptByCenter(const Party& p, const std::vector<double>& values,
                                            const std::vector<int>& center, size_t numCenters) {
    auto parts = splitByCenter(values, center, numCenters, p.slots);
    Ciphertext<DCRTPoly> sum;
    for (size_t c = 0; c < numCenters; ++c) {
        Plaintext pt = p.cc->MakeCKKSPackedPlaintext(parts[c], 1, 0, nullptr, p.slots);
        auto ct = p.cc->Encrypt(p.keys.publicKey, pt);
        sum = (c == 0) ? ct : p.cc->EvalAdd(sum, ct);
    }
    return sum;
}

static Plaintext indicatorPlaintext(const Party& p, size_t nEff) {
    return p.cc->MakeCKKSPackedPlaintext(indicatorVector(nEff, p.slots), 1, 0, nullptr, p.slots);
}

static Ciphertext<DCRTPoly> sumFirst(const CryptoContext<DCRTPoly>& cc, Ciphertext<DCRTPoly> c, size_t length) {
    auto s = c;
    for (size_t i = 1; i < length; i *= 2) s = cc->EvalAdd(s, cc->EvalRotate(s, i));
    return s;
}

static Ciphertext<DCRTPoly> sumAll(const Party& p, Ciphertext<DCRTPoly> c) {
    return sumFirst(p.cc, c, p.slots);
}

static double dec(const Party& p, const Ciphertext<DCRTPoly>& c) {
    return decryptScalar(p.cc, p.keys.secretKey, c);
}

// Pearson correlation (Section 3.4.1)
static double pearsonCtxt(const Party& p, Ciphertext<DCRTPoly> xC, Ciphertext<DCRTPoly> yC,
                          size_t nEff, const Plaintext& indP) {
    auto cc = p.cc;
    auto meanX = cc->EvalMult(sumAll(p, xC), 1.0 / static_cast<double>(nEff));
    auto meanY = cc->EvalMult(sumAll(p, yC), 1.0 / static_cast<double>(nEff));
    auto xCen = cc->EvalMult(cc->EvalSub(xC, meanX), indP);
    auto yCen = cc->EvalMult(cc->EvalSub(yC, meanY), indP);

    auto xyC = sumAll(p, cc->EvalMult(xCen, yCen));
    auto x2C = sumAll(p, cc->EvalMult(xCen, xCen));
    auto y2C = sumAll(p, cc->EvalMult(yCen, yCen));

    // random masks; release covariance and variances (sums divided by n), 1/n cancels in r
    const double u1 = sampleMask(), u2 = sampleMask(), u3 = u1 / u2;
    const double invN = 1.0 / static_cast<double>(nEff);
    const double m1 = dec(p, cc->EvalMult(xyC, u1 * invN));
    const double m2 = dec(p, cc->EvalMult(x2C, u2 * u2 * invN));
    const double m3 = dec(p, cc->EvalMult(y2C, u3 * u3 * invN));
    if (g_verbose) {
        std::ostringstream os;   // separate stream, so the format of cout is not changed
        os << std::setprecision(6) << "    masked release (Pearson): " << m1 << ", " << m2 << ", " << m3
           << "  ->  r = " << m1 / (std::sqrt(m2) * std::sqrt(m3));
        std::cout << os.str() << std::endl;
    }
    return m1 / (std::sqrt(m2) * std::sqrt(m3));   // masks cancel: u1/(u2 u3) = 1
}

// Wilcoxon rank-sum (Section 3.4.2), same padding as in mpc-statistics.cpp
static double wilcoxonCtxt(const Party& p, Ciphertext<DCRTPoly> xC, Ciphertext<DCRTPoly> gC,
                           size_t nEff, size_t matrixSize, usint compareDepth) {
    auto cc = p.cc;
    // every row of the n x n layout holds the rank vector and gC is zero outside
    // the first matrixSize slots, so we do not need to mask the other rows
    Ciphertext<DCRTPoly> rankVec = rank(xC, matrixSize, -1.0, 1.0, depth2degree(compareDepth));
    if (p.bootstrap) rankVec = cc->EvalBootstrap(rankVec);

    auto s1C = sumFirst(cc, cc->EvalMult(rankVec, gC), matrixSize);
    auto n1C = sumFirst(cc, gC, matrixSize);

    const double n = static_cast<double>(nEff);
    const double pad = static_cast<double>(matrixSize - nEff);
    auto d1C = cc->EvalSub(s1C, cc->EvalMult(n1C, pad + (n + 1.0) / 2.0));
    auto d2C = cc->EvalMult(cc->EvalMult(n1C, cc->EvalSub(n, n1C)), (n + 1.0) / 12.0);

    // masked release (u D1/n^2, u^2 D2/n^4), the powers of n cancel in Z
    const double u = sampleMask();
    const double m1 = dec(p, cc->EvalMult(d1C, u / (n * n)));
    const double m2 = dec(p, cc->EvalMult(d2C, u * u / (n * n * n * n)));
    if (g_verbose) {
        std::ostringstream os;
        os << std::setprecision(6) << "    masked release (Wilcoxon): " << m1 << ", " << m2
           << "  ->  z = " << m1 / std::sqrt(m2);
        std::cout << os.str() << std::endl;
    }
    return m1 / std::sqrt(m2);
}

// Chi-squared (Section 3.4.3), closed form with one masked pair
static double chi2Ctxt(const Party& p, Ciphertext<DCRTPoly> tC, Ciphertext<DCRTPoly> lC,
                       size_t nEff, const Plaintext& indP) {
    auto cc = p.cc;
    auto notT = cc->EvalSub(indP, tC);
    auto notL = cc->EvalSub(indP, lC);
    auto aC = sumAll(p, cc->EvalMult(notT, notL));
    auto bC = sumAll(p, cc->EvalMult(notT, lC));
    auto cC = sumAll(p, cc->EvalMult(tC, notL));
    auto dC = sumAll(p, cc->EvalMult(tC, lC));

    // work with the proportions a/N, ... (see chi2MPC in mpc-statistics.cpp)
    const double N = static_cast<double>(nEff);
    const double invN = 1.0 / N;
    auto paC = cc->EvalMult(aC, invN), pbC = cc->EvalMult(bC, invN);
    auto pcC = cc->EvalMult(cC, invN), pdC = cc->EvalMult(dC, invN);
    auto diff = cc->EvalSub(cc->EvalMult(paC, pdC), cc->EvalMult(pbC, pcC));   // level 3
    auto sqC  = cc->EvalMult(diff, diff);                                        // level 4
    auto qC   = cc->EvalMult(cc->EvalMult(cc->EvalAdd(paC, pbC), cc->EvalAdd(pcC, pdC)),
                             cc->EvalMult(cc->EvalAdd(paC, pcC), cc->EvalAdd(pbC, pdC)));   // level 4

    // masked release (v (ad-bc)^2/N^4, v r0 r1 c0 c1/N^4), the ratio is phi^2 = chi2 / N
    const double v = sampleMask();
    const double mP = dec(p, cc->EvalMult(sqC, v));
    const double mQ = dec(p, cc->EvalMult(qC, v));
    if (g_verbose) {
        std::ostringstream os;
        os << std::setprecision(6) << "    masked release (chi2): " << mP << ", " << mQ
           << "  ->  chi2 = " << N * mP / mQ;
        std::cout << os.str() << std::endl;
    }
    return N * mP / mQ;   // chi2 = N * phi^2, the mask cancels
}

struct Options {
    size_t n = 512;
    size_t runs = 50;
    usint depth = 8;         // Pearson context, same as in mpc-statistics
    usint intBits = 1;       // Pearson: 1+35 bits
    usint decBits = 35;
    usint chi2Depth = 5;     // chi-squared context
    usint chi2IntBits = 10;  // chi-squared: 10+50 bits
    usint chi2DecBits = 50;
    usint compareDepth = 8;  // comparison depth of the ranking primitive
    usint rankIntBits = 1;   // ranking: 1+40 bits
    usint rankDecBits = 40;
    std::string csv;
    bool wilcoxon = true;
    bool bootstrap = false;
    bool verbose = false;
};

static void usage() {
    std::cout << "usage: statistics-test [--n N] [--runs R] [--depth D] [--intbits I] [--decbits F]\n"
                 "                       [--chi2-depth E] [--chi2-intbits I] [--chi2-decbits F]\n"
                 "                       [--compare-depth C] [--rank-intbits I] [--rank-decbits F]\n"
                 "                       [--csv FILE] [--no-wilcoxon] [--bootstrap] [--verbose]\n";
}

static Options parseArgs(int argc, char* argv[]) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        auto next = [&](const char* what) -> std::string {
            if (i + 1 >= argc) { std::cerr << "missing value for " << what << std::endl; usage(); std::exit(2); }
            return std::string(argv[++i]);
        };
        if (a == "--n") o.n = std::stoul(next("--n"));
        else if (a == "--runs") o.runs = std::stoul(next("--runs"));
        else if (a == "--depth") o.depth = static_cast<usint>(std::stoul(next("--depth")));
        else if (a == "--chi2-depth") o.chi2Depth = static_cast<usint>(std::stoul(next("--chi2-depth")));
        else if (a == "--chi2-intbits") o.chi2IntBits = static_cast<usint>(std::stoul(next("--chi2-intbits")));
        else if (a == "--chi2-decbits") o.chi2DecBits = static_cast<usint>(std::stoul(next("--chi2-decbits")));
        else if (a == "--compare-depth") o.compareDepth = static_cast<usint>(std::stoul(next("--compare-depth")));
        else if (a == "--intbits") o.intBits = static_cast<usint>(std::stoul(next("--intbits")));
        else if (a == "--decbits") o.decBits = static_cast<usint>(std::stoul(next("--decbits")));
        else if (a == "--rank-intbits") o.rankIntBits = static_cast<usint>(std::stoul(next("--rank-intbits")));
        else if (a == "--rank-decbits") o.rankDecBits = static_cast<usint>(std::stoul(next("--rank-decbits")));
        else if (a == "--csv") o.csv = next("--csv");
        else if (a == "--no-wilcoxon") o.wilcoxon = false;
        else if (a == "--bootstrap") o.bootstrap = true;
        else if (a == "--verbose") o.verbose = true;
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else { std::cerr << "unknown argument: " << a << std::endl; usage(); std::exit(2); }
    }
    return o;
}

int main(int argc, char* argv[]) {
    Options opt = parseArgs(argc, argv);
    g_verbose = opt.verbose;

    // one "center" in the single-party setting
    Cohort cohort = opt.csv.empty() ? syntheticCohort(opt.n, 1) : loadCohortCsv(opt.csv, 1);

    size_t maxN = 0;
    for (const auto& a : cohort.analyses) maxN = std::max(maxN, a.nEff);
    const size_t slots = nextPow2(maxN);
    const bool runWilcoxon = opt.wilcoxon && slots <= 64;

    std::cout << "Single-party FHE, " << (opt.csv.empty() ? "synthetic cohort" : "cohort file " + opt.csv)
              << ", vector length " << slots << ", runs = " << opt.runs << std::endl;
    if (opt.wilcoxon && !runWilcoxon)
        std::cout << "Wilcoxon skipped: vector length " << slots << " exceeds the ranking primitive's limit (64)." << std::endl;

    // one context per test: A Pearson, C chi-squared, B ranking
    Party A = makeParty(opt.intBits, opt.decBits, opt.depth, slots, slots, false, false);
    Party C = makeParty(opt.chi2IntBits, opt.chi2DecBits, opt.chi2Depth, slots, slots, false, false);
    Party B;
    if (runWilcoxon) {
        usint depthB = opt.bootstrap ? opt.compareDepth : static_cast<usint>(opt.compareDepth + 4);
        B = makeParty(opt.rankIntBits, opt.rankDecBits, depthB, slots, slots * slots, true, opt.bootstrap);
    }

    std::cout << std::fixed;
    std::cout << "SUMMARY single-party slots=" << slots << " depth=" << opt.depth
              << " bits=" << opt.intBits << "+" << opt.decBits
              << " chi2-depth=" << opt.chi2Depth << " chi2-bits=" << opt.chi2IntBits << "+" << opt.chi2DecBits
              << " compare-depth=" << opt.compareDepth
              << " wilcoxon-depth=" << (runWilcoxon ? (opt.bootstrap ? opt.compareDepth : opt.compareDepth + 4) : 0)
              << " rank-bits=" << opt.rankIntBits << "+" << opt.rankDecBits << std::endl;
    for (const auto& a : cohort.analyses) {
        if (a.kind == "wilcoxon" && !runWilcoxon) continue;
        const Party& p = (a.kind == "wilcoxon") ? B : (a.kind == "chi2" ? C : A);

        auto xC = encryptByCenter(p, padTo(a.x, slots), a.center, 1);
        auto yC = encryptByCenter(p, padTo(a.y, slots), a.center, 1);
        Plaintext indP = indicatorPlaintext(p, a.nEff);

        double plain = (a.kind == "pearson") ? pearsonPlain(a.x, a.y)
                     : (a.kind == "wilcoxon") ? wilcoxonPlain(a.x, a.y) : chi2Plain(a.x, a.y);

        double total = 0.0, enc = 0.0;
        for (size_t r = 0; r < opt.runs; ++r) {
            std::pair<double, double> tr;
            if (a.kind == "pearson")       tr = timed([&] { return pearsonCtxt(p, xC, yC, a.nEff, indP); });
            else if (a.kind == "wilcoxon") tr = timed([&] { return wilcoxonCtxt(p, xC, yC, a.nEff, slots, opt.compareDepth); });
            else                           tr = timed([&] { return chi2Ctxt(p, xC, yC, a.nEff, indP); });
            total += tr.first;
            enc = tr.second;
        }
        double dev = (a.kind == "pearson") ? std::fabs(enc - plain) : 100.0 * std::fabs(enc - plain) / std::fabs(plain);
        std::cout << "  " << std::left << std::setw(60) << a.name << std::right
                  << " n=" << std::setw(4) << a.nEff
                  << " plain=" << std::setw(10) << std::setprecision(4) << plain
                  << " enc=" << std::setw(10) << std::setprecision(4) << enc
                  << " dev=" << std::setw(8) << std::setprecision(4) << dev << (a.kind == "pearson" ? " (abs)" : " (%)  ")
                  << " ms=" << std::setw(9) << std::setprecision(2) << 1000.0 * total / opt.runs << std::endl;
    }
    return 0;
}
