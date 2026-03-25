#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "statistics.h"
#include "ranking.h"

#include <chrono>
#include <numeric>
#include <cmath>
#include <iostream>

using namespace lbcrypto;

// Simple helper to measure runtime of a callable
template <class F>
auto timed(F&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    auto res = f();
    auto end = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration<double>(end - start).count();
    return std::make_pair(t, res);
}

// Plaintext Pearson correlation
static double pearsonPlain(const std::vector<double>& x, const std::vector<double>& y) {
    double meanX = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    double meanY = std::accumulate(y.begin(), y.end(), 0.0) / y.size();
    double num = 0.0, denX = 0.0, denY = 0.0;
    for (size_t i = 0; i < x.size(); ++i) {
        double dx = x[i] - meanX;
        double dy = y[i] - meanY;
        num += dx * dy;
        denX += dx * dx;
        denY += dy * dy;
    }
    return num / (std::sqrt(denX) * std::sqrt(denY));
}

// Ciphertext Pearson correlation
static double pearsonCtxt(const CryptoContext<DCRTPoly>& cc,
                           const KeyPair<DCRTPoly>& keys,
                           Ciphertext<DCRTPoly> xC,
                           Ciphertext<DCRTPoly> yC,
                           size_t n) {

    auto meanXC = meanVector(cc, keys.publicKey, xC, n);
    auto meanYC = meanVector(cc, keys.publicKey, yC, n);

    auto xCen = centerVector(cc, keys.publicKey, xC, n);
    auto yCen = centerVector(cc, keys.publicKey, yC, n);

    auto xyC = dotProduct(cc, keys.publicKey, xCen, yCen, n);
    xyC = sum (cc, keys.publicKey, xyC, n);
    auto x2C = dotProduct(cc, keys.publicKey, xCen, xCen, n);
    x2C = sum (cc, keys.publicKey, x2C, n);
    auto y2C = dotProduct(cc, keys.publicKey, yCen, yCen, n);
    y2C = sum (cc, keys.publicKey, y2C, n);

    double u1 = 1.23; // or use random
    double u2 = 2.34;
    double u3 = u1 / u2;

    auto maskedXdotY = cc->EvalMult(xyC, u1);
    auto maskedLenX  = cc->EvalMult(x2C, u2*u2);
    auto maskedLenY  = cc->EvalMult(y2C, u3*u3);

    // After masking
    double masked_dot = decryptScalar(cc, keys.secretKey, maskedXdotY);
    double masked_x2  = decryptScalar(cc, keys.secretKey, maskedLenX);
    double masked_y2  = decryptScalar(cc, keys.secretKey, maskedLenY);

    double numerator   = masked_dot;
    double denominator = std::sqrt(masked_x2) * std::sqrt(masked_y2);
    double correction  = u1 / (u2 * u3);

    return correction * numerator / denominator;
}

// Plaintext Wilcoxon rank-sum statistic (Z)
static double wilcoxonPlain(const std::vector<double>& x, const std::vector<double>& g) {
    auto ranks = rank(x); // fractional ranking
    double S1 = 0.0;
    double n1 = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        if (g[i] > 0.5) {
            S1 += ranks[i];
            n1 += 1.0;
        }
    }
    double n0 = static_cast<double>(x.size()) - n1;
    double U1 = S1 - ( n1 * (n1 + 1) ) / 2.0;
    return (U1 - n0 * n1 / 2.0) / std::sqrt(n0 * n1 * (n0 + n1 + 1) / 12.0);
}

// Plaintext chi-squared for 2x2 contingency table
static double chi2Plain(const std::vector<double>& t, const std::vector<double>& l) {
    double a=0,b=0,c=0,d=0;
    for (size_t i = 0; i < t.size(); ++i) {
        a += (1 - t[i]) * (1 - l[i]);
        b += (1 - t[i]) * l[i];
        c += t[i] * (1 - l[i]);
        d += t[i] * l[i];
    }
    double r0 = a + b;
    double r1 = c + d;
    double c0 = a + c;
    double c1 = b + d;
    double N  = r0 + r1;
    double ea = r0 * c0 / N;
    double eb = r0 * c1 / N;
    double ec = r1 * c0 / N;
    double ed = r1 * c1 / N;
    return (std::pow(a-ea,2)/ea) + (std::pow(b-eb,2)/eb) +
           (std::pow(c-ec,2)/ec) + (std::pow(d-ed,2)/ed);
}

// Ciphertext Wilcoxon rank-sum statistic
static double wilcoxonCtxt(const CryptoContext<DCRTPoly>& cc,
                           const KeyPair<DCRTPoly>& keys,
                           Ciphertext<DCRTPoly> xC,
                           Ciphertext<DCRTPoly> gC,
                           size_t n,
                           usint compareDepth) {
    Ciphertext<DCRTPoly> rankC = rank(xC, n, -1.0, 1.0, depth2degree(compareDepth));
    auto s1C = dotProduct(cc, keys.publicKey, rankC, gC, n);
    auto n1C = sum(cc, keys.publicKey, gC, n);

    double S1 = decryptScalar(cc, keys.secretKey, s1C);
    double n1 = decryptScalar(cc, keys.secretKey, n1C);
    double n0 = static_cast<double>(n) - n1;

    auto minusOneC = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(std::vector<double>(n, -1.0)));

    // Encrypt the size n as a vector filled with n
    std::vector<double> sizeVector(n, static_cast<double>(n));
    auto sizeC = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(sizeVector));

    auto minusn1C = cc->EvalMult(n1C, minusOneC);
    auto n0C = cc->EvalAdd(sizeC, minusn1C);

    double U1 = S1 - (n1 * (n1 + 1) ) / 2.0;
    return (U1 - (n0 * n1 ) / 2.0) / std::sqrt((n0 * n1 * (n0 + n1 + 1) ) / 12.0);
}

// Ciphertext chi-squared
static double chi2Ctxt(const CryptoContext<DCRTPoly>& cc,
                       const KeyPair<DCRTPoly>& keys,
                       Ciphertext<DCRTPoly> tC,
                       Ciphertext<DCRTPoly> lC,
                       Ciphertext<DCRTPoly> oneC,
                       size_t n) {
    auto aC = sum(cc, keys.publicKey, (oneC - tC) * (oneC - lC), n);
    auto bC = sum(cc, keys.publicKey, (oneC - tC) * lC, n);
    auto cC = sum(cc, keys.publicKey, tC * (oneC - lC), n);
    auto dC = sum(cc, keys.publicKey, tC * lC, n);

    auto r0C = aC + bC;
    auto r1C = cC + dC;
    auto c0C = aC + cC;
    auto c1C = bC + dC;
    auto nC = aC+bC+cC+dC;

    std::vector<double> sizeVector(n, static_cast<double>(1.0/n));
    auto onedivdedbysizeVector = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(sizeVector));

    double a = decryptScalar(cc, keys.secretKey, aC);
    double b = decryptScalar(cc, keys.secretKey, bC);
    double c = decryptScalar(cc, keys.secretKey, cC);
    double d = decryptScalar(cc, keys.secretKey, dC);

    double r0 = a + b;
    double r1 = c + d;
    double c0 = a + c;
    double c1 = b + d;
    double N  = r0 + r1;

    double ea = r0 * c0 / N;
    double eb = r0 * c1 / N;
    double ec = r1 * c0 / N;
    double ed = r1 * c1 / N;
    return (std::pow(a-ea,2)/ea) + (std::pow(b-eb,2)/eb) +
           (std::pow(c-ec,2)/ec) + (std::pow(d-ed,2)/ed);
}

int main() {
    // Sample data
    std::vector<double> x  = {0.2, 0.4, 0.6, 0.8};
    std::vector<double> y  = {0.1, 0.5, 0.3, 0.9};

    std::vector<double> g = { 1,0,1,0};

    std::vector<double> t  = {0, 0, 1, 1};
    std::vector<double> l  = {0, 1, 0, 1};

    const size_t n = x.size();
    const usint compareDepth = 3;
    const size_t runs = 50; // battery size
    const size_t slots = n;
    std::cout<<"generating context..."<<std::endl;
    CryptoContext<DCRTPoly> cc = generateCryptoContext(
        1, 20, compareDepth , slots, false, 0, true);

    std::cout<<"context generated"<<std::endl;

    std::vector<int32_t> indices = getRotationIndices(n);
    std::cout<<"generating keys..."<<std::endl;
    KeyPair<DCRTPoly> keys = keyGeneration(cc, indices, slots, false, true);

    auto enc = [&](const std::vector<double>& v){
        Plaintext pt = cc->MakeCKKSPackedPlaintext(v, 1, 0, nullptr, slots);
        return cc->Encrypt(keys.publicKey, pt);
    };

    Ciphertext<DCRTPoly> xC = enc(x);
    Ciphertext<DCRTPoly> yC = enc(y);
    Ciphertext<DCRTPoly> gC = enc(g);
    Ciphertext<DCRTPoly> tC = enc(t);
    Ciphertext<DCRTPoly> lC = enc(l);
    std::vector<double> ones(n,1.0);
    Ciphertext<DCRTPoly> oneC = enc(ones);

    auto [ent, rpt] = timed([&]() {
        auto encOne = cc->Encrypt(keys.publicKey, cc->MakeCKKSPackedPlaintext(std::vector<double>(n, 1.0)));
        return encOne;
    });
    std::cout<<"encryption time: "<<ent*1000<<"ms"<<std::endl;

    std::cout<<"Number of values: "<<n<<std::endl;

    // Pearson
    double tPlain=0, tCtxt=0, resPlain=0, resCtxt=0;

    for(size_t i=0;i<runs;i++) {
        auto [tp, rp] = timed([&]{ return pearsonPlain(x,y); });
        tPlain += tp; resPlain = rp;
        auto [tc, rc] = timed([&]{ return pearsonCtxt(cc, keys, xC, yC, n); });
        tCtxt += tc; resCtxt = rc;
    }
    std::cout << "Pearson plaintext result: " << resPlain
              << " avg time: " << (tPlain / runs) * 1000 << "ms" << std::endl;
    std::cout << "Pearson ciphertext result: " << resCtxt
              << " avg time: " << (tCtxt / runs) * 1000 << "ms" << std::endl;

    std::cout<<std::endl<<std::endl;
    return 0;

    // Wilcoxon
    tPlain=tCtxt=0; resPlain=resCtxt=0;
    for(size_t i=0;i<runs;i++) {
        auto [tp, rp] = timed([&]{ return wilcoxonPlain(x,g); });
        tPlain += tp; resPlain = rp;
        auto [tc, rc] = timed([&]{ return wilcoxonCtxt(cc, keys, xC, gC, n, compareDepth); });
        tCtxt += tc; resCtxt = rc;
    }
    std::cout << "Wilcoxon plaintext Z: " << resPlain
              << " avg time: " << (tPlain / runs) * 1000 << "ms" << std::endl;
    std::cout << "Wilcoxon ciphertext Z: " << resCtxt
              << " avg time: " << (tCtxt / runs) * 1000 << "ms" << std::endl;

    std::cout<<std::endl<<std::endl;

    // Chi-squared
    tPlain=tCtxt=0; resPlain=resCtxt=0;
    for(size_t i=0;i<runs;i++) {
        auto [tp, rp] = timed([&]{ return chi2Plain(t,l); });
        tPlain += tp; resPlain = rp;
        auto [tc, rc] = timed([&]{ return chi2Ctxt(cc, keys, tC, lC, oneC, n); });
        tCtxt += tc; resCtxt = rc;
    }

    std::cout << "Chi-squared plaintext: " << resPlain
              << " avg time: " << (tPlain / runs) * 1000 << "ms" << std::endl;
    std::cout << "Chi-squared ciphertext: " << resCtxt
              << " avg time: " << (tCtxt / runs) * 1000 << "ms" << std::endl;

    return 0;
}
