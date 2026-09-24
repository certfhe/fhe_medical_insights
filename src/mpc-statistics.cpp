// Threshold HE version (2 or 3 parties) of the three tests used in the paper:
// Pearson correlation, Wilcoxon rank-sum and chi-squared.
//
// Usage:
//   mpc-statistics [--parties 2|3] [--n N] [--runs R] [--csv FILE] [--no-wilcoxon] [--verbose]
//                  [--depth D] [--intbits I] [--decbits F]                   (Pearson context)
//                  [--chi2-depth E] [--chi2-intbits I] [--chi2-decbits F]    (chi-squared context)
//                  [--compare-depth C] [--wilcoxon-depth B] [--rank-intbits I] [--rank-decbits F]
//   mpc-statistics 2party|3party    (old form, still works)
//
// Without --csv a synthetic cohort of size N is used (fixed seed), as for the
// benchmark tables. With --csv the five analyses of Table 1 are computed from a
// cohort file (format in the README). Each emulated center encrypts its own rows
// and the server adds the ciphertexts.
//
// Each test has its own CKKS context (128-bit security, FLEXIBLEAUTO, HYBRID):
//   Pearson      depth 8,  1+35 bits, ring dim 16384
//   chi-squared  depth 5, 10+50 bits, ring dim 16384 (with 25-bit scaling it fits
//                in 8192, but then the error is too large for weak associations)
//   ranking      depth 14, 1+40 bits, ring dim 32768, interactive bootstrapping
//
// All tests end with the client-aided step (Sections 3.4-3.5 of the paper): the
// server multiplies the aggregates by fresh random masks (mask.h), only the
// masked values are threshold decrypted, and the analyst gets the statistic as
// a ratio in which the masks cancel. The released values are divided by powers
// of n so that they stay small (the CKKS decoder fails when a decrypted value
// carries a large error).

#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "statistics.h"
#include "ranking.h"
#include "ranking-mpc.h"
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

// Simple helper to measure runtime of a callable
template <class F>
auto timed(F&& f) {
    auto start = std::chrono::high_resolution_clock::now();
    auto res = f();
    auto end = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration<double>(end - start).count();
    return std::make_pair(t, res);
}

static bool g_verbose = false;

struct MPCKeys {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> party1Keys;
    KeyPair<DCRTPoly> party2Keys;
    KeyPair<DCRTPoly> party3Keys;
    PublicKey<DCRTPoly> jointPublicKey;
    usint integralPrecision;    // Bits for integer part
    usint decimalPrecision;     // Bits for fractional part
    size_t numParties = 3;      // 2 or 3 parties
    size_t slots = 0;           // batch size of this context
};

// MPC Decryption helper (slot 0)
double decryptMPCScalar(const MPCKeys& mpcKeys, const Ciphertext<DCRTPoly>& ciphertext) {
    // Party 1 performs lead decryption
    auto ciphertextPartial1 = mpcKeys.cc->MultipartyDecryptLead({ciphertext}, mpcKeys.party1Keys.secretKey);

    // Party 2 performs main decryption
    auto ciphertextPartial2 = mpcKeys.cc->MultipartyDecryptMain({ciphertext}, mpcKeys.party2Keys.secretKey);

    std::vector<Ciphertext<DCRTPoly>> partialCiphertextVec;
    partialCiphertextVec.push_back(ciphertextPartial1[0]);
    partialCiphertextVec.push_back(ciphertextPartial2[0]);

    if (mpcKeys.numParties == 3) {
        // Party 3 performs main decryption
        auto ciphertextPartial3 = mpcKeys.cc->MultipartyDecryptMain({ciphertext}, mpcKeys.party3Keys.secretKey);
        partialCiphertextVec.push_back(ciphertextPartial3[0]);
    }

    Plaintext plaintextResult;
    mpcKeys.cc->MultipartyDecryptFusion(partialCiphertextVec, &plaintextResult);

    return plaintextResult->GetRealPackedValue()[0];
}

// MPC Decryption helper (full vector)
std::vector<double> decryptMPC(const MPCKeys& mpcKeys, const Ciphertext<DCRTPoly>& ciphertext) {
    auto ciphertextPartial1 = mpcKeys.cc->MultipartyDecryptLead({ciphertext}, mpcKeys.party1Keys.secretKey);
    auto ciphertextPartial2 = mpcKeys.cc->MultipartyDecryptMain({ciphertext}, mpcKeys.party2Keys.secretKey);

    std::vector<Ciphertext<DCRTPoly>> partialCiphertextVec;
    partialCiphertextVec.push_back(ciphertextPartial1[0]);
    partialCiphertextVec.push_back(ciphertextPartial2[0]);

    if (mpcKeys.numParties == 3) {
        auto ciphertextPartial3 = mpcKeys.cc->MultipartyDecryptMain({ciphertext}, mpcKeys.party3Keys.secretKey);
        partialCiphertextVec.push_back(ciphertextPartial3[0]);
    }

    Plaintext plaintextResult;
    mpcKeys.cc->MultipartyDecryptFusion(partialCiphertextVec, &plaintextResult);

    return plaintextResult->GetRealPackedValue();
}

// Interactive MPC bootstrapping (TCKKS IntMPBoot) for 2 or 3 parties
static Ciphertext<DCRTPoly> mpcBootstrap(const MPCKeys& mpcKeys, Ciphertext<DCRTPoly> ct) {

    auto cc = mpcKeys.cc;

    // Adjust towers/scale for interactive bootstrapping
    auto inCtxt = cc->IntMPBootAdjustScale(ct);

    // Generate the common random polynomial using the joint public key
    auto a = cc->IntMPBootRandomElementGen(mpcKeys.jointPublicKey);

    // Prepare c1 (remove c0 component)
    auto c1 = inCtxt->Clone();
    c1->GetElements().erase(c1->GetElements().begin());

    // Each party produces masked decryption and re-encryption shares
    auto sharesP1 = cc->IntMPBootDecrypt(mpcKeys.party1Keys.secretKey, c1, a);
    auto sharesP2 = cc->IntMPBootDecrypt(mpcKeys.party2Keys.secretKey, c1, a);

    std::vector<std::vector<Ciphertext<DCRTPoly>>> sharesPairVec{sharesP1, sharesP2};
    if (mpcKeys.numParties == 3) {
        auto sharesP3 = cc->IntMPBootDecrypt(mpcKeys.party3Keys.secretKey, c1, a);
        sharesPairVec.push_back(sharesP3);
    }

    // Aggregate and re-encrypt to obtain refreshed ciphertext
    auto aggregated = cc->IntMPBootAdd(sharesPairVec);
    return cc->IntMPBootEncrypt(mpcKeys.jointPublicKey, aggregated, a, inCtxt);
}

// MPC key generation for 2 or 3 parties.
// slots is the batch size (n*n for the ranking context, n otherwise); the rotation
// keys for the n x n matrix layout are only needed for the ranking.
MPCKeys generateMPCKeys(const usint compareDepth,
                       const usint integralPrecision,
                       const usint decimalPrecision,
                       const size_t matrixSize,
                       const size_t slots,
                       const size_t numParties,
                       const bool needMatrixKeys) {
    std::cout << "Setting up " << numParties << "-Party MPC keys (slots = " << slots
              << ", depth = " << compareDepth << ")..." << std::endl;

    // Calculate total precision
    const usint totalPrecision = integralPrecision + decimalPrecision;

    // Create parameters with multiparty mode enabled
    CCParams<CryptoContextCKKSRNS> parameters;
    parameters.SetMultiplicativeDepth(compareDepth);
    parameters.SetScalingModSize(decimalPrecision);        // Fractional bits
    parameters.SetFirstModSize(totalPrecision);            // Total precision
    parameters.SetBatchSize(slots);
    parameters.SetSecurityLevel(HEStd_128_classic);
    parameters.SetSecretKeyDist(UNIFORM_TERNARY);
    parameters.SetScalingTechnique(ScalingTechnique::FLEXIBLEAUTO);
    parameters.SetKeySwitchTechnique(KeySwitchTechnique::HYBRID);
    parameters.SetInteractiveBootCompressionLevel(COMPRESSION_LEVEL::SLACK);

    CryptoContext<DCRTPoly> cc = GenCryptoContext(parameters);
    cc->Enable(PKE);
    cc->Enable(KEYSWITCH);
    cc->Enable(LEVELEDSHE);
    cc->Enable(ADVANCEDSHE);
    cc->Enable(MULTIPARTY);

    if (g_verbose) {
        const BigInteger ciphertextModulus = cc->GetModulus();
        const usint ciphertextModulusBitsize = ciphertextModulus.GetLengthForBase(2);
        const usint ringDimension = cc->GetRingDimension();
        const auto elementParameters = cc->GetCryptoParameters()->GetElementParams()->GetParams();
        std::vector<usint> moduliChainBitsize(elementParameters.size());
        for (size_t i = 0; i < elementParameters.size(); i++)
            moduliChainBitsize[i] = elementParameters[i]->GetModulus().GetLengthForBase(2);

        // sizes for Table 7: a fresh ciphertext has 2*N*L words,
        // a key-switching key 2*dnum*(L+|P|)*N words
        const auto rnsParams = std::dynamic_pointer_cast<CryptoParametersRNS>(cc->GetCryptoParameters());
        const size_t numQ = elementParameters.size();
        const size_t numP = (rnsParams->GetParamsP() ? rnsParams->GetParamsP()->GetParams().size() : 0);
        const uint32_t dnum = rnsParams->GetNumPartQ();
        const double mib = 1024.0 * 1024.0;
        const double ctBytes = 2.0 * ringDimension * numQ * 8.0;
        const double keyBytes = 2.0 * dnum * (numQ + numP) * ringDimension * 8.0;

        std::cout << "CKKS PARAMETERS" << std::endl
                  << "Integral Bit Precision        : " << integralPrecision << std::endl
                  << "Decimal Bit Precision         : " << decimalPrecision << std::endl
                  << "Ciphertext Modulus Precision  : " << ciphertextModulusBitsize << std::endl
                  << "Ring Dimension                : " << ringDimension << std::endl
                  << "Max Slots                     : " << ringDimension / 2 << std::endl
                  << "Slots                         : " << slots << std::endl
                  << "Multiplicative Depth          : " << parameters.GetMultiplicativeDepth() << std::endl
                  << "Moduli Chain Bitsize          : " << moduliChainBitsize << std::endl
                  << "CRT primes in Q / in P        : " << numQ << " / " << numP << std::endl
                  << "Key-switching digits (dnum)   : " << dnum << std::endl
                  << "Fresh ciphertext (2 N L words): " << ctBytes / mib << " MiB" << std::endl
                  << "One key-switching key         : " << keyBytes / mib << " MiB" << std::endl
                  << std::endl;
    }

    MPCKeys mpcKeys;
    mpcKeys.cc = cc;
    mpcKeys.integralPrecision = integralPrecision;
    mpcKeys.decimalPrecision = decimalPrecision;
    mpcKeys.numParties = numParties;
    mpcKeys.slots = slots;

    // Rotation indices: matrix layout of the ranking (if needed) and +/- powers of two for the slot sums
    std::set<int32_t> indexSet;
    if (needMatrixKeys)
        for (int32_t idx : getRotationIndices(matrixSize)) indexSet.insert(idx);
    for (size_t i = 1; i < slots; i *= 2) {
        indexSet.insert(static_cast<int32_t>(i));
        indexSet.insert(-static_cast<int32_t>(i));
    }
    std::vector<int32_t> indices(indexSet.begin(), indexSet.end());

    if (g_verbose) {
        // public evaluation material: relinearisation, rotation and EvalSum keys
        size_t log2Slots = 0;
        for (size_t s = slots; s > 1; s /= 2) ++log2Slots;
        const auto rnsParams = std::dynamic_pointer_cast<CryptoParametersRNS>(cc->GetCryptoParameters());
        const size_t numQ = cc->GetCryptoParameters()->GetElementParams()->GetParams().size();
        const size_t numP = (rnsParams->GetParamsP() ? rnsParams->GetParamsP()->GetParams().size() : 0);
        const double keyMiB = 2.0 * rnsParams->GetNumPartQ() * (numQ + numP) * cc->GetRingDimension() * 8.0 / (1024.0 * 1024.0);
        const size_t numKeys = 1 + indices.size() + log2Slots;
        std::cout << "Evaluation keys               : 1 relinearisation + " << indices.size() << " rotation + "
                  << log2Slots << " EvalSum keys = " << numKeys << " keys, about "
                  << numKeys * keyMiB << " MiB of public material" << std::endl << std::endl;
    }

    // Round 1: Party 1 generates initial keys
    mpcKeys.party1Keys = cc->KeyGen();
    cc->EvalMultKeyGen(mpcKeys.party1Keys.secretKey);
    cc->EvalSumKeyGen(mpcKeys.party1Keys.secretKey);
    cc->EvalRotateKeyGen(mpcKeys.party1Keys.secretKey, indices);

    auto evalSumKeys1 = std::make_shared<std::map<usint, EvalKey<DCRTPoly>>>(
        cc->GetEvalSumKeyMap(mpcKeys.party1Keys.secretKey->GetKeyTag()));

    // Party 2 contributes
    mpcKeys.party2Keys = cc->MultipartyKeyGen(mpcKeys.party1Keys.publicKey);
    cc->EvalMultKeyGen(mpcKeys.party2Keys.secretKey);
    cc->EvalSumKeyGen(mpcKeys.party2Keys.secretKey);
    cc->EvalRotateKeyGen(mpcKeys.party2Keys.secretKey, indices);

    // Combine evaluation keys (start with party1)
    auto evalMultKey1 = cc->KeySwitchGen(mpcKeys.party1Keys.secretKey, mpcKeys.party1Keys.secretKey);
    auto evalMultKey2 = cc->MultiKeySwitchGen(mpcKeys.party2Keys.secretKey, mpcKeys.party2Keys.secretKey, evalMultKey1);

    auto evalSumKeys2 = cc->MultiEvalSumKeyGen(mpcKeys.party2Keys.secretKey, evalSumKeys1, mpcKeys.party2Keys.publicKey->GetKeyTag());
    auto evalSumKeys12 = cc->MultiAddEvalSumKeys(evalSumKeys1, evalSumKeys2, mpcKeys.party2Keys.publicKey->GetKeyTag());

    // Combine rotation keys for MPC (party1 + party2)
    auto evalRotateKeys1 = std::make_shared<std::map<usint, EvalKey<DCRTPoly>>>(
        cc->GetEvalAutomorphismKeyMap(mpcKeys.party1Keys.secretKey->GetKeyTag()));
    auto evalRotateKeys2 = cc->MultiEvalAtIndexKeyGen(mpcKeys.party2Keys.secretKey, evalRotateKeys1, indices, mpcKeys.party2Keys.publicKey->GetKeyTag());
    auto evalRotateKeys12 = cc->MultiAddEvalAutomorphismKeys(evalRotateKeys1, evalRotateKeys2, mpcKeys.party2Keys.publicKey->GetKeyTag());

    if (numParties == 2) {
        cc->InsertEvalSumKey(evalSumKeys12);
        cc->InsertEvalAutomorphismKey(evalRotateKeys12);

        auto evalMultAB = cc->MultiAddEvalKeys(evalMultKey1, evalMultKey2, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultAAB = cc->MultiMultEvalKey(mpcKeys.party1Keys.secretKey, evalMultAB, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultBAB = cc->MultiMultEvalKey(mpcKeys.party2Keys.secretKey, evalMultAB, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMultAAB, evalMultBAB, evalMultAAB->GetKeyTag());
        cc->InsertEvalMultKey({evalMultFinal});

        // Use last party's public key as the joint encryption key (matches eval keys' tag)
        mpcKeys.jointPublicKey = mpcKeys.party2Keys.publicKey;
        std::cout << "2-Party MPC key generation completed." << std::endl;
        return mpcKeys;
    }

    // Party 3 contributes
    mpcKeys.party3Keys = cc->MultipartyKeyGen(mpcKeys.party2Keys.publicKey);
    cc->EvalMultKeyGen(mpcKeys.party3Keys.secretKey);
    cc->EvalSumKeyGen(mpcKeys.party3Keys.secretKey);
    cc->EvalRotateKeyGen(mpcKeys.party3Keys.secretKey, indices);

    auto evalMultKey3 = cc->MultiKeySwitchGen(mpcKeys.party3Keys.secretKey, mpcKeys.party3Keys.secretKey, evalMultKey1);
    auto evalMultAB = cc->MultiAddEvalKeys(evalMultKey1, evalMultKey2, mpcKeys.party2Keys.publicKey->GetKeyTag());
    auto evalMultABC = cc->MultiAddEvalKeys(evalMultAB, evalMultKey3, mpcKeys.party3Keys.publicKey->GetKeyTag());

    auto evalSumKeys3 = cc->MultiEvalSumKeyGen(mpcKeys.party3Keys.secretKey, evalSumKeys1, mpcKeys.party3Keys.publicKey->GetKeyTag());
    auto evalSumKeys123 = cc->MultiAddEvalSumKeys(evalSumKeys12, evalSumKeys3, mpcKeys.party3Keys.publicKey->GetKeyTag());
    cc->InsertEvalSumKey(evalSumKeys123);

    auto evalRotateKeys3 = cc->MultiEvalAtIndexKeyGen(mpcKeys.party3Keys.secretKey, evalRotateKeys1, indices, mpcKeys.party3Keys.publicKey->GetKeyTag());
    auto evalRotateKeys123 = cc->MultiAddEvalAutomorphismKeys(evalRotateKeys12, evalRotateKeys3, mpcKeys.party3Keys.publicKey->GetKeyTag());
    cc->InsertEvalAutomorphismKey(evalRotateKeys123);

    auto evalMultAABC = cc->MultiMultEvalKey(mpcKeys.party1Keys.secretKey, evalMultABC, mpcKeys.party3Keys.publicKey->GetKeyTag());
    auto evalMultBABC = cc->MultiMultEvalKey(mpcKeys.party2Keys.secretKey, evalMultABC, mpcKeys.party3Keys.publicKey->GetKeyTag());
    auto evalMultCABC = cc->MultiMultEvalKey(mpcKeys.party3Keys.secretKey, evalMultABC, mpcKeys.party3Keys.publicKey->GetKeyTag());
    auto evalMultAB_ABC = cc->MultiAddEvalMultKeys(evalMultAABC, evalMultBABC, evalMultAABC->GetKeyTag());
    auto evalMultFinal  = cc->MultiAddEvalMultKeys(evalMultAB_ABC, evalMultCABC, evalMultAB_ABC->GetKeyTag());
    cc->InsertEvalMultKey({evalMultFinal});

    mpcKeys.jointPublicKey = mpcKeys.party3Keys.publicKey;
    std::cout << "3-Party MPC key generation completed." << std::endl;
    return mpcKeys;
}

// ---------------------------------------------------------------------------
// Encryption at the centers
// ---------------------------------------------------------------------------

// Every emulated center encrypts only its own rows (zeros elsewhere) with the
// joint public key, the server adds the ciphertexts.
static Ciphertext<DCRTPoly> encryptByCenter(const MPCKeys& mpcKeys,
                                            const std::vector<double>& values,
                                            const std::vector<int>& center,
                                            size_t numCenters) {
    auto parts = splitByCenter(values, center, numCenters, mpcKeys.slots);
    Ciphertext<DCRTPoly> sum;
    for (size_t c = 0; c < numCenters; ++c) {
        Plaintext pt = mpcKeys.cc->MakeCKKSPackedPlaintext(parts[c], 1, 0, nullptr, mpcKeys.slots);
        auto ct = mpcKeys.cc->Encrypt(mpcKeys.jointPublicKey, pt);
        sum = (c == 0) ? ct : mpcKeys.cc->EvalAdd(sum, ct);
    }
    return sum;
}

static Plaintext indicatorPlaintext(const MPCKeys& mpcKeys, size_t nEff) {
    return mpcKeys.cc->MakeCKKSPackedPlaintext(indicatorVector(nEff, mpcKeys.slots), 1, 0, nullptr, mpcKeys.slots);
}

// ---------------------------------------------------------------------------
// Homomorphic building blocks
// ---------------------------------------------------------------------------

// Sum of the first `length` slots, result in slot 0
static Ciphertext<DCRTPoly> sumMPC(const CryptoContext<DCRTPoly>& cc, Ciphertext<DCRTPoly> c, size_t length) {
    auto sum = c;
    for (size_t i = 1; i < length; i *= 2) {
        auto rotated = cc->EvalRotate(sum, i);
        sum = cc->EvalAdd(sum, rotated);
    }
    return sum;
}

// Sum over all slots of the batch
static Ciphertext<DCRTPoly> sumAllSlotsMPC(const MPCKeys& mpcKeys, Ciphertext<DCRTPoly> c) {
    return sumMPC(mpcKeys.cc, c, mpcKeys.slots);
}

// Mean of the nEff real entries (the padding is zero)
static Ciphertext<DCRTPoly> meanVectorMPC(const MPCKeys& mpcKeys, Ciphertext<DCRTPoly> c, size_t nEff) {
    auto sum = sumAllSlotsMPC(mpcKeys, c);
    return mpcKeys.cc->EvalMult(sum, 1.0 / static_cast<double>(nEff));
}

// ---------------------------------------------------------------------------
// Pearson correlation (Section 3.4.1)
// ---------------------------------------------------------------------------

// MPC Pearson correlation, returns r as the analyst reconstructs it
static double pearsonMPC(const MPCKeys& mpcKeys,
                         Ciphertext<DCRTPoly> xC,
                         Ciphertext<DCRTPoly> yC,
                         size_t nEff,
                         const Plaintext& indP) {
    auto cc = mpcKeys.cc;

    // Step 1: center the data, the public indicator sets the padding back to zero
    auto meanXC = meanVectorMPC(mpcKeys, xC, nEff);
    auto meanYC = meanVectorMPC(mpcKeys, yC, nEff);
    auto xCen = cc->EvalMult(cc->EvalSub(xC, meanXC), indP);
    auto yCen = cc->EvalMult(cc->EvalSub(yC, meanYC), indP);

    // Step 2: XdotY, LenX, LenY (slot 0)
    auto xyC = sumAllSlotsMPC(mpcKeys, cc->EvalMult(xCen, yCen));
    auto x2C = sumAllSlotsMPC(mpcKeys, cc->EvalMult(xCen, xCen));
    auto y2C = sumAllSlotsMPC(mpcKeys, cc->EvalMult(yCen, yCen));

    // Step 3: random masks u1, u2 and u3 = u1/u2 (Section 3.5). We release the
    // covariance and variances (sums divided by n) so that the values stay small;
    // 1/n cancels in r and goes into the same multiplication as the mask.
    const double u1 = sampleMask();
    const double u2 = sampleMask();
    const double u3 = u1 / u2;
    const double invN = 1.0 / static_cast<double>(nEff);

    auto maskedXdotY = cc->EvalMult(xyC, u1 * invN);
    auto maskedLenX  = cc->EvalMult(x2C, u2 * u2 * invN);
    auto maskedLenY  = cc->EvalMult(y2C, u3 * u3 * invN);

    // Step 4: threshold decryption of the masked values
    const double m1 = decryptMPCScalar(mpcKeys, maskedXdotY);
    const double m2 = decryptMPCScalar(mpcKeys, maskedLenX);
    const double m3 = decryptMPCScalar(mpcKeys, maskedLenY);

    if (g_verbose) {
        std::ostringstream os;   // separate stream, so the format of cout is not changed
        os << std::setprecision(6) << "    masked release (Pearson): " << m1 << ", " << m2 << ", " << m3
           << "  ->  r = " << m1 / (std::sqrt(m2) * std::sqrt(m3));
        std::cout << os.str() << std::endl;
    }

    // r = m1 / (sqrt(m2) sqrt(m3)), the masks cancel because u1 / (u2 u3) = 1
    return m1 / (std::sqrt(m2) * std::sqrt(m3));
}

// ---------------------------------------------------------------------------
// Wilcoxon rank-sum test (Section 3.4.2)
// ---------------------------------------------------------------------------

// MPC Wilcoxon rank-sum statistic.
// x is scaled to [0.2, 1] and padded with zeros, so the padding gets the lowest
// ranks and every real rank is shifted by p = slots - nEff (corrected below).
// g is the 0/1 group label, also zero padded.
static double wilcoxonMPC(const MPCKeys& mpcKeys,
                          Ciphertext<DCRTPoly> xC,
                          Ciphertext<DCRTPoly> gC,
                          size_t nEff,
                          size_t matrixSize,
                          usint compareDepth) {
    auto cc = mpcKeys.cc;

    // approximate ranks, in the first matrixSize slots
    Ciphertext<DCRTPoly> rankC = rankMPC(mpcKeys, xC, matrixSize, -1.0, 1.0, depth2degree(compareDepth));
    Ciphertext<DCRTPoly> rankVec = maskRow(rankC, matrixSize, 0);

    // MPC interactive bootstrapping on rankVec before using it
    rankVec = mpcBootstrap(mpcKeys, rankVec);

    // S1 = sum_i R_i g_i and n1 = sum_i g_i (slot 0)
    auto s1C = sumMPC(cc, cc->EvalMult(rankVec, gC), matrixSize);
    auto n1C = sumMPC(cc, gC, matrixSize);

    const double n = static_cast<double>(nEff);
    const double p = static_cast<double>(matrixSize - nEff);   // number of padded entries

    // D1 = S1_true - n1 (n+1)/2 = S1_padded - n1 (p + (n+1)/2)     (= U1 - n0 n1 / 2)
    auto d1C = cc->EvalSub(s1C, cc->EvalMult(n1C, p + (n + 1.0) / 2.0));
    // D2 = n1 (n - n1) (n+1) / 12
    auto n0C = cc->EvalSub(n, n1C);
    auto d2C = cc->EvalMult(cc->EvalMult(n1C, n0C), (n + 1.0) / 12.0);

    // masked release (u D1/n^2, u^2 D2/n^4) with one random mask u > 0,
    // the powers of n keep the values small and cancel in Z
    const double u = sampleMask();
    auto maskedD1 = cc->EvalMult(d1C, u / (n * n));
    auto maskedD2 = cc->EvalMult(d2C, u * u / (n * n * n * n));

    const double m1 = decryptMPCScalar(mpcKeys, maskedD1);
    const double m2 = decryptMPCScalar(mpcKeys, maskedD2);

    if (g_verbose) {
        std::ostringstream os;
        os << std::setprecision(6) << "    masked release (Wilcoxon): " << m1 << ", " << m2
           << "  ->  z = " << m1 / std::sqrt(m2);
        std::cout << os.str() << std::endl;
    }

    // Z = D1 / sqrt(D2), u cancels and u > 0 keeps the sign
    return m1 / std::sqrt(m2);
}

// ---------------------------------------------------------------------------
// Chi-squared test (Section 3.4.3)
// ---------------------------------------------------------------------------

// MPC chi-squared for the 2x2 table.
// t = technique, l = label (0/1, zero padded); indP is the public indicator of
// the real slots, so (ind - t)(ind - l) etc. are zero on the padding.
static double chi2MPC(const MPCKeys& mpcKeys,
                      Ciphertext<DCRTPoly> tC,
                      Ciphertext<DCRTPoly> lC,
                      size_t nEff,
                      const Plaintext& indP) {
    auto cc = mpcKeys.cc;

    // cell counts a, b, c, d (slot 0)
    auto notT = cc->EvalSub(indP, tC);
    auto notL = cc->EvalSub(indP, lC);
    auto aC = sumAllSlotsMPC(mpcKeys, cc->EvalMult(notT, notL));
    auto bC = sumAllSlotsMPC(mpcKeys, cc->EvalMult(notT, lC));
    auto cC = sumAllSlotsMPC(mpcKeys, cc->EvalMult(tC, notL));
    auto dC = sumAllSlotsMPC(mpcKeys, cc->EvalMult(tC, lC));

    // chi2 = N (ad - bc)^2 / (r0 r1 c0 c1). We compute with the proportions a/N, ...
    // so that everything stays below 1. (Dividing by N^4 only at the end does not
    // work: OpenFHE rounds the constant v/N^4 to a few units of the scaling factor.)
    const double N = static_cast<double>(nEff);
    const double invN = 1.0 / N;
    auto paC = cc->EvalMult(aC, invN);                      // level 2: p_a = a/N
    auto pbC = cc->EvalMult(bC, invN);
    auto pcC = cc->EvalMult(cC, invN);
    auto pdC = cc->EvalMult(dC, invN);

    auto adC  = cc->EvalMult(paC, pdC);                     // level 3
    auto bcC  = cc->EvalMult(pbC, pcC);
    auto diff = cc->EvalSub(adC, bcC);                      // (ad - bc) / N^2
    auto sqC  = cc->EvalMult(diff, diff);                   // level 4: (ad - bc)^2 / N^4

    auto r0C = cc->EvalAdd(paC, pbC);                       // margins as proportions
    auto r1C = cc->EvalAdd(pcC, pdC);
    auto c0C = cc->EvalAdd(paC, pcC);
    auto c1C = cc->EvalAdd(pbC, pdC);
    auto qC  = cc->EvalMult(cc->EvalMult(r0C, r1C), cc->EvalMult(c0C, c1C));   // level 4: r0 r1 c0 c1 / N^4

    // masked release (v (ad-bc)^2/N^4, v r0 r1 c0 c1/N^4) with one random mask v > 0.
    // The ratio is phi^2 = chi2 / N and both values are at most v/16 (level 5).
    const double v = sampleMask();
    auto maskedP = cc->EvalMult(sqC, v);
    auto maskedQ = cc->EvalMult(qC, v);

    const double mP = decryptMPCScalar(mpcKeys, maskedP);
    const double mQ = decryptMPCScalar(mpcKeys, maskedQ);

    if (g_verbose) {
        std::ostringstream os;
        os << std::setprecision(6) << "    masked release (chi2): " << mP << ", " << mQ
           << "  ->  chi2 = " << N * mP / mQ;
        std::cout << os.str() << std::endl;
    }

    // chi2 = N * phi^2, the mask cancels
    return N * mP / mQ;
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

// Extra levels of the ranking context on top of the comparison depth: 3 for the
// matrix masks before the bootstrapping and a few spare towers for
// IntMPBootAdjustScale (with 12 or 13 in total it fails, 14 works).
static const usint kWilcoxonExtraLevels = 6;

struct Options {
    size_t parties = 3;
    size_t n = 512;          // synthetic cohort size
    size_t runs = 50;
    usint depth = 8;         // Pearson context (4 levels are used)
    usint intBits = 1;       // Pearson: 1+35 bits, ring dim 16384
    usint decBits = 35;
    usint chi2Depth = 5;     // chi-squared context (5 levels are used)
    usint chi2IntBits = 10;  // chi-squared: 10+50 bits, ring dim 16384
    usint chi2DecBits = 50;
    usint compareDepth = 8;  // comparison depth of the ranking (Wilcoxon)
    usint wilcoxonDepth = 0; // ranking context depth; 0 = compareDepth + kWilcoxonExtraLevels
    usint rankIntBits = 1;   // ranking: 1+40 bits, ring dim 32768
    usint rankDecBits = 40;
    std::string csv;
    bool wilcoxon = true;
    bool verbose = false;
};

static void usage() {
    std::cout << "usage: mpc-statistics [--parties 2|3] [--n N] [--runs R] [--depth D] [--intbits I] [--decbits F]\n"
                 "                      [--chi2-depth E] [--chi2-intbits I] [--chi2-decbits F]\n"
                 "                      [--compare-depth C] [--wilcoxon-depth B] [--rank-intbits I] [--rank-decbits F]\n"
                 "                      [--csv FILE] [--no-wilcoxon] [--verbose]\n"
                 "       mpc-statistics 2party|3party            (legacy form)\n";
}

static Options parseArgs(int argc, char* argv[]) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        std::string a(argv[i]);
        auto next = [&](const char* what) -> std::string {
            if (i + 1 >= argc) { std::cerr << "missing value for " << what << std::endl; usage(); std::exit(2); }
            return std::string(argv[++i]);
        };
        if (a == "2party" || a == "2") o.parties = 2;
        else if (a == "3party" || a == "3") o.parties = 3;
        else if (a == "--parties") o.parties = std::stoul(next("--parties"));
        else if (a == "--n") o.n = std::stoul(next("--n"));
        else if (a == "--runs") o.runs = std::stoul(next("--runs"));
        else if (a == "--depth") o.depth = static_cast<usint>(std::stoul(next("--depth")));
        else if (a == "--chi2-depth") o.chi2Depth = static_cast<usint>(std::stoul(next("--chi2-depth")));
        else if (a == "--chi2-intbits") o.chi2IntBits = static_cast<usint>(std::stoul(next("--chi2-intbits")));
        else if (a == "--chi2-decbits") o.chi2DecBits = static_cast<usint>(std::stoul(next("--chi2-decbits")));
        else if (a == "--compare-depth") o.compareDepth = static_cast<usint>(std::stoul(next("--compare-depth")));
        else if (a == "--wilcoxon-depth") o.wilcoxonDepth = static_cast<usint>(std::stoul(next("--wilcoxon-depth")));
        else if (a == "--intbits") o.intBits = static_cast<usint>(std::stoul(next("--intbits")));
        else if (a == "--decbits") o.decBits = static_cast<usint>(std::stoul(next("--decbits")));
        else if (a == "--rank-intbits") o.rankIntBits = static_cast<usint>(std::stoul(next("--rank-intbits")));
        else if (a == "--rank-decbits") o.rankDecBits = static_cast<usint>(std::stoul(next("--rank-decbits")));
        else if (a == "--csv") o.csv = next("--csv");
        else if (a == "--no-wilcoxon") o.wilcoxon = false;
        else if (a == "--verbose") o.verbose = true;
        else if (a == "-h" || a == "--help") { usage(); std::exit(0); }
        else { std::cerr << "unknown argument: " << a << std::endl; usage(); std::exit(2); }
    }
    if (o.parties != 2 && o.parties != 3) { std::cerr << "--parties must be 2 or 3" << std::endl; std::exit(2); }
    if (o.wilcoxonDepth == 0) o.wilcoxonDepth = o.compareDepth + kWilcoxonExtraLevels;
    return o;
}

struct Result {
    std::string name;
    std::string kind;
    size_t nEff;
    double plain;
    double enc;
    double avgMs;
};

static void printResult(const Result& r) {
    std::cout << std::fixed;
    std::cout << "=== " << r.name << " ===" << std::endl;
    std::cout << "  n_eff = " << r.nEff << std::endl;
    std::cout << "  plaintext = " << std::setprecision(6) << r.plain << std::endl;
    std::cout << "  encrypted = " << std::setprecision(6) << r.enc << std::endl;
    if (r.kind == "pearson")
        std::cout << "  deviation (absolute, in r) = " << std::setprecision(6) << std::fabs(r.enc - r.plain) << std::endl;
    else
        std::cout << "  deviation (relative, %)    = " << std::setprecision(3)
                  << 100.0 * std::fabs(r.enc - r.plain) / std::fabs(r.plain) << std::endl;
    std::cout << "  avg time (ms)              = " << std::setprecision(2) << r.avgMs << std::endl << std::endl;
}

int main(int argc, char* argv[]) {
    Options opt = parseArgs(argc, argv);
    g_verbose = opt.verbose;

    // --- data -----------------------------------------------------------------
    Cohort cohort = opt.csv.empty() ? syntheticCohort(opt.n, opt.parties)
                                    : loadCohortCsv(opt.csv, opt.parties);   // rows -> center (i mod parties)
    const size_t numCenters = cohort.numCenters;

    size_t maxN = 0;
    for (const auto& a : cohort.analyses) maxN = std::max(maxN, a.nEff);
    const size_t slots = nextPow2(maxN);   // padded vector length (public)

    std::cout << "Selected mode: " << opt.parties << "-party threshold HE, "
              << (opt.csv.empty() ? "synthetic cohort" : "cohort file " + opt.csv)
              << ", " << numCenters << " emulated centers, vector length " << slots
              << ", runs = " << opt.runs << std::endl;
    for (const auto& a : cohort.analyses)
        std::cout << "  analysis: " << a.name << "  (n_eff = " << a.nEff << ")" << std::endl;

    // Wilcoxon only up to 64 values (the ranking needs n^2 slots)
    const bool runWilcoxon = opt.wilcoxon && slots <= 64;
    if (opt.wilcoxon && !runWilcoxon)
        std::cout << "Wilcoxon skipped: vector length " << slots << " exceeds the ranking primitive's limit (64)." << std::endl;

    // --- keys ------------------------------------------------------------------
    // Context A: batch size = vector length, depth --depth (Pearson).
    MPCKeys keysA = generateMPCKeys(opt.depth, opt.intBits, opt.decBits, slots, slots, opt.parties, false);
    // Context C: batch size = vector length, depth --chi2-depth, 10+50 bits (chi-squared).
    MPCKeys keysC = generateMPCKeys(opt.chi2Depth, opt.chi2IntBits, opt.chi2DecBits, slots, slots, opt.parties, false);
    // Context B: batch size = (vector length)^2 with matrix rotation keys (Wilcoxon ranking).
    MPCKeys keysB;
    if (runWilcoxon)
        keysB = generateMPCKeys(opt.wilcoxonDepth, opt.rankIntBits, opt.rankDecBits, slots, slots * slots, opt.parties, true);

    // --- analyses --------------------------------------------------------------
    std::vector<Result> results;
    for (const auto& a : cohort.analyses) {
        if (a.kind == "wilcoxon" && !runWilcoxon) continue;
        const MPCKeys& keys = (a.kind == "wilcoxon") ? keysB : (a.kind == "chi2" ? keysC : keysA);

        std::vector<double> x = padTo(a.x, slots);
        std::vector<double> y = padTo(a.y, slots);
        Ciphertext<DCRTPoly> xC = encryptByCenter(keys, x, a.center, numCenters);
        Ciphertext<DCRTPoly> yC = encryptByCenter(keys, y, a.center, numCenters);
        Plaintext indP = indicatorPlaintext(keys, a.nEff);

        double plain = 0.0, enc = 0.0, total = 0.0;
        if (a.kind == "pearson")       plain = pearsonPlain(a.x, a.y);
        else if (a.kind == "wilcoxon") plain = wilcoxonPlain(a.x, a.y);
        else                           plain = chi2Plain(a.x, a.y);

        std::cout << "Running " << a.name << " ..." << std::endl;
        for (size_t r = 0; r < opt.runs; ++r) {
            std::pair<double, double> tr;
            if (a.kind == "pearson")
                tr = timed([&] { return pearsonMPC(keys, xC, yC, a.nEff, indP); });
            else if (a.kind == "wilcoxon")
                tr = timed([&] { return wilcoxonMPC(keys, xC, yC, a.nEff, slots, opt.compareDepth); });
            else
                tr = timed([&] { return chi2MPC(keys, xC, yC, a.nEff, indP); });   // x = technique, y = label
            total += tr.first;
            enc = tr.second;
        }
        results.push_back({a.name, a.kind, a.nEff, plain, enc, 1000.0 * total / opt.runs});
        printResult(results.back());
    }

    // --- summary, one line per analysis -----------------------------------------
    std::cout << "SUMMARY parties=" << opt.parties << " slots=" << slots << " depth=" << opt.depth
              << " bits=" << opt.intBits << "+" << opt.decBits
              << " chi2-depth=" << opt.chi2Depth << " chi2-bits=" << opt.chi2IntBits << "+" << opt.chi2DecBits
              << " compare-depth=" << opt.compareDepth << " wilcoxon-depth=" << (runWilcoxon ? opt.wilcoxonDepth : 0)
              << " rank-bits=" << opt.rankIntBits << "+" << opt.rankDecBits << std::endl;
    std::cout << std::fixed;
    for (const auto& r : results) {
        double dev = (r.kind == "pearson") ? std::fabs(r.enc - r.plain)
                                           : 100.0 * std::fabs(r.enc - r.plain) / std::fabs(r.plain);
        std::cout << "  " << std::left << std::setw(60) << r.name << std::right
                  << " n=" << std::setw(4) << r.nEff
                  << " plain=" << std::setw(10) << std::setprecision(4) << r.plain
                  << " enc=" << std::setw(10) << std::setprecision(4) << r.enc
                  << " dev=" << std::setw(8) << std::setprecision(4) << dev << (r.kind == "pearson" ? " (abs)" : " (%)  ")
                  << " ms=" << std::setw(9) << std::setprecision(2) << r.avgMs << std::endl;
    }
    return 0;
}
