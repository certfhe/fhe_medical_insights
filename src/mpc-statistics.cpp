#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "statistics.h"
#include "ranking.h"
#include "ranking-mpc.h"

#include <chrono>
#include <numeric>
#include <cmath>
#include <iostream>
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

struct MPCKeys {
    CryptoContext<DCRTPoly> cc;
    KeyPair<DCRTPoly> party1Keys;
    KeyPair<DCRTPoly> party2Keys;
    KeyPair<DCRTPoly> party3Keys;
    PublicKey<DCRTPoly> jointPublicKey;
    usint integralPrecision;    // Bits for integer part
    usint decimalPrecision;     // Bits for fractional part
    size_t numParties = 3;      // 2 or 3 parties
};

// MPC Decryption helper
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

// MPC Decryption helper
std::vector<double> decryptMPC(const MPCKeys& mpcKeys, const Ciphertext<DCRTPoly>& ciphertext) {
    // Party 1 performs lead decryption
    auto ciphertextPartial1 = mpcKeys.cc->MultipartyDecryptLead({ciphertext}, mpcKeys.party1Keys.secretKey);

    // Party 2 performs main decryption
    auto ciphertextPartial2 = mpcKeys.cc->MultipartyDecryptMain({ciphertext}, mpcKeys.party2Keys.secretKey);

    // Combine partial decryptions
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

    return plaintextResult->GetRealPackedValue();
}

// Interactive MPC bootstrapping (TCKKS IntMPBoot) for 3 parties
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

// MPC Key Generation for 2 parties with decimal/fractional precision
MPCKeys generateMPCKeys(const usint compareDepth,
                       const usint integralPrecision,
                       const usint decimalPrecision,
                       const size_t n,
                       const size_t slots,
                       const size_t numParties) {
    std::cout << "Setting up " << numParties << "-Party MPC keys with precision..." << std::endl;
    std::cout << "Integral precision: " << integralPrecision << " bits" << std::endl;
    std::cout << "Decimal precision: " << decimalPrecision << " bits" << std::endl;

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

    // Print detailed crypto context information
    const BigInteger ciphertextModulus = cc->GetModulus();
    const usint ciphertextModulusBitsize = ciphertextModulus.GetLengthForBase(2);
    const usint ringDimension = cc->GetRingDimension();
    const usint maxNumSlots = ringDimension / 2;
    const auto elementParameters = cc->GetCryptoParameters()->GetElementParams()->GetParams();
    std::vector<NativeInteger> moduliChain(elementParameters.size());
    std::vector<usint> moduliChainBitsize(elementParameters.size());
    for (size_t i = 0; i < elementParameters.size(); i++) {
        moduliChain[i] = elementParameters[i]->GetModulus();
        moduliChainBitsize[i] = moduliChain[i].GetLengthForBase(2);
    }

    std::ostringstream logMessage;
    logMessage << "CKKS PARAMETERS"                                                              << std::endl;
    logMessage << "Integral Bit Precision        : " << integralPrecision                        << std::endl;
    logMessage << "Decimal Bit Precision         : " << decimalPrecision                         << std::endl;
    logMessage << "Ciphertext Modulus Precision  : " << ciphertextModulusBitsize                 << std::endl;
    logMessage << "Ring Dimension                : " << ringDimension                            << std::endl;
    logMessage << "Max Slots                     : " << maxNumSlots                              << std::endl;
    logMessage << "Slots                         : " << slots                                    << std::endl;
    logMessage << "Multiplicative Depth          : " << parameters.GetMultiplicativeDepth()      << std::endl;
    logMessage << "Security Level                : " << parameters.GetSecurityLevel()            << std::endl;
    logMessage << "Secret Key Distribution       : " << parameters.GetSecretKeyDist()            << std::endl;
    logMessage << "Scaling Technique             : " << parameters.GetScalingTechnique()         << std::endl;
    logMessage << "Encryption Technique          : " << parameters.GetEncryptionTechnique()      << std::endl;
    logMessage << "Multiplication Technique      : " << parameters.GetMultiplicationTechnique()  << std::endl;
    logMessage << "Moduli Chain Bitsize          : " << moduliChainBitsize                       << std::endl;
    logMessage << std::endl;

    std::cout << logMessage.str();

    MPCKeys mpcKeys;
    mpcKeys.cc = cc;
    mpcKeys.integralPrecision = integralPrecision;
    mpcKeys.decimalPrecision = decimalPrecision;
    mpcKeys.numParties = numParties;

    std::cout << "Party 1: Generating initial keys..." << std::endl;

    // Round 1: Party 1 generates initial keys
    mpcKeys.party1Keys = cc->KeyGen();

    // Generate evaluation keys for Party 1 (same as original keyGeneration)
    cc->EvalMultKeyGen(mpcKeys.party1Keys.secretKey);
    cc->EvalSumKeyGen(mpcKeys.party1Keys.secretKey);

    // Generate rotation keys for Party 1 (needed for ranking operations)
    std::vector<int32_t> indices = getRotationIndices(n);
    cc->EvalRotateKeyGen(mpcKeys.party1Keys.secretKey, indices);

    auto evalSumKeys1 = std::make_shared<std::map<usint, EvalKey<DCRTPoly>>>(
        cc->GetEvalSumKeyMap(mpcKeys.party1Keys.secretKey->GetKeyTag()));

    std::cout << "Party 2: Contributing to key generation..." << std::endl;
    // Party 2 contributes
    mpcKeys.party2Keys = cc->MultipartyKeyGen(mpcKeys.party1Keys.publicKey);

    // Generate evaluation keys for Party 2
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
        // Insert 2-party eval sum and rotation keys
        cc->InsertEvalSumKey(evalSumKeys12);
        cc->InsertEvalAutomorphismKey(evalRotateKeys12);

        std::cout << "Finalizing joint keys (2-party)..." << std::endl;
        // Finalize eval mult for 2 parties
        auto evalMultAB = cc->MultiAddEvalKeys(evalMultKey1, evalMultKey2, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultAAB = cc->MultiMultEvalKey(mpcKeys.party1Keys.secretKey, evalMultAB, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultBAB = cc->MultiMultEvalKey(mpcKeys.party2Keys.secretKey, evalMultAB, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultFinal = cc->MultiAddEvalMultKeys(evalMultAAB, evalMultBAB, evalMultAAB->GetKeyTag());
        cc->InsertEvalMultKey({evalMultFinal});

        // Use last party's public key as the joint encryption key (matches eval keys' tag)
        mpcKeys.jointPublicKey = mpcKeys.party2Keys.publicKey;

        std::cout << "2-Party MPC key generation completed!" << std::endl;
        std::cout << "Total precision: " << totalPrecision << " bits" << std::endl;
        return mpcKeys;
    } else {
        std::cout << "Party 3: Contributing to key generation..." << std::endl;
        // Party 3 contributes
        mpcKeys.party3Keys = cc->MultipartyKeyGen(mpcKeys.party2Keys.publicKey);
        cc->EvalMultKeyGen(mpcKeys.party3Keys.secretKey);
        cc->EvalSumKeyGen(mpcKeys.party3Keys.secretKey);
        cc->EvalRotateKeyGen(mpcKeys.party3Keys.secretKey, indices);

        // For multiparty, derive evalMult share for party3 using the original base key share
        auto evalMultKey3 = cc->MultiKeySwitchGen(mpcKeys.party3Keys.secretKey, mpcKeys.party3Keys.secretKey, evalMultKey1);
        auto evalMultAB = cc->MultiAddEvalKeys(evalMultKey1, evalMultKey2, mpcKeys.party2Keys.publicKey->GetKeyTag());
        auto evalMultABC = cc->MultiAddEvalKeys(evalMultAB, evalMultKey3, mpcKeys.party3Keys.publicKey->GetKeyTag());

        // Merge party3 eval sum
        auto evalSumKeys3 = cc->MultiEvalSumKeyGen(mpcKeys.party3Keys.secretKey, evalSumKeys1, mpcKeys.party3Keys.publicKey->GetKeyTag());
        auto evalSumKeys123 = cc->MultiAddEvalSumKeys(evalSumKeys12, evalSumKeys3, mpcKeys.party3Keys.publicKey->GetKeyTag());
        cc->InsertEvalSumKey(evalSumKeys123);

        // Merge party3 rotations
        auto evalRotateKeys3 = cc->MultiEvalAtIndexKeyGen(mpcKeys.party3Keys.secretKey, evalRotateKeys1, indices, mpcKeys.party3Keys.publicKey->GetKeyTag());
        auto evalRotateKeys123 = cc->MultiAddEvalAutomorphismKeys(evalRotateKeys12, evalRotateKeys3, mpcKeys.party3Keys.publicKey->GetKeyTag());
        cc->InsertEvalAutomorphismKey(evalRotateKeys123);

        std::cout << "Finalizing joint keys (3-party)..." << std::endl;
        // Finalize eval mult for 3 parties (compute shares for each and add)
        auto evalMultAABC = cc->MultiMultEvalKey(mpcKeys.party1Keys.secretKey, evalMultABC, mpcKeys.party3Keys.publicKey->GetKeyTag());
        auto evalMultBABC = cc->MultiMultEvalKey(mpcKeys.party2Keys.secretKey, evalMultABC, mpcKeys.party3Keys.publicKey->GetKeyTag());
        auto evalMultCABC = cc->MultiMultEvalKey(mpcKeys.party3Keys.secretKey, evalMultABC, mpcKeys.party3Keys.publicKey->GetKeyTag());
        auto evalMultAB_ABC = cc->MultiAddEvalMultKeys(evalMultAABC, evalMultBABC, evalMultAABC->GetKeyTag());
        auto evalMultFinal  = cc->MultiAddEvalMultKeys(evalMultAB_ABC, evalMultCABC, evalMultAB_ABC->GetKeyTag());
        cc->InsertEvalMultKey({evalMultFinal});

        // Use last party's public key as the joint encryption key (matches eval keys' tag)
        mpcKeys.jointPublicKey = mpcKeys.party3Keys.publicKey;

        std::cout << "3-Party MPC key generation completed!" << std::endl;
        std::cout << "Total precision: " << totalPrecision << " bits" << std::endl;
        return mpcKeys;
    }
}

// MPC-compatible version that avoids fresh encryptions and supports MPC decryption
Ciphertext<DCRTPoly> sumMPC(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length, const MPCKeys& mpcKeys) {
    // Sum all slots using rotations (same logic)
    auto sum = c;
    for (size_t i = 1; i < length; i *= 2) {
        auto rotated = cc->EvalRotate(sum, i);
        sum = cc->EvalAdd(sum, rotated);
    }

    return sum;
}

Ciphertext<DCRTPoly> dotProductMPC(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> a, Ciphertext<DCRTPoly> b, size_t length, const MPCKeys& mpcKeys) {
    auto prod = cc->EvalMult(a, b);
    return sumMPC(cc, pk, prod, length, mpcKeys);
}

// MPC-compatible version that avoids fresh encryptions and supports MPC decryption
Ciphertext<DCRTPoly> meanVectorMPC(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length, const MPCKeys& mpcKeys) {
    // Sum all slots using rotations (same logic)
    auto sum = c;
    for (size_t i = 1; i < length; i *= 2) {
        auto rotated = cc->EvalRotate(sum, i);
        sum = cc->EvalAdd(sum, rotated);
    }

    // Use scalar multiplication instead of fresh encryption for MPC compatibility
    double divisor = 1.0 / static_cast<double>(length);
    return cc->EvalMult(sum, divisor);
}

Ciphertext<DCRTPoly> centerVectorMPC(const CryptoContext<DCRTPoly>& cc,
                                    const PublicKey<DCRTPoly>& pk,
                                    Ciphertext<DCRTPoly> c,
                                    size_t length,
                                    const MPCKeys& mpcKeys,
                                    const Ciphertext<DCRTPoly>& maskC) {
    // Step 1: Compute mean (mean value is in slot 0, garbage in other slots)
    auto meanC = meanVectorMPC(cc, pk, c, length, mpcKeys);
    cc->RescaleInPlace(meanC);
    cc->RescaleInPlace(c);

    return cc->EvalSub(c, meanC);
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

// MPC Pearson correlation
static double pearsonMPC(const MPCKeys& mpcKeys,
                        Ciphertext<DCRTPoly> xC,
                        Ciphertext<DCRTPoly> yC,
                        size_t n,
                        const Ciphertext<DCRTPoly>& maskC) {

    auto meanXC = meanVectorMPC(mpcKeys.cc, mpcKeys.jointPublicKey, xC, n, mpcKeys);
    auto meanYC = meanVectorMPC(mpcKeys.cc, mpcKeys.jointPublicKey, yC, n, mpcKeys);

    auto xCen = centerVectorMPC(mpcKeys.cc, mpcKeys.jointPublicKey, xC, n, mpcKeys, maskC);
    auto yCen = centerVectorMPC(mpcKeys.cc, mpcKeys.jointPublicKey, yC, n, mpcKeys, maskC);

    auto xyC = dotProductMPC(mpcKeys.cc, mpcKeys.jointPublicKey, xCen, yCen, n, mpcKeys);
    xyC = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, xyC, n, mpcKeys);
    auto x2C = dotProductMPC(mpcKeys.cc, mpcKeys.jointPublicKey, xCen, xCen, n, mpcKeys);
    x2C = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, x2C, n, mpcKeys);
    auto y2C = dotProductMPC(mpcKeys.cc, mpcKeys.jointPublicKey, yCen, yCen, n, mpcKeys);
    y2C = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, y2C, n, mpcKeys);

    // Apply masking for security
    double u1 = 1.23;
    double u2 = 2.34;
    double u3 = u1 / u2;

    auto maskedXdotY = mpcKeys.cc->EvalMult(xyC, u1);
    auto maskedLenX = mpcKeys.cc->EvalMult(x2C, u2*u2);
    auto maskedLenY = mpcKeys.cc->EvalMult(y2C, u3*u3);

    // MPC Decryption
    double masked_dot = decryptMPCScalar(mpcKeys, maskedXdotY);
    double masked_x2 = decryptMPCScalar(mpcKeys, maskedLenX);
    double masked_y2 = decryptMPCScalar(mpcKeys, maskedLenY);

    double numerator = masked_dot;
    double denominator = std::sqrt(masked_x2) * std::sqrt(masked_y2);
    double correction = u1 / (u2 * u3);

    return correction * numerator / denominator;
}

// Plaintext Wilcoxon rank-sum statistic
static double wilcoxonPlain(const std::vector<double>& x, const std::vector<double>& g) {
    auto ranks = rank(x);
    double S1 = 0.0;
    double n1 = 0.0;

    for (size_t i = 0; i < x.size(); ++i) {
        if (g[i] > 0.5) {
            S1 += ranks[i];
            n1 += 1.0;
        }
    }

    double n0 = static_cast<double>(x.size()) - n1;
    double U1 = S1 - (n1 * (n1 + 1)) / 2.0;
    return (U1 - n0 * n1 / 2.0) / std::sqrt(n0 * n1 * (n0 + n1 + 1) / 12.0);
}

// MPC Wilcoxon rank-sum statistic
static double wilcoxonMPC(const MPCKeys& mpcKeys,
                         Ciphertext<DCRTPoly> xC,
                         Ciphertext<DCRTPoly> gC,
                         size_t n,
                         usint compareDepth) {
    Ciphertext<DCRTPoly> rankC = rankMPC(mpcKeys, xC, n, -1.0, 1.0, depth2degree(compareDepth));

    // SumRows in rankMPC produces a matrix where each row i contains rank[i]
    // replicated across columns. Extract the first row to obtain a row vector
    // packed in the first n slots.
    Ciphertext<DCRTPoly> rankVec = maskRow(rankC, n, 0);

    // MPC interactive bootstrapping on rankVec before using it
    rankVec = mpcBootstrap(mpcKeys, rankVec);

    auto s1C = dotProductMPC(mpcKeys.cc, mpcKeys.jointPublicKey, rankVec, gC, n, mpcKeys);
    auto n1C = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, gC, n,mpcKeys);

    double S1 = decryptMPCScalar(mpcKeys, s1C);
    double n1 = decryptMPCScalar(mpcKeys, n1C);
    double n0 = static_cast<double>(n) - n1;

    double U1 = S1 - (n1 * (n1 + 1)) / 2.0;
    return (U1 - (n0 * n1) / 2.0) / std::sqrt((n0 * n1 * (n0 + n1 + 1)) / 12.0);
}

// Plaintext chi-squared for 2x2 contingency table
static double chi2Plain(const std::vector<double>& t, const std::vector<double>& l) {
    double a=0, b=0, c=0, d=0;
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
    double N = r0 + r1;
    double ea = r0 * c0 / N;
    double eb = r0 * c1 / N;
    double ec = r1 * c0 / N;
    double ed = r1 * c1 / N;
    return (std::pow(a-ea,2)/ea) + (std::pow(b-eb,2)/eb) +
           (std::pow(c-ec,2)/ec) + (std::pow(d-ed,2)/ed);
}

// MPC chi-squared
static double chi2MPC(const MPCKeys& mpcKeys,
                     Ciphertext<DCRTPoly> tC,
                     Ciphertext<DCRTPoly> lC,
                     Ciphertext<DCRTPoly> oneC,
                     size_t n) {
    auto aC = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, (oneC - tC) * (oneC - lC), n,mpcKeys);
    auto bC = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, (oneC - tC) * lC, n, mpcKeys);
    auto cC = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, tC * (oneC - lC), n, mpcKeys);
    auto dC = sumMPC(mpcKeys.cc, mpcKeys.jointPublicKey, tC * lC, n,mpcKeys);

    // MPC Decryption
    double a = decryptMPCScalar(mpcKeys, aC);
    double b = decryptMPCScalar(mpcKeys, bC);
    double c = decryptMPCScalar(mpcKeys, cC);
    double d = decryptMPCScalar(mpcKeys, dC);

    double r0 = a + b;
    double r1 = c + d;
    double c0 = a + c;
    double c1 = b + d;
    double N = r0 + r1;
    double ea = r0 * c0 / N;
    double eb = r0 * c1 / N;
    double ec = r1 * c0 / N;
    double ed = r1 * c1 / N;
    return (std::pow(a-ea,2)/ea) + (std::pow(b-eb,2)/eb) +
           (std::pow(c-ec,2)/ec) + (std::pow(d-ed,2)/ed);
}

int main(int argc, char* argv[]) {

    // Sample data (512 elements)
    std::vector<double> x = {0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9};
    std::vector<double> y = {0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9,0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9, 0.1, 0.3, 0.7, 0.9};

    // Wilcoxon group vector
    std::vector<double> g = {1, 0, 0, 1};

    // Chi-squared treatment and label vectors (512 elements)
    std::vector<double> t = {1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1  };
    std::vector<double> l = {0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1 , 1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1  } ;

    const size_t n = x.size();
    const usint compareDepth = 8;
    const size_t runs = 50;
    const size_t slots = n;

    // Mode selection: default 3-party, allow CLI arg: "2party" or "3party"
    size_t numParties = 3;
    if (argc > 1) {
        std::string modeArg(argv[1]);
        if (modeArg == "2party" || modeArg == "2") {
            numParties = 2;
        } else if (modeArg == "3party" || modeArg == "3") {
            numParties = 3;
        } else {
            std::cout << "Unknown mode '" << modeArg << "'. Using default 3-party." << std::endl;
        }
    }
    std::cout << "Selected mode: " << numParties << "-party" << std::endl;

    // Precision parameters
    const usint integralPrecision = 1;
    const usint decimalPrecision = 25;
    // Generate MPC keys with precision parameters
    auto mpcKeys = generateMPCKeys(compareDepth, integralPrecision, decimalPrecision, n, slots, numParties);

    // Encryption function using joint public key
    auto enc = [&](const std::vector<double>& v) {
        double scale = 1;  // Fixed scale
        Plaintext pt = mpcKeys.cc->MakeCKKSPackedPlaintext(v, scale, 0, nullptr, slots);
        return mpcKeys.cc->Encrypt(mpcKeys.jointPublicKey, pt);
    };

    std::cout << "Encrypting data..." << std::endl;
    Ciphertext<DCRTPoly> xC = enc(x);
    Ciphertext<DCRTPoly> yC = enc(y);
    Ciphertext<DCRTPoly> gC = enc(g);
    Ciphertext<DCRTPoly> tC = enc(t);
    Ciphertext<DCRTPoly> lC = enc(l);
    std::vector<double> ones(n, 1.0);
    Ciphertext<DCRTPoly> oneC = enc(ones);
    oneC*=oneC; // get rid of the issue with the floor noise which in the case of MPC breaks up the encryption

    // Encrypt the mask
    std::vector<double> mask(n, 0.0);
    mask[0] = 1.0; // [1, 0, 0, 0]
    Ciphertext<DCRTPoly> maskC = enc(mask);

    // Verify encryption
    auto dec = [&](const Ciphertext<DCRTPoly>& c) {
        return decryptMPCScalar(mpcKeys, c);
    };

    std::cout << "Decrypting oneC..." << std::endl;
    std::cout<<"Level: "<<oneC->GetLevel()<<"\nScaling Factor: "<<oneC->GetScalingFactor()<<std::endl;
    std::cout << "oneC: " << dec(oneC) << std::endl;

    std::cout << "Number of values: " << n << std::endl;
    std::cout << "Running " << runs << " iterations for timing..." << std::endl << std::endl;

    double tPlain = 0, tMPC = 0, resPlain = 0, resMPC = 0;

    // Pearson Correlation
    for(size_t i = 0; i < runs; i++) {
        auto [tp, rp] = timed([&]{ return pearsonPlain(x, y); });
        tPlain += tp; resPlain = rp;
        auto [tm, rm] = timed([&]{ return pearsonMPC(mpcKeys, xC, yC, n, maskC); });
        tMPC += tm; resMPC = rm;
    }

    std::cout << "=== Pearson Correlation Results ===" << std::endl;
    std::cout << "Plaintext result: " << resPlain
              << " avg time: " << (tPlain / runs) * 1000 << "ms" << std::endl;
    std::cout << "MPC result: " << resMPC
              << " avg time: " << (tMPC / runs) * 1000 << "ms" << std::endl << std::endl;

    // Wilcoxon Test
    tPlain = tMPC = 0; resPlain = resMPC = 0;
    for(size_t i = 0; i < runs; i++) {
        auto [tp, rp] = timed([&]{ return wilcoxonPlain(x, g); });
        tPlain += tp; resPlain = rp;
        auto [tm, rm] = timed([&]{ return wilcoxonMPC(mpcKeys, xC, gC, n, compareDepth); });
        tMPC += tm; resMPC = rm;
    }

    std::cout << "=== Wilcoxon Test Results ===" << std::endl;
    std::cout << "Plaintext Z: " << resPlain
              << " avg time: " << (tPlain / runs) * 1000 << "ms" << std::endl;
    std::cout << "MPC Z: " << resMPC
              << " avg time: " << (tMPC / runs) * 1000 << "ms" << std::endl << std::endl;

    // Chi-squared Test
    tPlain = tMPC = 0; resPlain = resMPC = 0;
    for(size_t i = 0; i < runs; i++) {
        auto [tp, rp] = timed([&]{ return chi2Plain(t, l); });
        tPlain += tp; resPlain = rp;
        auto [tm, rm] = timed([&]{ return chi2MPC(mpcKeys, tC, lC, oneC, n); });
        tMPC += tm; resMPC = rm;
    }

    std::cout << "=== Chi-squared Test Results ===" << std::endl;
    std::cout << "Plaintext: " << resPlain
              << " avg time: " << (tPlain / runs) * 1000 << "ms" << std::endl;
    std::cout << "MPC: " << resMPC
              << " avg time: " << (tMPC / runs) * 1000 << "ms" << std::endl << std::endl;

    return 0;
}
