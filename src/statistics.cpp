#include "statistics.h"

Ciphertext<DCRTPoly> sumVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length) {
    // Create a vector of ones with the same length
    std::vector<double> ones(length, 1.0);
    auto onesC = cc->Encrypt(pk, cc->MakeCKKSPackedPlaintext(ones));
    
    // Multiply the input vector with ones and sum
    auto prod = cc->EvalMult(c, onesC);
    return prod;
}

Ciphertext<DCRTPoly> sum(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length) {
      auto sum = c;    
    for (size_t i = 1; i < length; i *= 2) {
        auto rotated = cc->EvalRotate(sum, i);
        sum = cc->EvalAdd(sum, rotated);
    }
    return sum;

}

Ciphertext<DCRTPoly> meanVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length) {
    // Sum all slots using rotations
    auto sum = c;    
    for (size_t i = 1; i < length; i *= 2) {
        auto rotated = cc->EvalRotate(sum, i);
        sum = cc->EvalAdd(sum, rotated);
    }

    // Use scalar multiplication instead of fresh encryption for better performance
    double divisor = 1.0 / static_cast<double>(length);
    return cc->EvalMult(sum, divisor);
}

Ciphertext<DCRTPoly> dotProduct(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> a, Ciphertext<DCRTPoly> b, size_t length) {
    auto prod = cc->EvalMult(a, b);
    return sum(cc, pk, prod, length);
}

double decryptScalar(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& sk, Ciphertext<DCRTPoly> c) {
    Plaintext p;
    cc->Decrypt(sk, c, &p);
    p->SetLength(1);
    return p->GetRealPackedValue()[0];
}

Ciphertext<DCRTPoly> centerVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length) {
    auto meanC = meanVector(cc, pk, c, length);
    return cc->EvalSub(c, meanC);
}

