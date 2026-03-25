#pragma once
#include "utils-basics.h"
#include "utils-ptxt.h"

using namespace lbcrypto;

// Sum all slots of a ciphertext vector. Result is replicated across all slots.
Ciphertext<DCRTPoly> sumVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length);

// Compute the average of all slots. Result is replicated across all slots.
Ciphertext<DCRTPoly> meanVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length);

// Element-wise product then sum over all slots.
Ciphertext<DCRTPoly> dotProduct(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> a, Ciphertext<DCRTPoly> b, size_t length);

// Decrypt ciphertext that encodes a single value replicated across slots.
Ciphertext<DCRTPoly> centerVector(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length);

double decryptScalar(const CryptoContext<DCRTPoly>& cc, const PrivateKey<DCRTPoly>& sk, Ciphertext<DCRTPoly> c);

Ciphertext<DCRTPoly> sum(const CryptoContext<DCRTPoly>& cc, const PublicKey<DCRTPoly>& pk, Ciphertext<DCRTPoly> c, size_t length);