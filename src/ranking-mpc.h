#pragma once

#include "openfhe.h"
#include "ranking.h"

using namespace lbcrypto;

// Forward declaration of MPCKeys from mpc-statistics.cpp
struct MPCKeys;

/**
 * @brief Computes the rank of elements in a ciphertext vector using MPC.
 * 
 * This function computes the rank of elements in a ciphertext vector `c` using
 * multi-party computation. The result is returned as an encrypted ciphertext.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c The ciphertext vector for which to compute the rank.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param cmpGt Flag indicating whether to compute standard (true) or
 * fractional rank (false), default is false.
 * @return Ciphertext<DCRTPoly> The encrypted ciphertext containing the computed ranks.
 */
Ciphertext<DCRTPoly> rankMPC(
    const MPCKeys& mpcKeys,
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const bool cmpGt = false
);

/**
 * @brief Computes the rank of elements in a vector stored in multiple
 * ciphertexts using MPC.
 * 
 * This function computes the rank of elements in multiple ciphertext vectors 
 * using multi-party computation.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param cmpGt Flag indicating whether to compute standard (true) or
 * fractional rank (false), default is false.
 * @param complOpt Flag indicating whether to use the complementary comparison
 * optimization.
 * @return std::vector<Ciphertext<DCRTPoly>> A vector containing encrypted
 * ciphertexts representing the computed ranks.
 */
std::vector<Ciphertext<DCRTPoly>> rankMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const bool cmpGt = false,
    const bool complOpt = true
);

/**
 * @brief Computes the rank of elements in a ciphertext vector using MPC
 * with tie correction.
 * 
 * This function computes the rank of elements in a ciphertext vector `c`,
 * while handling ties at the same time using multi-party computation.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c The ciphertext vector for which to compute the rank.
 * @param vectorLength The length of the vector.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @param parallel Flag indicating whether to compute the correction offset in
 * parallel or not, default is false.
 * @return Ciphertext<DCRTPoly> The encrypted ciphertext containing the computed ranks.
 */
Ciphertext<DCRTPoly> rankWithCorrectionMPC(
    const MPCKeys& mpcKeys,
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const bool parallel = false
);

/**
 * @brief Computes the rank of elements in multiple ciphertext vectors using MPC
 * with tie correction.
 * 
 * This function computes the rank of elements in multiple ciphertext vectors,
 * while handling ties at the same time using multi-party computation.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param leftBoundC The left bound for comparison's approximation.
 * @param rightBoundC The right bound for comparison's approximation.
 * @param degreeC The degree of the comparison's approximation.
 * @return std::vector<Ciphertext<DCRTPoly>> A vector containing encrypted
 * ciphertexts representing the computed ranks.
 */
std::vector<Ciphertext<DCRTPoly>> rankWithCorrectionMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC
);

/**
 * @brief Computes the rank of elements using FG approximation with MPC.
 * 
 * This function computes the rank of elements in multiple ciphertext vectors
 * using the fg approximation of the sign function and multi-party computation.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param dg The composition degree of g (reduce the input gap).
 * @param df The composition degree of f (reduce the output error).
 * @param cmpGt Flag indicating whether to compute standard (true) or
 * fractional rank (false), default is false.
 * @param complOpt Flag indicating whether to use the complementary comparison
 * optimization.
 * @return std::vector<Ciphertext<DCRTPoly>> A vector containing encrypted
 * ciphertexts representing the computed ranks.
 */
std::vector<Ciphertext<DCRTPoly>> rankFGMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const uint32_t dg,
    const uint32_t df,
    const bool cmpGt = false,
    const bool complOpt = true
);

/**
 * @brief Computes the rank with tie correction using FG approximation and MPC.
 * 
 * This function computes the rank of elements in a ciphertext vector,
 * while handling ties using the fg approximation and multi-party computation.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c The ciphertext vector for which to compute the rank.
 * @param vectorLength The length of the vector.
 * @param dg The composition degree of g (reduce the input gap).
 * @param df The composition degree of f (reduce the output error).
 * @return Ciphertext<DCRTPoly> The encrypted ciphertext containing the computed ranks.
 */
Ciphertext<DCRTPoly> rankWithCorrectionFGMPC(
    const MPCKeys& mpcKeys,
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const uint32_t dg,
    const uint32_t df
);

/**
 * @brief Computes the rank with tie correction using FG approximation and MPC
 * for multiple vectors.
 * 
 * This function computes the rank of elements in multiple ciphertext vectors,
 * while handling ties using the fg approximation and multi-party computation.
 * 
 * @param mpcKeys The MPC keys structure containing joint cryptographic context.
 * @param c A vector containing ciphertexts, each representing a portion of the
 * input vector.
 * @param subVectorLength The length of each sub-vector stored across the
 * ciphertexts.
 * @param dg The composition degree of g (reduce the input gap).
 * @param df The composition degree of f (reduce the output error).
 * @return std::vector<Ciphertext<DCRTPoly>> A vector containing encrypted
 * ciphertexts representing the computed ranks.
 */
std::vector<Ciphertext<DCRTPoly>> rankWithCorrectionFGMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const uint32_t dg,
    const uint32_t df
);

// Note: decryptMPC and decryptMPCScalar functions are defined in mpc-statistics.cpp