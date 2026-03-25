#include "ranking-mpc.h"
#include "utils-basics.h"
#include "utils-eval.h"
#include "utils-matrices.h"
#include "utils-ptxt.h"
#include <cassert>
#include <omp.h>

// Note: Using decryptMPC and decryptMPCScalar functions from mpc-statistics.cpp

Ciphertext<DCRTPoly> rankMPC(
    const MPCKeys& mpcKeys,
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const bool cmpGt
)
{
    // Same logic as original rank function but return encrypted result
    if (!cmpGt)
    {
        c = compare(
            replicateRow(c, vectorLength),
            replicateColumn(transposeRow(c, vectorLength, true), vectorLength),
            leftBoundC, rightBoundC, degreeC
        );
    }
    else
    {
        c = compareGt(
            replicateRow(c, vectorLength),
            replicateColumn(transposeRow(c, vectorLength, true), vectorLength),
            leftBoundC, rightBoundC, degreeC,
            0.005
        );
    }

    // Sum across columns to obtain per-row ranks, and keep only the first row
    // so that ranks are packed in the first `vectorLength` slots.
    c = sumRows(c, vectorLength, /*maskOutput=*/true);
    c = c + (!cmpGt ? 0.5 : 1.0);

    // Note: MPC scaling should be handled by the calling function
    // since we can't access mpcKeys.cc here due to forward declaration

    return c;
}

std::vector<Ciphertext<DCRTPoly>> rankMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const bool cmpGt,
    const bool complOpt
)
{
    // Call the original rank function to get encrypted results and return them
    return rank(c, subVectorLength, leftBoundC, rightBoundC, degreeC, cmpGt, complOpt);
}

Ciphertext<DCRTPoly> rankWithCorrectionMPC(
    const MPCKeys& mpcKeys,
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC,
    const bool parallel
)
{
    // Call the original rankWithCorrection function and return encrypted result
    return rankWithCorrection(c, vectorLength, leftBoundC, rightBoundC, degreeC, parallel);
}

std::vector<Ciphertext<DCRTPoly>> rankWithCorrectionMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const double leftBoundC,
    const double rightBoundC,
    const uint32_t degreeC
)
{
    // Call the original rankWithCorrection function and return encrypted results
    return rankWithCorrection(c, subVectorLength, leftBoundC, rightBoundC, degreeC);
}

std::vector<Ciphertext<DCRTPoly>> rankFGMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const uint32_t dg,
    const uint32_t df,
    const bool cmpGt,
    const bool complOpt
)
{
    // Call the original rankFG function and return encrypted results
    return rankFG(c, subVectorLength, dg, df, cmpGt, complOpt);
}

Ciphertext<DCRTPoly> rankWithCorrectionFGMPC(
    const MPCKeys& mpcKeys,
    Ciphertext<DCRTPoly> c,
    const size_t vectorLength,
    const uint32_t dg,
    const uint32_t df
)
{
    // Call the original rankWithCorrectionFG function and return encrypted result
    return rankWithCorrectionFG(c, vectorLength, dg, df);
}

std::vector<Ciphertext<DCRTPoly>> rankWithCorrectionFGMPC(
    const MPCKeys& mpcKeys,
    const std::vector<Ciphertext<DCRTPoly>>& c,
    const size_t subVectorLength,
    const uint32_t dg,
    const uint32_t df
)
{
    // Call the original rankWithCorrectionFG function and return encrypted results
    return rankWithCorrectionFG(c, subVectorLength, dg, df);
}