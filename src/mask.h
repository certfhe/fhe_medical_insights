#pragma once

#include <cmath>
#include <random>

// Random mask for the client-aided release (Section 3.5 of the paper):
// u = exp(sigma * Z) with Z ~ N(0,1) truncated to |Z| <= 4, sigma = 1.
// A new mask is drawn every time and never stored or reused. u > 0, so signs
// are kept. The leakage bound in the paper depends on sigma, so keep sigma = 1.
inline double sampleMask(double sigma = 1.0)
{
    static std::random_device rd;   // OS entropy, the analyst must not be able to predict the mask
    std::normal_distribution<double> gauss(0.0, 1.0);
    double z;
    do {
        z = gauss(rd);
    } while (std::fabs(z) > 4.0);   // truncate at 4 sigma
    return std::exp(sigma * z);
}
