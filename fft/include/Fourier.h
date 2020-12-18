#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

using ComplexPoly = std::vector<std::complex<double>>;

struct PolyPair {
    ComplexPoly first;
    ComplexPoly second;
};

class FourierHandler {
private:
    ComplexPoly initial_;
    bool reverse_ = false;

    ComplexPoly GetSecondPowerFFT(ComplexPoly&& current);
    ComplexPoly GetSecondPowerPoly(const ComplexPoly& current);

    size_t getNextPowerOfTwo(size_t number) noexcept;

    PolyPair GetSplitted(const ComplexPoly& current);
    ComplexPoly Merge(const ComplexPoly& oddResult, const ComplexPoly& evenResult);

public:
    [[maybe_unused]] explicit FourierHandler(const ComplexPoly& initial) : initial_(initial) {}
    explicit FourierHandler(ComplexPoly&& initial) : initial_(std::move(initial)) {}

    ComplexPoly GetFFT();
    ComplexPoly GetReverseFFT(const ComplexPoly& fft);
};
