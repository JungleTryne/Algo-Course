#include <Fourier.h>

/* Fast Fourier Transform algorithm. Mishin Danila */

PolyPair FourierHandler::GetSplitted(const ComplexPoly &current)  {
    size_t size = current.size();

    ComplexPoly odd(size >> 1u);
    ComplexPoly even(size >> 1u);

    for(size_t i = 0; i < (size >> 1u); ++i) {
        even[i] = current[2 * i];
        odd[i] = current[2 * i + 1];
    }

    return {std::move(odd), std::move(even)};
}

ComplexPoly FourierHandler::Merge(const ComplexPoly &oddResult, const ComplexPoly &evenResult) {
    size_t size = oddResult.size() + evenResult.size();

    ComplexPoly finalResult(size);

    double stepAngle = !reverse_ ? (2 * M_PI / size) : (-2 * M_PI / size);
    std::complex<double> stepW{cos(stepAngle), sin(stepAngle)};
    std::complex<double> currentW{cos(0), sin(0)};

    for(size_t i = 0; i < (size >> 1u); ++i) {
        size_t firstPos = i;                    //i
        size_t secondPos = i + (size >> 1u);    //i + n/2

        auto ithFirstResult  = evenResult[i] + currentW * oddResult[i];
        auto ithSecondResult = evenResult[i] - currentW * oddResult[i];

        finalResult[firstPos]  = ithFirstResult;
        finalResult[secondPos] = ithSecondResult;
        currentW *= stepW;
    }

    return finalResult;
}

ComplexPoly FourierHandler::GetSecondPowerFFT(ComplexPoly&& current) {
    size_t size = current.size();
    assert(size > 0 && (size & (size - 1)) == 0); //checking if size = 2^n

    if(size == 1) {
        return std::move(current);
    }

    auto [odd, even] = GetSplitted(current);
    assert(odd.size() == size >> 1u);

    auto oddResult = GetSecondPowerFFT(std::move(odd));
    auto evenResult = GetSecondPowerFFT(std::move(even));
    assert(odd.size() == even.size());

    return Merge(oddResult, evenResult);
}

ComplexPoly FourierHandler::GetFFT() {
    reverse_ = false;
    auto poly = GetSecondPowerPoly(initial_);
    return GetSecondPowerFFT(std::move(poly));
}

ComplexPoly FourierHandler::GetReverseFFT(const ComplexPoly& fft) {
    reverse_ = true;
    auto secondPower = GetSecondPowerPoly(fft);
    auto result = GetSecondPowerFFT(std::move(secondPower));

    std::transform(result.begin(), result.end(), result.begin(), [&](auto num) -> auto {
        return num / static_cast<double>(result.size());
    });

    return result;
}

ComplexPoly FourierHandler::GetSecondPowerPoly(const ComplexPoly& initial) {
    size_t size = initial.size();
    size_t new_size = size;

    if((size & (size - 1)) != 0) {
        new_size = getNextPowerOfTwo(size);
    }

    ComplexPoly polyToProceed(new_size);
    std::copy(initial.begin(), initial.end(), polyToProceed.begin());
    return polyToProceed;
}

size_t FourierHandler::getNextPowerOfTwo(size_t number) noexcept {
    size_t counter = 0;
    while(number != 0) {
        number >>= 1u;
        ++counter;
    }
    return 1u << counter;
}
