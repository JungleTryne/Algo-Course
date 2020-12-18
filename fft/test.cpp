#include <Fourier.h>

void Test() {
    ComplexPoly testing = {1, 2, 3, 4, 5, 6, 7, 8};
    FourierHandler handler(testing);

    ComplexPoly resultFFT = handler.GetFFT();
    ComplexPoly reverseFFT = handler.GetReverseFFT(resultFFT);

    for(size_t i = 0; i < testing.size(); ++i) {
        assert(std::abs(reverseFFT[i] - testing[i]) < 1e-5);
    }

    std::cout << "Fourier reverse test: " << "OK" << std::endl;

    ComplexPoly correctFFT = {
            {36, 0},
            {-4, -9.657},
            {-4, -4},
            {-4, -1.657},
            {-4, 0},
            {-4, 1.657},
            {-4, 4},
            {-4, 9.657}
    };

    for(size_t i = 0; i < correctFFT.size(); ++i) {
        assert(std::abs(resultFFT[i] - correctFFT[i]) < 1e-3);
    }

    std::cout << "Fourier direct test: " << "OK" << std::endl;
}

int main() {
    Test();
    return 0;
}