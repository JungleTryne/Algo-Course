#include <Fourier.h>
#include <WavReader.h>

int main() {
    WavCompressor compressor("files/thomas.wav");
    compressor.Compress("files/thomas_com.wav");
    return 0;
}

