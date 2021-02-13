#pragma once

#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <vector>

#include <Fourier.h>

using FilePtr = std::unique_ptr<FILE, decltype(&fclose)>;

struct WavHeader
{
    /* Resource Interchange File Format Symbols -> 4 bytes */
    uint8_t  chunkID[4];

    /* Rest file size without chunkID and chunckSize */
    uint32_t chunkSize;

    /* Info about the format. Consists of WAVE symbols */
    uint8_t  format[4];

    /* Chunks info: fmt symbols and PCM format*/
    uint8_t  subchunk1Id[4];
    uint32_t subchunk1Size;

    uint16_t audioFormat;
    uint16_t numChannels;

    uint32_t sampleRate;
    uint32_t byteRate;

    uint16_t blockAlign;
    uint16_t bitsPerSample;
    uint8_t  subchunk2Id[4];
    uint32_t subchunk2Size;

    /*..data next..*/
};

class WavCompressor {
private:
    const double compressionRate = 0.1;
    std::vector<uint8_t> data;
    WavHeader header_;

    void CompressData();
public:
    explicit WavCompressor(const std::string& fileName);
    void Compress(const std::string& fileName);
};
