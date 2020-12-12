#include <WavReader.h>

#include <fstream>

WavCompressor::WavCompressor(const std::string& fileName)
{
    FilePtr filePointer{fopen(fileName.c_str(), "r"), &fclose};
    size_t countOfHeader = 1; //we need to read only one header
    fread(&header_, sizeof(WavHeader), countOfHeader, filePointer.get());
    data.resize(header_.subchunk2Size);

    /* unsafe - using C pointers*/
    {
        uint8_t * raw_data = new uint8_t[header_.subchunk2Size];
        size_t countOfData = 1; //we need to read only one block of data
        fread(raw_data, header_.subchunk2Size, countOfData, filePointer.get());

        std::copy(raw_data, raw_data + header_.subchunk2Size, data.begin());

        delete[] raw_data;
    }
}

void WavCompressor::Compress(const std::string &fileName) {
    CompressData();

    FilePtr filePointer{fopen(fileName.c_str(), "w+"), &fclose};

    size_t countOfHeader = 1; //we need to write only one header
    fwrite(&header_, sizeof(WavHeader), countOfHeader, filePointer.get());

    size_t countOfData = 1; //we need to write only one block of data
    fwrite(data.data(), data.size(), countOfData, filePointer.get());
}

void WavCompressor::CompressData() {
    ComplexPoly initial(data.size());

    for(size_t i = 0; i < data.size(); ++i) {
        initial[i] = {static_cast<double>(data[i]), 0};
    }

    FourierHandler handler(initial);
    auto result = handler.GetFFT();

    size_t border =  static_cast<size_t>(compressionRate * static_cast<double>(result.size()));
    for(size_t i = result.size() - border; i < result.size(); ++i) {
        result[i] = {0, 0};
    }

    auto new_data = handler.GetReverseFFT(result);
    for(size_t i = 0; i < data.size(); ++i) {
        data[i] = new_data[i].real();
    }
}
