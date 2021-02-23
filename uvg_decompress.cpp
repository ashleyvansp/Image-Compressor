/* uvg_decompress.cpp

   Starter code for Assignment 4 (in C++). This program
    - Reads a height/width value from the input file
    - Reads YCbCr data from the file, with the Y plane
      in full w x h resolution and the other two planes
      in half resolution.
    - Upscales the Cb and Cr planes to full resolution and
      transforms them to RGB.
    - Writes the result as a BMP image
     (Using bitmap_image.hpp, originally from 
      http://partow.net/programming/bitmap/index.html)

   B. Bird - 07/01/2020
*/

#include <iostream>
#include <fstream>
#include <array>
#include <string>
#include <cassert>
#include <cstdint>
#include "input_stream.hpp"
#include "bitmap_image.hpp"
#include "uvg_common.hpp"

std::vector<std::vector<int>> luminance(8);
std::vector<std::vector<int>> chrominance(8);

int read_delta(InputBitStream &stream){
    int start = stream.read_bit();
    if (start == 0){
        return 0;
    }else{
        int sign = stream.read_bit() == 0 ? 1 : -1;
        int value = 0;
        while (stream.read_bit() == 1){
            value++;
        }
        return value * sign;
    }
}

std::vector<std::vector<int>> get_quantized(InputBitStream &stream, int width, int height){
    auto quantized = create_2d_vector<int>(height, width);
    int y = 0;
    int x = 0;
    while (y <= height - 8){
        while (x <= width - 8){
            quantized.at(y).at(x) = stream.read_u16();
            int prev = quantized.at(y).at(x);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    if (j != 0 || i != 0){
                        quantized.at(y + j).at(x + i) = prev + read_delta(stream);
                        prev = quantized.at(y + j).at(x + i);
                    }
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
    return quantized;
}

std::vector<std::vector<int>> inv_dct(const std::vector<std::vector<int>> &q){
    std::vector<std::vector<int>> dct_vals = create_2d_vector<int>(8,8);
    float ci, cj, sum;
    for (int y = 0; y < 8; y++){
        for (int x = 0; x < 8; x++){
            sum = 0;
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    ci = i == 0 ? 1. / sqrt(2) : 1.;
                    cj = j == 0 ? 1. / sqrt(2) : 1.;
                    sum += ci * cj * q.at(j).at(i) * cos(M_PI * i * (2.0 * x + 1)/(2*8)) * cos(M_PI * j * (2.0 * y + 1)/(2*8));
                }
            }
            dct_vals.at(y).at(x) = round(1.0/4 * sum);
            if (dct_vals.at(y).at(x) > 255)
                dct_vals.at(y).at(x) = 255;
            if (dct_vals.at(y).at(x) < 0)
                dct_vals.at(y).at(x) = 0;
        }    
    }

    return dct_vals;
}

void inv_quantized(std::vector<std::vector<int>> &vals, const std::vector<std::vector<int>> &factors, unsigned int width, unsigned int height, float quality){
    unsigned int x = 0;
    unsigned int y = 0;
    while (y <= height - 8){
        while (x <= width - 8){
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    vals.at(y + j).at(x + i) *= (factors.at(j).at(i) / quality);
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
}

void init(){
    luminance = create_2d_vector<int>(8,8);
    chrominance = create_2d_vector<int>(8,8);

    luminance = {
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55}, 
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };
    
    for (int i = 0; i < 8; i++){
        for (int j = 0; j < 8; j++){
            chrominance.at(i).at(j) = 99;
        }
    }
    chrominance.at(0).at(0) = 17;
    chrominance.at(0).at(1) = 18;
    chrominance.at(0).at(2) = 24;
    chrominance.at(0).at(3) = 47;
    chrominance.at(1).at(0) = 18;
    chrominance.at(1).at(1) = 21;
    chrominance.at(1).at(2) = 26;
    chrominance.at(1).at(3) = 66;
    chrominance.at(2).at(0) = 24;
    chrominance.at(2).at(1) = 26;
    chrominance.at(2).at(2) = 56;
    chrominance.at(3).at(0) = 47;
    chrominance.at(3).at(1) = 66;
}

int main(int argc, char** argv){
    init();
    if (argc < 3){
        std::cerr << "Usage: " << argv[0] << " <input file> <output BMP>" << std::endl;
        return 1;
    }
    std::string input_filename {argv[1]};
    std::string output_filename {argv[2]};

    std::ifstream input_file{input_filename,std::ios::binary};
    InputBitStream input_stream {input_file};

    unsigned int height = input_stream.read_u32();
    unsigned int width = input_stream.read_u32();
    unsigned int quality_bits = input_stream.read_bits(2);
    float quality = quality_bits == 0 ? 0.2 : (quality_bits == 1 ? 1 : 2);

    unsigned int padded_height = height + (8 - (height % 8));
    unsigned int padded_width = width + (8 - (width % 8));

    unsigned int scaled_height = (height+1)/2;
    unsigned int scaled_width = (width+1)/2;
    unsigned int padded_scaled_height = scaled_height + 8 - (scaled_height % 8);
    unsigned int padded_scaled_width = scaled_width + 8 - (scaled_width % 8);

    // Y values first
    auto y_quant = get_quantized(input_stream, padded_width, padded_height);
    inv_quantized(y_quant, luminance, padded_width, padded_height, quality);

    // Then scale colour planes
    auto cb_quant = get_quantized(input_stream, padded_scaled_width, padded_scaled_height);
    inv_quantized(cb_quant, chrominance, padded_scaled_width, padded_scaled_height, quality);

    auto cr_quant = get_quantized(input_stream, padded_scaled_width, padded_scaled_height);
    inv_quantized(cr_quant, chrominance, padded_scaled_width, padded_scaled_height, quality);
    
    // Invert the DCT for the Y vals
    unsigned int x = 0;
    unsigned int y = 0;
    while (y <= padded_height - 8){
        while (x <= padded_width - 8){
            std::vector<std::vector<int>> vals = create_2d_vector<int>(8,8);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    vals.at(j).at(i) = y_quant.at(y + j).at(x + i);
                }
            }
            std::vector<std::vector<int>> inverted = inv_dct(vals);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    y_quant.at(y + j).at(x + i) = inverted.at(j).at(i);
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }

    // Invert the DCT for the colour planes
    x = 0;
    y = 0;
    while (y <= padded_scaled_height - 8){
        while (x <= padded_scaled_width - 8){
            std::vector<std::vector<int>> cb_vals = create_2d_vector<int>(8,8);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    cb_vals.at(j).at(i) = cb_quant.at(y + j).at(x + i);
                }
            }
            std::vector<std::vector<int>> cb_inverted = inv_dct(cb_vals);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    cb_quant.at(y + j).at(x + i) = cb_inverted.at(j).at(i);
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }

    x = 0;
    y = 0;
    while (y <= padded_scaled_height - 8){
        while (x <= padded_scaled_width - 8){
            std::vector<std::vector<int>> cr_vals = create_2d_vector<int>(8,8);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    cr_vals.at(j).at(i) = cr_quant.at(y + j).at(x + i);
                }
            }
            std::vector<std::vector<int>> cr_inverted = inv_dct(cr_vals);
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    cr_quant.at(y + j).at(x + i) = cr_inverted.at(j).at(i);
                }
            }
            x += 8;
        }
        y += 8;
        x = 0;
    }
    
    auto Y = create_2d_vector<unsigned char>(height,width);
    auto Cb_scaled = create_2d_vector<unsigned char>((height+1)/2,(width+1)/2);
    auto Cr_scaled = create_2d_vector<unsigned char>((height+1)/2,(width+1)/2);

    for (unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Y.at(y).at(x) = y_quant.at(y).at(x);

    for (unsigned int y = 0; y < (height+1)/2; y++)
        for (unsigned int x = 0; x < (width+1)/2; x++)
            Cb_scaled.at(y).at(x) = cb_quant.at(y).at(x);

    for (unsigned int y = 0; y < (height+1)/2; y++){
        for (unsigned int x = 0; x < (width+1)/2; x++){
            Cr_scaled.at(y).at(x) = cr_quant.at(y).at(x);
        }
    }
    
    auto imageYCbCr = create_2d_vector<PixelYCbCr>(height,width);
    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            imageYCbCr.at(y).at(x) = {
                Y.at(y).at(x),
                Cb_scaled.at(y/2).at(x/2),
                Cr_scaled.at(y/2).at(x/2)
            };
        }
    }
    
    input_stream.flush_to_byte();
    input_file.close();

    bitmap_image output_image {width,height};

    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            auto pixel_rgb = imageYCbCr.at(y).at(x).to_rgb();
            auto [r,g,b] = pixel_rgb;
            output_image.set_pixel(x,y,r,g,b);
        }
    }
    
    output_image.save_image(output_filename);

    return 0;
}