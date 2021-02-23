/* uvg_compress.cpp

   Starter code for Assignment 4 (in C++). This program
    - Reads an input image in BMP format
     (Using bitmap_image.hpp, originally from 
      http://partow.net/programming/bitmap/index.html)
    - Transforms the image from RGB to YCbCr (i.e. "YUV").
    - Downscales the Cb and Cr planes by a factor of two
      (producing the same resolution that would result
       from 4:2:0 subsampling, but using interpolation
       instead of ignoring some samples)
    - Writes each colour plane (Y, then Cb, then Cr)
      in 8 bits per sample to the output file.

   B. Bird - 07/01/2020
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <math.h>
#include <cstdint>
#include "output_stream.hpp"
#include "bitmap_image.hpp"
#include "uvg_common.hpp"

std::vector<std::vector<int>> luminance(8);
std::vector<std::vector<int>> chrominance(8);

//A simple downscaling algorithm using averaging.
std::vector<std::vector<unsigned char> > scale_down(std::vector<std::vector<unsigned char> > source_image, unsigned int source_width, unsigned int source_height, int factor){

    unsigned int scaled_height = (source_height+factor-1)/factor;
    unsigned int scaled_width = (source_width+factor-1)/factor;

    //Note that create_2d_vector automatically initializes the array to all-zero
    auto sums = create_2d_vector<unsigned int>(scaled_height,scaled_width);
    auto counts = create_2d_vector<unsigned int>(scaled_height,scaled_width);

    for(unsigned int y = 0; y < source_height; y++)
        for (unsigned int x = 0; x < source_width; x++){
            sums.at(y/factor).at(x/factor) += source_image.at(y).at(x);
            counts.at(y/factor).at(x/factor)++;
        }

    auto result = create_2d_vector<unsigned char>(scaled_height,scaled_width);
    for(unsigned int y = 0; y < scaled_height; y++)
        for (unsigned int x = 0; x < scaled_width; x++)
            result.at(y).at(x) = (unsigned char)((sums.at(y).at(x)+0.5)/counts.at(y).at(x));
    return result;
}

// Applies the 2-dimensional discrete cosine transform to the "matrix" @g
std::vector<std::vector<float>> dct(const std::vector<std::vector<unsigned char>> &g){
    std::vector<std::vector<float>> dct_vals = create_2d_vector<float>(8,8);
    float cx, cy, sum;
    for (int y = 0; y < 8; y++){
        for (unsigned int x = 0; x < 8; x++){
            cx = x == 0 ? 1. / sqrt(2) : 1.;
            cy = y == 0 ? 1. / sqrt(2) : 1.;
            sum = 0;
            for (int j = 0; j < 8; j++){
                for (int i = 0; i < 8; i++){
                    sum += cx * cy * g.at(j).at(i) * cos(M_PI * x * (2.0 * i + 1)/(2*8)) * cos(M_PI * y * (2.0 * j + 1)/(2*8));
                }
            }
            dct_vals.at(y).at(x) = 1.0/4 * sum;
        }    
    }

    return dct_vals;
}

// Quantizes the matrix @dct_vals by dividing elementwise by each value in @quant_vals
std::vector<std::vector<int>> quantize(const std::vector<std::vector<float>> &dct_vals, const std::vector<std::vector<int>> &quant_vals, float quality){
    std::vector<std::vector<int>> q = create_2d_vector<int>(8,8);
    for (int j = 0; j < 8; j++){
        for (int i = 0; i < 8; i++){
            q.at(j).at(i) = std::round(dct_vals.at(j).at(i) / (quant_vals.at(j).at(i) / quality));
        }
    }
    return q;
}

// 10 -> positive
// 11 -> negative
// Followed by @difference 1's and a 0 to indicate a stop
void delta(OutputBitStream &stream, int difference){
    if (difference != 0){
        stream.push_bits(difference > 0 ? 1 : 3, 2);
        for (int i = 0; i < abs(difference); i++){
            stream.push_bit(1);
        }
    }
    stream.push_bit(0);
}

void write_codes(OutputBitStream &output_stream, const std::vector<std::vector<int>> &q){
    int prev_val = q.at(0).at(0);
    output_stream.push_u16(prev_val);
    for (int j = 0; j < 8; j++){
        for (int i = 0; i < 8; i++){
            if (j != 0 || i != 0){
                int difference = q.at(j).at(i) - prev_val;
                delta(output_stream, difference);
                prev_val = q.at(j).at(i);
            }
        }
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

    // TODO - encode both the DC and first AC val not in unary

    init();
    if (argc < 4){
        std::cerr << "Usage: " << argv[0] << " <low/medium/high> <input BMP> <output file>" << std::endl;
        return 1;
    }
    std::string quality{argv[1]};
    std::string input_filename {argv[2]};
    std::string output_filename {argv[3]};
    float quality_val = quality == "low" ? 0.2 : (quality == "medium" ? 1 : 2);

    bitmap_image input_image {input_filename};

    unsigned int height = input_image.height();
    unsigned int padded_height = height + (8 - (height % 8));
    unsigned int width = input_image.width();
    unsigned int padded_width = width + (8 - (width % 8));
    //Read the entire image into a 2d array of PixelRGB objects 
    //(Notice that height is the outer dimension, so the pixel at coordinates (x,y) 
    // must be accessed as imageRGB.at(y).at(x)).
    std::vector<std::vector<PixelYCbCr>> imageYCbCr = create_2d_vector<PixelYCbCr>(padded_height, padded_width);

    for(unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < width; x++){
            auto [r,g,b] = input_image.get_pixel(x,y);
            PixelRGB rgb_pixel {r,g,b};
            imageYCbCr.at(y).at(x) = rgb_pixel.to_ycbcr();            
        }
    }

    // Pad 
    unsigned int right_padding = 8 - (width % 8);
    unsigned int lower_padding = 8 - (height % 8);
    for (unsigned int y = 0; y < height; y++){
        for (unsigned int x = 0; x < right_padding; x++){
            imageYCbCr.at(y).at(x + width) = imageYCbCr.at(y).at(width - 1);
        }
    }

    for (unsigned int y = 0; y < lower_padding; y++){
        for (unsigned int x = 0; x < width; x++){
            imageYCbCr.at(y + height).at(x) = imageYCbCr.at(height - 1).at(x);
        }
    }

    std::ofstream output_file{output_filename,std::ios::binary};
    OutputBitStream output_stream {output_file};

    //Placeholder: Use a simple bitstream containing the height/width (in 32 bits each)
    //followed by the entire set of values in each colour plane (in row major order).

    output_stream.push_u32(height);
    output_stream.push_u32(width);
    output_stream.push_bits(quality == "low" ? 0 : (quality == "medium" ? 1 : 2), 2);

    // //Write the Y values 
    // for(unsigned int y = 0; y < height; y++)
    //     for (unsigned int x = 0; x < width; x++)
    //         output_stream.push_byte(imageYCbCr.at(y).at(x).Y);

    //Extract the Cb plane into its own array 
    auto Cb = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Cb.at(y).at(x) = imageYCbCr.at(y).at(x).Cb;
    auto Cb_scaled = scale_down(Cb,width,height,2);

    //Extract the Cr plane into its own array 
    auto Cr = create_2d_vector<unsigned char>(height,width);
    for(unsigned int y = 0; y < height; y++)
        for (unsigned int x = 0; x < width; x++)
            Cr.at(y).at(x) = imageYCbCr.at(y).at(x).Cr;
    auto Cr_scaled = scale_down(Cr,width,height,2);

    // Pad the scaled values
    unsigned int scaled_height = (height + 1)/2;
    unsigned int scaled_width = (width + 1)/2;
    right_padding = 8 - (scaled_width % 8);
    lower_padding = 8 - (scaled_height % 8);
    unsigned int padded_scaled_height = scaled_height + lower_padding;
    unsigned int padded_scaled_width = scaled_width + right_padding;

    Cb_scaled.resize(padded_scaled_height);
    Cr_scaled.resize(padded_scaled_height);
    for (unsigned int i = 0; i < padded_scaled_height; i++){
        Cb_scaled.at(i).resize(padded_scaled_width);
        Cr_scaled.at(i).resize(padded_scaled_width);
    }

    for (unsigned int y = 0; y < scaled_height; y++){
        for (unsigned int x = 0; x < right_padding; x++){
            Cb_scaled.at(y).at(x + scaled_width) = Cb_scaled.at(y).at(scaled_width - 1);
            Cr_scaled.at(y).at(x + scaled_width) = Cr_scaled.at(y).at(scaled_width - 1);
        }
    }
    for (unsigned int y = 0; y < lower_padding; y++){
        for (unsigned int x = 0; x < scaled_width; x++){
            Cb_scaled.at(y + scaled_height).at(x) = Cb_scaled.at(scaled_height - 1).at(x);
            Cr_scaled.at(y + scaled_height).at(x) = Cr_scaled.at(scaled_height - 1).at(x);
        }
    }

    // Use 8 by 8 blocks
    // If file dimensions are not divisible by 8 then must pad
    std::vector<std::vector<unsigned char>> y_values = create_2d_vector<unsigned char>(8,8);
    std::vector<std::vector<unsigned char>> cb_values = create_2d_vector<unsigned char>(8,8);
    std::vector<std::vector<unsigned char>> cr_values = create_2d_vector<unsigned char>(8,8);
    
    // Y values first
    unsigned int x = 0;
    unsigned int y = 0;
    while (y <= padded_height - 8){
        while (x <= padded_width - 8){
            for (unsigned int j = y; j < y + 8; j++){
                for (unsigned int i = x; i < x + 8; i++){
                    y_values.at(j % 8).at(i % 8) = imageYCbCr.at(j).at(i).Y;
                }
            }
            std::vector<std::vector<float>> y_dct =  dct(y_values);
            std::vector<std::vector<int>> y_quant = quantize(y_dct, luminance, quality_val);

            write_codes(output_stream, y_quant);
            x += 8;
        }
        y += 8;
        x = 0;
    }

    // Then colour planes
    x = 0;
    y = 0;
    while (y <= padded_scaled_height - 8){
        while (x <= padded_scaled_width - 8){
            for (unsigned int j = y; j < y + 8; j++){
                for (unsigned int i = x; i < x + 8; i++){
                    cb_values.at(j % 8).at(i % 8) = Cb_scaled.at(j).at(i);
                }
            }
            std::vector<std::vector<float>> cb_dct =  dct(cb_values);
            std::vector<std::vector<int>> cb_quant = quantize(cb_dct, chrominance, quality_val);

            write_codes(output_stream, cb_quant);

            x += 8;
        }
        y += 8;
        x = 0;
    }
    x = 0;
    y = 0;
    while (y <= padded_scaled_height - 8){
        while (x <= padded_scaled_width - 8){
            for (unsigned int j = y; j < y + 8; j++){
                for (unsigned int i = x; i < x + 8; i++){
                    cr_values.at(j % 8).at(i % 8) = Cr_scaled.at(j).at(i);
                }
            }
            std::vector<std::vector<float>> cr_dct =  dct(cr_values);
            std::vector<std::vector<int>> cr_quant = quantize(cr_dct, chrominance, quality_val);

            write_codes(output_stream, cr_quant);

            x += 8;
        }
        y += 8;
        x = 0;
    }
    output_stream.flush_to_byte();
    output_file.close();

    return 0;
}