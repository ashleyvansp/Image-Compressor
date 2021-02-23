# Image Compression Program

Applies downscaling, quantization, and the discrete cosine transform to achieve lossy compression on still image bitmap files.
Supports three quality modes (low, medium, or high), allowing the user to specify the tradeoff between decompressed image quality and compressed file size.

## Usage
After compiling with Make, run the following to compress a bitmap image:

    ./uvg_compress < quality > input_image.bmp compressed_image.uvg

The following command can then decompress the file:

    ./uvg_decompress compressed_image.uvg decompressed_image.bmp