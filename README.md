# SDF

Sweep-and-update Euclidean distance transform of an anti-aliased image for contour texturing.

The code is based on [edtaa3func.c](http://contourtextures.wikidot.com/) by Stefan Gustavson and improves the original in terms of memory usage and execution time.

The algorithms first traverse the image and uses gradient direction and the edge function from edtaa3 to find an approximated point on the contour of the input image.
After this pass the distance at the edge pixels are known, and the code proceeds to update the rest of the distance field using sweep-and-update until the distance field convergences (or max passes run).

The additional memory required by the code is 3 floats (`px`, `py`, `distance`) per pixel,
compared to 5 doubles (or floats) of edtaa3.
The sweep-and-update is done using squared distances and contour points calculated only in the first pass.
This is greatly reduces the amount of computation (especially square roots) per distance update.

The code produces comparable, but probably not as accurate distance fields as the original code.

The code is intended to be used to calculate distance fields for [contour texturing](http://contourtextures.wikidot.com/).

## Usage

```cpp
int sdfBuildDistanceField(unsigned char* out, int outstride, float radius,
						  const unsigned char* img, int width, int height, int stride);
```
The output distance field is encoded as bytes, where 0 = radius (outside) and 255 = -radius (inside). Input and output can be the same buffer.
* _out_ - Output of the distance transform, one byte per pixel.
* _outstride_ - Bytes per row on output image. 
* _radius_ - The radius of the distance field narrow band in pixels.
* _img_ - Input image, one byte per pixel.
* _width_ - Width if the image. 
* _height_ - Height if the image. 
* _stride_ - Bytes per row on input image.

White (255) pixels are treated as object pixels, zero pixels are treated as background. An attempt is made to treat anti-aliased edges correctly.
The input image must have pixels in the range `[0,255]`,
and the anti-aliased image should be a box-filter sampling of the ideal, crisp edge. If the antialias region is more than 1 pixel wide, the result from this transform will be inaccurate. Pixels at image border are not calculated and are set to 0.
(Explanation borrowed from the original eedtaa3func.c)

```cpp
void sdfBuildDistanceFieldNoAlloc(unsigned char* out, int outstride, float radius,
								  const unsigned char* img, int width, int height, int stride,
								  unsigned char* temp);
```
Same as distXform, but does not allocate any memory. The `temp` array should be enough to fit `width * height * sizeof(float) * 3` bytes.

```cpp
void sdfCoverageToDistanceField(unsigned char* out, int outstride,
								const unsigned char* img, int width, int height, int stride);
```
The output distance field is encoded as bytes, where `0 = sqrt(2)` (outside) and `255 = -sqrt(2)` (inside).
Input and output must be different buffers.
* _out_ - Output of the distance transform, one byte per pixel.
* _outstride_ - Bytes per row on output image. 
* _radius_ - The radius of the distance field narrow band in pixels.
* _img_ - Input image, one byte per pixel.
* _width_ - Width if the image. 
* _height_ - Height if the image. 
* _stride_ - Bytes per row on input image.

This function converts the anti-aliased image where each pixel represents coverage (box-filter sampling of the ideal, crisp edge) to a distance field with narrow band radius of `sqrt(2)`.
This is the fastest way (often up to 10x faster than `sdfBuildDistanceField`) to turn an image into contour atexture.
This function is good if you don't need the distance field for effects (i.e. fat outline or drop-shadow).

The code is single header file only.
Use following code once in your project to compile the implementation.

```cpp
#define SDF_IMPLEMENTATION
#include "sdf.h"
```

## License

MIT License

## Fork

This is a C++17 fork of the orignal [C version](https://github.com/memononen/SDF).
