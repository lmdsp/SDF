/*
 Copyright (C) 2014 Mikko Mononen (memon@inside.org)
 Copyright (C) 2009-2012 Stefan Gustavson (stefan.gustavson@gmail.com)

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

 /// Fork by Lorcán Mc Donagh 2019 https://www.lmdsp.com

#pragma once

#ifndef NANO_SDF_H
#define NANO_SDF_H

namespace nanosdf {

/// Internal precision.
using precision_t = float;

/// Pixel data type.
using pixel_t = uint8_t;

/// Coordinate type.
using coord_t = int32_t;

/// Sweep-and-update Euclidean distance transform of an anti-aliased image for contour textures.
/// Based on edtaa3func.c by Stefan Gustavson.
///
/// White (255) pixels are treated as object pixels, zero pixels are treated as background.
/// An attempt is made to treat anti-aliased edges correctly. The input image must have
/// pixels in the range [0,255], and the anti-aliased image should be a box-filter
/// sampling of the ideal, crisp edge. If the anti-alias region is more than 1 pixel wide,
/// the result from this transform will be inaccurate.
/// Pixels at image border are not calculated and are set to 0.
///
/// The output distance field is encoded as bytes, where 0 = radius (outside) and 255 = -radius (inside).
/// Input and output can be the same buffer.
///   out - Output of the distance transform, one byte per pixel.
///   outstride - Bytes per row on output image. 
///   radius - The radius of the distance field narrow band in pixels.
///   img - Input image, one byte per pixel.
///   width - Width if the image. 
///   height - Height if the image. 
///   stride - Bytes per row on input image. 
int sdfBuildDistanceField(pixel_t* out, coord_t outstride, precision_t radius, const pixel_t* img, coord_t width, coord_t height, coord_t stride);

/// Same as distXform, but does not allocate any memory.
/// The 'temp' array should be enough to fit width * height * sizeof(float) * 3 bytes.
void sdfBuildDistanceFieldNoAlloc(pixel_t* out, coord_t outstride, precision_t radius, const pixel_t* img, coord_t width, coord_t height, coord_t stride, pixel_t* temp);

/// This function converts the anti-aliased image where each pixel represents coverage (box-filter
/// sampling of the ideal, crisp edge) to a distance field with narrow band radius of sqrt(2).
/// This is the fastest way to turn anti-aliased image to contour texture. This function is good
/// if you don't need the distance field for effects (i.e. fat outline or drop-shadow).
/// Input and output buffers must be different.
///   out - Output of the distance transform, one byte per pixel.
///   outstride - Bytes per row on output image. 
///   img - Input image, one byte per pixel.
///   width - Width if the image. 
///   height - Height if the image. 
///   stride - Bytes per row on input image. 
void sdfCoverageToDistanceField(pixel_t* out, coord_t outstride, const pixel_t* img, coord_t width, coord_t height, coord_t stride);

#endif // NANO_SDF_H

#ifdef SDF_IMPLEMENTATION

#include <cmath>
#include <cstdlib>

namespace constants {

/// Maximum number of distance transform passes
constexpr int max_passes = 10;	

/// Controls how much smaller the neighbor value must be to consider, too small slack increase iteration count.
constexpr precision_t black = 0.001f;

/// sqrt(2)
constexpr precision_t sqrt_2 = 1.4142136f;

/// 1 / sqrt(2)
constexpr precision_t sqrt_2_inv = 1.0f / sqrt_2;

/// Big value used to initialize the distance field.
constexpr precision_t big_value = 1e+37f;

}

namespace detail {

[[nodiscard]] constexpr float sdf__clamp01(float x) noexcept
{
	return (x < 0.0f) ? 0.0f : ((x > 1.0f) ? 1.0f : x);
}

}

void sdfCoverageToDistanceField(pixel_t* out, coord_t outstride, const pixel_t* img, coord_t width, coord_t height, coord_t stride)
{
	using constants::sqrt_2;

	// Zero out borders
	for (coord_t x = 0; x < width; ++x)
	{
		out[x] = 0;
	}

	for (coord_t y = 1; y < height; ++y)
	{
		out[y * outstride] = 0;
		out[width - 1 + y * outstride] = 0;
	}

	for (coord_t x = 0; x < width; ++x)
	{
		out[x + (height - 1) * outstride] = 0;
	}

	for (coord_t y = 1; y < height - 1; ++y)
	{
		for (coord_t x = 1; x < width-1; ++x) 
		{
			const coord_t k = x + y * stride;

			// Skip flat areas.
			if (img[k] == 255)
			{
				out[x + y * outstride] = 255;

				continue;
			}

			if (img[k] == 0)
			{
				// Special handling for cases where full opaque pixels are next to full transparent pixels.
				// See: https://github.com/memononen/SDF/issues/2
				const bool he = (img[k - 1] == 255) || (img[k + 1] == 255);
				const bool ve = (img[k - stride] == 255) || (img[k + stride] == 255);

				if (!he && !ve)
				{
					out[x + y * outstride] = 0;

					continue;
				}
			}

			const precision_t px_ms_m1 = static_cast<precision_t>(img[k - stride - 1]);
			const precision_t px_ms_p1 = static_cast<precision_t>(img[k - stride + 1]);
			const precision_t px_ps_m1 = static_cast<precision_t>(img[k + stride - 1]);
			const precision_t px_ps_p1 = static_cast<precision_t>(img[k + stride + 1]);

			precision_t gx = -px_ms_m1 - sqrt_2 * static_cast<precision_t>(img[k - 1]) - px_ps_m1 + px_ms_p1 + sqrt_2 * static_cast<precision_t>(img[k + 1]) + px_ps_p1;
			precision_t gy = -px_ms_m1 - sqrt_2 * static_cast<precision_t>(img[k - stride]) - px_ms_p1 + px_ps_m1 + sqrt_2 * static_cast<precision_t>(img[k + stride]) + px_ps_p1;

			const precision_t a = static_cast<precision_t>(img[k]) / 255.0f;

			gx = std::fabs(gx);
			gy = std::fabs(gy);

			precision_t d = 0.0f;

			if ((gx < 0.0001f) || (gy < 0.000f))
			{
				d = (0.5f - a) * sqrt_2;
			}
			else
			{
				const precision_t glen2 = gx * gx + gy * gy;
				const precision_t gnorm = 1.0f / std::sqrt(glen2);

				gx *= gnorm;
				gy *= gnorm;

				if (gx < gy) 
				{
					std::swap(gx, gy);
				}

				const precision_t a1 = 0.5f * gy / gx;

				if (a < a1)
				{ 
					// 0 <= a < a1
					d = 0.5f * (gx + gy) - std::sqrt(2.0f * gx * gy * a);
				}
				else if (a < (1.0 - a1))
				{ 
					// a1 <= a <= 1 - a1
					d = (0.5f - a) * gx;
				}
				else
				{ 
					// 1 - a1 < a <= 1
					d = -0.5f * (gx + gy) + std::sqrt(2.0f * gx * gy * (1.0f - a));
				}
			}

			d *= constants::sqrt_2_inv;

			out[x + y * outstride] = static_cast<pixel_t>(detail::sdf__clamp01(0.5f - d) * 255.0f);
		}
	}
}

namespace detail {

[[nodiscard]] precision_t sdf_edgedf(precision_t gx, precision_t gy, precision_t a) noexcept
{
	if ((gx == 0) || (gy == 0)) 
	{
		// Either A) gu or gv are zero, or B) both
		// Linear approximation is A) correct or B) a fair guess
		return 0.5f - a;
	}
	
	// Everything is symmetric w.r.t sign and transposition,
	// so move to first octant (gx >= 0, gy >= 0, gx >= gy) to
	// avoid handling all possible edge directions.
	gx = std::fabs(gx);
	gy = std::fabs(gy);

	if (gx < gy) 
	{
		std::swap(gx, gy);
	}

	const precision_t a1 = 0.5f * gy / gx;

	if (a < a1) 
	{ 
		// 0 <= a < a1
		return 0.5f * (gx + gy) - std::sqrt(2.0f * gx * gy * a);
	}
	 
	if (a < (1.0 - a1))
	{ 
		// a1 <= a <= 1 - a1
		return (0.5f - a) * gx;
	}

	// 1 - a1 < a <= 1
	return -0.5f * (gx + gy) + std::sqrt(2.0f * gx * gy * (1.0f - a));
}

}

struct Point 
{
	precision_t x;
	precision_t y;

	[[nodiscard]] precision_t distance_squared(Point other) const noexcept
	{
		const precision_t dx = other.x - x;
		const precision_t dy = other.y - y;

		return dx * dx + dy * dy;
	}
};

void sdfBuildDistanceFieldNoAlloc(pixel_t* out, coord_t outstride, float radius, const pixel_t* img, coord_t width, coord_t height, coord_t stride, pixel_t* temp)
{
	using constants::sqrt_2;

	const coord_t pixel_count = width * height;

	precision_t* const __restrict tdist = reinterpret_cast<precision_t*>(&temp[0]);
	Point* const __restrict tpt = reinterpret_cast<Point*>(&temp[pixel_count * sizeof(precision_t)]);
	
	// Initialize buffers
	for (coord_t i = 0; i < pixel_count; ++i)
	{
		tpt[i].x = 0;
		tpt[i].y = 0;
		tdist[i] = constants::big_value;
	}

	// Calculate position of the anti-aliased pixels and distance to the boundary of the shape.
	const coord_t height_m1 = height - 1;
	const coord_t width_m1 = width - 1;

	for (coord_t y = 1; y < height_m1; ++y)
	{
		for (coord_t x = 1; x < width_m1; x++) 
		{
			 const coord_t k = x + y * stride;

			// Skip flat areas.
			if (255 == img[k])
			{
				continue;
			}

			if (0 == img[k])
			{
				// Special handling for cases where full opaque pixels are next to full transparent pixels.
				// See: https://github.com/memononen/SDF/issues/2
				const bool he = (img[k - 1] == 255) || (img[k + 1] == 255);
				const bool ve = (img[k - stride] == 255) || (img[k + stride] == 255);

				if (!he && !ve)
				{
					continue;
				}
			}

			// Calculate gradient direction
			precision_t gx = -(float)img[k - stride - 1] - sqrt_2 * (float)img[k - 1] - (float)img[k + stride - 1] + (float)img[k - stride + 1] + sqrt_2 * (float)img[k + 1] + (float)img[k + stride + 1];
			precision_t gy = -(float)img[k - stride - 1] - sqrt_2 * (float)img[k - stride] - (float)img[k - stride + 1] + (float)img[k + stride - 1] + sqrt_2 * (float)img[k + stride] + (float)img[k + stride + 1];

			if ((std::fabs(gx) < 0.001f) && (std::fabs(gy) < 0.001f))
			{
				continue;
			}

			precision_t glen = gx * gx + gy * gy;
			
			if (glen > 0.0001f) 
			{
				glen = 1.0f / std::sqrt(glen);
				gx *= glen;
				gy *= glen;
			}

			// Find nearest point on contour.
			const coord_t tk = x + y * width;
			const precision_t d = detail::sdf_edgedf(gx, gy, static_cast<precision_t>(img[k]) / 255.0f);

			tpt[tk].x = x + gx * d;
			tpt[tk].y = y + gy * d;

			const Point c = { static_cast<precision_t>(x), static_cast<precision_t>(y) };

			tdist[tk] = c.distance_squared(tpt[tk]);
		}
	}

	using constants::black;

	// Calculate distance transform using sweep-and-update.
	for (int pass = 0; pass < constants::max_passes; ++pass)
	{
		int changed = 0;

		// Bottom-left to top-right.
		for (coord_t y = 1; y < height_m1; ++y) 
		{
			for (coord_t x = 1; x < width_m1; ++x) 
			{
				const coord_t k = x + y * width;
				const Point c = { static_cast<precision_t>(x), static_cast<precision_t>(y) };

				Point pt{};
				precision_t pd = tdist[k];
				bool ch = false;
			
				// (-1, -1)
				{
					const coord_t kn = k - 1 - width;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				// (0, -1)
				{
					const coord_t kn = k - width;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				// (1,-1)
				{
					const coord_t kn = k + 1 - width;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				// (-1, 0)
				{
					const coord_t kn = k - 1;

					if (tdist[kn] < tdist[k])
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				if (ch) 
				{
					tpt[k] = pt;
					tdist[k] = pd;

					++changed;
				}
			}
		}

		// Top-right to bottom-left.
		const coord_t height_m2 = height - 2;
		const coord_t width_m2 = width - 2;

		for (coord_t y = height_m2; y > 0 ; --y)
		{
			for (coord_t x = width_m2; x > 0; --x) 
			{
				const coord_t k = x + y * width;
				const Point c = { static_cast<precision_t>(x), static_cast<precision_t>(y) };

				Point pt{};
				bool ch = false;
				precision_t pd = tdist[k];

				// (1, 0)
				{
					const coord_t kn = k + 1;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				// (-1,1)
				{
					const coord_t kn = k - 1 + width;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				// (0, 1)
				{
					const coord_t kn = k + width;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				// (1, 1)
				{
					const coord_t kn = k + 1 + width;

					if (tdist[kn] < pd)
					{
						const precision_t d = c.distance_squared(tpt[kn]);

						if ((d + black) < pd)
						{
							pt = tpt[kn];
							pd = d;
							ch = true;
						}
					}
				}

				if (ch) 
				{
					tpt[k] = pt;
					tdist[k] = pd;
					++changed;
				}
			}
		}

		if (0 == changed)
		{
			break;
		}
	}

	// Map to good range.
	const precision_t scale = 1.0f / radius;

	for (coord_t y = 0; y < height; ++y)
	{
		for (coord_t x = 0; x < width; ++x)
		{
			precision_t d = std::sqrt(tdist[x + y * width]) * scale;

			if (img[x + y * stride] > 127)
			{
				d = -d;
			}

			out[x + y * outstride] = static_cast<pixel_t>(detail::sdf__clamp01(0.5f - d * 0.5f) * 255.0f);
		}
	}
}

int sdfBuildDistanceField(unsigned char* out, int outstride, float radius,
						  const unsigned char* img, int width, int height, int stride)
{
	unsigned char* temp = (unsigned char*)malloc(width*height*sizeof(float)*3);
	if (temp == NULL) return 0;
	sdfBuildDistanceFieldNoAlloc(out, outstride, radius, img, width, height, stride, temp);
	free(temp);
	return 1;
}

}

#endif // SDF_IMPLEMENTATION
