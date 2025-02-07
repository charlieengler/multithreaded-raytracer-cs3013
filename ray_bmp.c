#include <stdio.h>
#include <stdint.h>
#include <errno.h>

#include "ray_bmp.h"
	
#define bmp_file_header_size 14
#define bmp_info_header_size 40

int render_bmp(struct framebuffer_pt4 *fb, const char *output_filepath) {
	FILE *f;
	if ((f = fopen(output_filepath, "wb")) == NULL) {
		fprintf(stderr, "error opening BMP '%s' for writing: %d %s\n", output_filepath, errno, strerror(errno));
		exit(1);
	}

	const int filesize = bmp_file_header_size + 3 * fb->width * fb->height;

	uint8_t bmp_file_header[bmp_file_header_size] = {
		// 2	The header field used to identify the BMP and DIB file is 0x42 0x4D in hexadecimal, same as BM in ASCII.
		'B', 'M',
		// 4	The size of the BMP file in bytes
		filesize, filesize >> 8, filesize >> 16, filesize >> 24,
		// 2	Reserved; actual value depends on the application that creates the image, if created manually can be 0
		0, 0,
		// 2	Reserved; actual value depends on the application that creates the image, if created manually can be 0
		0, 0,
		// 4	The offset, i.e. starting address, of the byte where the bitmap image data (pixel array) can be found.
		bmp_file_header_size + bmp_info_header_size, 0, 0, 0,
	};
	uint8_t bmp_info_header[bmp_info_header_size] = {
		// 4	the size of this header, in bytes (40)
		40, 0, 0, 0,
		// 4	the bitmap width in pixels (signed integer)
		fb->width, fb->width >> 8, fb->width >> 16, fb->width >> 24,	
		// 4	the bitmap height in pixels (signed integer)
		fb->height, fb->height >> 8, fb->height >> 16, fb->height >> 24,	
		// 2	the number of color planes (must be 1)
		1, 0,
		// 2	the number of bits per pixel, which is the color depth of the image. Typical values are 1, 4, 8, 16, 24 and 32.
		24, 0,
		// 4	the compression method being used. See the next table for a list of possible values
		0, 0, 0, 0,	// 0 = no compression
				// 4	the image size. This is the size of the raw bitmap data; a dummy 0 can be given for BI_RGB bitmaps.
		0, 0, 0, 0,
		// 4	the horizontal resolution of the image. (pixel per metre, signed integer)
		0, 0, 0, 0,
		// 4	the vertical resolution of the image. (pixel per metre, signed integer)
		0, 0, 0, 0,
		// 4	the number of colors in the color palette, or 0 to default to 2n
		0, 0, 0, 0,
		// 4	the number of important colors used, or 0 when every color is important; generally ignored
		0, 0, 0, 0,
	};

	fwrite(bmp_file_header, sizeof(bmp_file_header), 1, f);
	fwrite(bmp_info_header, sizeof(bmp_info_header), 1, f);
	for (int y = fb->height - 1; y >= 0; y--) {
		int bytes_in_row = 0;
		for (int x = 0; x < fb->width; x++) {
			// bmp_source is source data that we convert to png
			const pt4 *p = framebuffer_pt4_get(fb, x, y);
			uint8_t px[3] = {
				color_double_to_u8(p->v[2]),
				color_double_to_u8(p->v[1]),
				color_double_to_u8(p->v[0]),
			};
			bytes_in_row += sizeof(px);
			fwrite(px, sizeof(px), 1, f);
		}
		while (bytes_in_row++ % 4) {
			uint8_t pad = 0;
			fwrite(&pad, sizeof(pad), 1, f);
		}
	}

	fclose(f);
	return 0;
}

