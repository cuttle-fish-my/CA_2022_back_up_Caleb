#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <algorithm>
#include <sys/time.h>
#include <time.h>
#include <immintrin.h>
#include <xmmintrin.h>
#include <omp.h>
//implement dynamic

#define STB_IMAGE_IMPLEMENTATION

#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image_write.h"

#define PI 3.14159


#define HORIZONTAL_BLOCK_SIZE 8

#define VERTICAL_BLOCK_SIZE 8

typedef struct FVec {
    int length;
    int min_length;
    int min_deta;
    float *data;
    float *sum;
} FVec;

typedef struct Image {
    int dimX, dimY, numChannels;
    float *data;
} Image;

inline int min(int a, int b) {
    return a < b ? a : b;
}

void normalize_FVec(FVec v) {
    int i, j;
    int ext = v.length / 2;
    v.sum[0] = v.data[ext];
    for (i = ext + 1, j = 1; i < v.length; i++, j++) {
        v.sum[j] = v.sum[j - 1] + v.data[i] * 2;
    }
}

float *get_pixel(Image img, int x, int y) {
    // the data array assigned as following:
    // |R,G,B|R,G,R|R,G,B|R,G,B|R,G,B|
    if (x < 0) {
        x = 0;
    }
    if (x >= img.dimX) {
        x = (int) (img.dimX - 1);
    }
    if (y < 0) {
        y = 0;
    }
    if (y >= img.dimY) {
        y = (int) (img.dimY - 1);
    }
    // return the address of a set of RGB value
    return img.data + img.numChannels * (y * img.dimX + x);
}

float gd(float a, float b, float x) {
    // @a is sigma:
    //      the standard deviation - influences how significantly the
    //      center pixel’s neighboring pixels affect the computations result.
    // @b is mu(center position)
    float c = (x - b) / a;
    return exp((-.5) * c * c) / (a * sqrt(2 * PI));
}

FVec make_gv(float a, float x0, float x1, int length, int min_length) {
    FVec v;
    v.length = length;
    v.min_length = min_length;
    if (v.min_length > v.length) {
        v.min_deta = 0;
    } else {
        v.min_deta = ((v.length - v.min_length) / 2);
    }
    v.data = malloc(length * sizeof(float));
    v.sum = malloc((length / 2 + 1) * sizeof(float));
    float step = (x1 - x0) / ((float) length);
    int offset = (int) length / 2;
    for (int i = 0; i < length; i++) {
        v.data[i] = gd(a, 0.0f, (float) (i - offset) * step);
    }
    normalize_FVec(v);
    return v;
}


Image img_sc(Image a) {
    Image b = a;
    b.data = malloc(b.dimX * b.dimY * b.numChannels * sizeof(float));
    return b;
}


Image gb_h(Image a, FVec gv) {
    Image b = img_sc(a);

    int ext = (int) gv.length / 2;
    int offset;
    int x, y, channel;
    float *pc;
    float sum[3];
    int i;
    __m256i index = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    __m256i zeros = _mm256_setzero_si256();
    __m256i dimYMinus1 = _mm256_set1_epi32(a.dimY - 1);
    __m256i dimY = _mm256_set1_epi32(a.dimY);
    __m256i dimXMinus1 = _mm256_set1_epi32(a.dimX - 1);
    __m256i dimX = _mm256_set1_epi32(a.dimX);
    __m256i numChannels = _mm256_set1_epi32(a.numChannels);

    int blockSizeV = VERTICAL_BLOCK_SIZE;
    int blockSizeH = HORIZONTAL_BLOCK_SIZE;

    int xb, yb;

#pragma omp parallel for default(none) schedule(dynamic) private(x, y, pc, i, sum, offset, channel, yb, xb) shared(gv, b, ext, a, index, zeros, dimYMinus1, dimX, dimXMinus1, dimY, numChannels, blockSizeH, blockSizeV)
    for (y = 0; y < a.dimY; y++) {
        for (x = 0; x < a.dimX; x++) {

            pc = get_pixel(b, x, y);
            int deta = min(min(a.dimY - y - 1, y), min(a.dimX - x - 1, x));
            deta = min(deta, gv.min_deta);
            __m256 gvSum = _mm256_set1_ps(gv.sum[ext - deta]);
            __m256i location[3];
            __m256 pixel_value_p[3];
            sum[0] = 0;
            sum[1] = 0;
            sum[2] = 0;
            __m256 sum_p[3] = {_mm256_setzero_ps(), _mm256_setzero_ps(), _mm256_setzero_ps()};

            for (i = deta; i < deta / 8 * 8 + 8; ++i) {
                offset = i - ext;
                sum[0] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[0];
                sum[1] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[1];
                sum[2] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[2];
            }

            __m256i y_p = _mm256_set1_epi32(y);
            y_p = _mm256_max_epi32(y_p, zeros);
            y_p = _mm256_min_epi32(y_p, dimYMinus1);
            y_p = _mm256_mullo_epi32(y_p, dimX);

            for (i = deta / 8 * 8 + 8; i < (gv.length - deta) / 8 * 8; i += 8) {
                __m256i x_p = _mm256_add_epi32(_mm256_set1_epi32(x + i - ext), index);

                x_p = _mm256_max_epi32(x_p, zeros);
                x_p = _mm256_min_epi32(x_p, dimXMinus1);
                __m256i index_p = _mm256_add_epi32(y_p, x_p);

                index_p = _mm256_mullo_epi32(index_p, numChannels);
                __m256 res_p = _mm256_div_ps(_mm256_loadu_ps(gv.data + i), gvSum);
                location[0] = _mm256_add_epi32(index_p, _mm256_set1_epi32(0));
                location[1] = _mm256_add_epi32(index_p, _mm256_set1_epi32(1));
                location[2] = _mm256_add_epi32(index_p, _mm256_set1_epi32(2));
                pixel_value_p[0] = _mm256_i32gather_ps(a.data, location[0], sizeof(float));
                pixel_value_p[1] = _mm256_i32gather_ps(a.data, location[1], sizeof(float));
                pixel_value_p[2] = _mm256_i32gather_ps(a.data, location[2], sizeof(float));
                sum_p[0] = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_p[0]), sum_p[0]);
                sum_p[1] = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_p[1]), sum_p[1]);
                sum_p[2] = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_p[2]), sum_p[2]);
            }

            for (i = (gv.length - deta) / 8 * 8; i < gv.length - deta; ++i) {
                offset = i - ext;
                sum[0] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[0];
                sum[1] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[1];
                sum[2] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[2];
            }
            for (channel = 0; channel < a.numChannels; channel++) {
                for (i = 0; i < 8; ++i) {
                    sum[channel] += sum_p[channel][i];
                }
            }
            pc[0] = sum[0];
            pc[1] = sum[1];
            pc[2] = sum[2];
        }
    }
    return b;
}


Image gb_v(Image a, FVec gv) {
    Image b = img_sc(a);

    int ext = (int) gv.length / 2;
    int offset;
    int x, y, channel;
    float *pc;
    float sum[3];
    int i;
    __m256i index = _mm256_set_epi32(7, 6, 5, 4, 3, 2, 1, 0);
    __m256i zeros = _mm256_setzero_si256();
    __m256i dimYMinus1 = _mm256_set1_epi32(a.dimY - 1);
    __m256i dimY = _mm256_set1_epi32(a.dimY);
    __m256i dimXMinus1 = _mm256_set1_epi32(a.dimX - 1);
    __m256i dimX = _mm256_set1_epi32(a.dimX);
    __m256i numChannels = _mm256_set1_epi32(a.numChannels);

    int blockSizeV = VERTICAL_BLOCK_SIZE;
    int blockSizeH = HORIZONTAL_BLOCK_SIZE;

    int xb, yb;

#pragma omp parallel for schedule(dynamic) default(none) private(x, y, pc, i, sum, offset, channel, xb, yb) shared(gv, b, ext, a, index, zeros, dimY, dimXMinus1, dimX, dimYMinus1, numChannels, blockSizeH, blockSizeV)
    for (x = 0; x < a.dimX; ++x) {
        for (y = 0; y < a.dimY; ++y) {
            pc = get_pixel(b, x, y);
            int deta = min(min(a.dimY - y - 1, y), min(a.dimX - x - 1, x));
            deta = min(deta, gv.min_deta);
            __m256 gvSum = _mm256_set1_ps(gv.sum[ext - deta]);
            __m256i location[3];
            __m256 pixel_value_p[3];
            sum[0] = 0;
            sum[1] = 0;
            sum[2] = 0;
            __m256 sum_p[3] = {_mm256_setzero_ps(), _mm256_setzero_ps(), _mm256_setzero_ps()};

            for (i = deta; i < deta / 8 * 8 + 8; ++i) {
                offset = i - ext;
                sum[0] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[0];
                sum[1] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[1];
                sum[2] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[2];
            }

            __m256i x_p = _mm256_set1_epi32(x);
            x_p = _mm256_max_epi32(x_p, zeros);
            x_p = _mm256_min_epi32(x_p, dimXMinus1);

            for (i = deta / 8 * 8 + 8; i < (gv.length - deta) / 8 * 8; i += 8) {

                __m256i y_p = _mm256_add_epi32(_mm256_set1_epi32(y + i - ext), index);
                y_p = _mm256_max_epi32(y_p, zeros);
                y_p = _mm256_min_epi32(y_p, dimYMinus1);
                y_p = _mm256_mullo_epi32(y_p, dimX);

                __m256i index_p = _mm256_add_epi32(y_p, x_p);
                index_p = _mm256_mullo_epi32(index_p, numChannels);

                __m256 res_p = _mm256_div_ps(_mm256_loadu_ps(gv.data + i), gvSum);

                location[0] = _mm256_add_epi32(index_p, _mm256_set1_epi32(0));
                location[1] = _mm256_add_epi32(index_p, _mm256_set1_epi32(1));
                location[2] = _mm256_add_epi32(index_p, _mm256_set1_epi32(2));

                pixel_value_p[0] = _mm256_i32gather_ps(a.data, location[0], sizeof(float));
                pixel_value_p[1] = _mm256_i32gather_ps(a.data, location[1], sizeof(float));
                pixel_value_p[2] = _mm256_i32gather_ps(a.data, location[2], sizeof(float));

                sum_p[0] = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_p[0]), sum_p[0]);
                sum_p[1] = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_p[1]), sum_p[1]);
                sum_p[2] = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_p[2]), sum_p[2]);
            }

            for (i = (gv.length - deta) / 8 * 8; i < gv.length - deta; ++i) {
                offset = i - ext;
                sum[0] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[0];
                sum[1] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[1];
                sum[2] += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[2];

            }
            for (channel = 0; channel < a.numChannels; channel++) {
                for (i = 0; i < 8; ++i) {
                    sum[channel] += sum_p[channel][i];
                }
            }
            pc[0] = sum[0];
            pc[1] = sum[1];
            pc[2] = sum[2];
        }
    }
    return b;
}

Image apply_gb(Image a, FVec gv) {
    Image b = gb_h(a, gv);
    Image c = gb_v(b, gv);
    free(b.data);
    return c;
}

int main(int argc, char **argv) {
    struct timeval start_time, stop_time, elapsed_time;
    gettimeofday(&start_time, NULL);
    if (argc < 6) {
        printf("Usage: ./gb.exe <inputjpg> <outputname> <float: a> <float: x0> <float: x1> <unsigned int: dim>\n");
        exit(0);
    }

    float a, x0, x1;
    int dim, min_dim;

    sscanf(argv[3], "%f", &a);
    sscanf(argv[4], "%f", &x0);
    sscanf(argv[5], "%f", &x1);
    sscanf(argv[6], "%u", &dim);
    sscanf(argv[7], "%u", &min_dim);

    FVec v = make_gv(a, x0, x1, dim, min_dim);
//    print_fvec(v);

    Image img;
    img.data = stbi_loadf(argv[1], (int *) &(img.dimX), (int *) &(img.dimY), (int *) &(img.numChannels), 0);
    Image imgOut;
    imgOut = apply_gb(img, v);
    stbi_write_jpg(argv[2], (int) imgOut.dimX, (int) imgOut.dimY, (int) imgOut.numChannels, imgOut.data, 90);
    gettimeofday(&stop_time, NULL);
    timersub(&stop_time, &start_time, &elapsed_time);
    printf("%f \n", (double) elapsed_time.tv_sec + (double) elapsed_time.tv_usec / 1000000.0);
    free(imgOut.data);
    free(v.data);
    free(v.sum);
    free(img.data);
    return 0;
}

