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

#define THREAD_NUMBER 20

#define HORIZONTAL_BLOCK_SIZE 20

#define VERTICAL_BLOCK_SIZE 20

typedef struct FVec {
    unsigned int length;
    unsigned int min_length;
    unsigned int min_deta;
    float *data;
    float *sum;
} FVec;

typedef struct Image {
    unsigned int dimX, dimY, numChannels;
    float *data;
} Image;

void normalize_FVec(FVec v) {
    unsigned int i, j;
    unsigned int ext = v.length / 2;
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

FVec make_gv(float a, float x0, float x1, unsigned int length, unsigned int min_length) {
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

//void print_fvec(FVec v) {
//    unsigned int i;
//    printf("\n");
//    for (i = 0; i < v.length; i++) {
//        printf("%f ", v.data[i]);
//    }
//    printf("\n");
//}

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
    float sum_1, sum_2, sum_3;
    int i;
    omp_set_num_threads(THREAD_NUMBER);

#pragma omp parallel for default(none) private(x, y, pc, i, sum_1, sum_2, sum_3, offset) shared(gv, channel, b, ext, a)
    for (y = 0; y < a.dimY; y++) {
        for (x = 0; x < a.dimX; x++) {
            pc = get_pixel(b, x, y);
            unsigned int deta = fmin(fmin(a.dimY - y - 1, y), fmin(a.dimX - x - 1, x));
            deta = fmin(deta, gv.min_deta);
//            for (channel = 0; channel < a.numChannels; channel++) {
            sum_1 = 0;
            sum_2 = 0;
            sum_3 = 0;
            __m256 sum_1_p = _mm256_setzero_ps();
            __m256 sum_2_p = _mm256_setzero_ps();
            __m256 sum_3_p = _mm256_setzero_ps();

            for (i = deta; i < deta / 8 * 8 + 8; ++i) {
                offset = i - ext;
                sum_1 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[0];
                sum_2 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[1];
                sum_3 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[2];
            }

            __m256i y_p = _mm256_set1_epi32(y);
            y_p = _mm256_max_epi32(y_p, _mm256_set1_epi32(0));
            y_p = _mm256_min_epi32(y_p, _mm256_set1_epi32(a.dimY - 1));
            y_p = _mm256_mullo_epi32(y_p, _mm256_set1_epi32(a.dimX));


            for (i = deta / 8 * 8 + 8; i < (gv.length - deta) / 8 * 8; i += 8) {
                offset = x + i - ext;
                __m256i x_p = _mm256_set_epi32(
                        offset + 7,
                        offset + 6,
                        offset + 5,
                        offset + 4,
                        offset + 3,
                        offset + 2,
                        offset + 1,
                        offset + 0
                );

                x_p = _mm256_max_epi32(x_p, _mm256_set1_epi32(0));
                x_p = _mm256_min_epi32(x_p, _mm256_set1_epi32(a.dimX - 1));
                __m256i index_p = _mm256_add_epi32(y_p, x_p);

                index_p = _mm256_mullo_epi32(index_p, _mm256_set1_epi32(a.numChannels));

                __m256 res_p = _mm256_div_ps(((__m256 *) gv.data)[i / 8], _mm256_set1_ps(gv.sum[ext - deta]));

                __m256i index_1_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(0));
                __m256i index_2_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(1));
                __m256i index_3_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(2));

                __m256 pixel_value_1_p = _mm256_i32gather_ps(a.data, index_1_p, sizeof(float));
                __m256 pixel_value_2_p = _mm256_i32gather_ps(a.data, index_2_p, sizeof(float));
                __m256 pixel_value_3_p = _mm256_i32gather_ps(a.data, index_3_p, sizeof(float));

                sum_1_p = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_1_p), sum_1_p);
                sum_2_p = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_2_p), sum_2_p);
                sum_3_p = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_3_p), sum_3_p);
            }

            for (i = (gv.length - deta) / 8 * 8; i < gv.length - deta; ++i) {
                offset = i - ext;
                sum_1 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[0];
                sum_2 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[1];
                sum_3 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x + offset, y)[2];

            }

            for (i = 0; i < 8; ++i) {
                sum_1 += sum_1_p[i];
                sum_2 += sum_2_p[i];
                sum_3 += sum_3_p[i];
            }
            pc[0] = sum_1;
            pc[1] = sum_2;
            pc[2] = sum_3;
//            }
        }
    }
    return b;
}


Image gb_v(Image a, FVec gv) {
    Image b = img_sc(a);

    int ext = (int) gv.length / 2;
    int offset;
    int x, y, channel, xb, yb;
    float *pc;
    float sum_1, sum_2, sum_3;
    int i;
//    int blockSize = VERTICAL_BLOCK_SIZE;
    omp_set_num_threads(THREAD_NUMBER);

#pragma omp parallel for default(none) private(x, y, xb, yb, pc, i, sum_1, sum_2, sum_3, offset) shared(gv, channel, b, ext, a)
    for (x = 0; x < a.dimX; ++x) {
        for (y = 0; y < a.dimY; ++y) {
//            for (channel = 0; channel < a.numChannels; ++channel) {
                pc = get_pixel(b, x, y);
                unsigned int deta = fmin(fmin(a.dimY - y - 1, y), fmin(a.dimX - x - 1, x));
                deta = fmin(deta, gv.min_deta);
                for (channel = 0; channel < a.numChannels; channel++) {
                    sum_1 = 0;
                    sum_2 = 0;
                    sum_3 = 0;
                    __m256 sum_1_p = _mm256_setzero_ps();
                    __m256 sum_2_p = _mm256_setzero_ps();
                    __m256 sum_3_p = _mm256_setzero_ps();

                    for (i = deta; i < deta / 8 * 8 + 8; ++i) {
                        offset = i - ext;
                        sum_1 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[0];
                        sum_2 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[1];
                        sum_3 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[2];
                    }

                    __m256i x_p = _mm256_set1_epi32(x);
                    x_p = _mm256_max_epi32(x_p, _mm256_set1_epi32(0));
                    x_p = _mm256_min_epi32(x_p, _mm256_set1_epi32(a.dimX - 1));

                    for (i = deta / 8 * 8 + 8; i < (gv.length - deta) / 8 * 8; i += 8) {
                        offset = y + i - ext;
                        __m256i y_p = _mm256_set_epi32(
                                offset + 7,
                                offset + 6,
                                offset + 5,
                                offset + 4,
                                offset + 3,
                                offset + 2,
                                offset + 1,
                                offset + 0
                        );

                        y_p = _mm256_max_epi32(y_p, _mm256_set1_epi32(0));
                        y_p = _mm256_min_epi32(y_p, _mm256_set1_epi32(a.dimY - 1));
                        y_p = _mm256_mullo_epi32(y_p, _mm256_set1_epi32(a.dimX));
                        __m256i index_p = _mm256_add_epi32(y_p, x_p);
                        index_p = _mm256_mullo_epi32(index_p, _mm256_set1_epi32(a.numChannels));
//                        index_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(channel));
                        __m256 res_p = _mm256_div_ps(((__m256 *) gv.data)[i / 8], _mm256_set1_ps(gv.sum[ext - deta]));

                        __m256i index_1_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(0));
                        __m256i index_2_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(1));
                        __m256i index_3_p = _mm256_add_epi32(index_p, _mm256_set1_epi32(2));

                        __m256 pixel_value_1_p = _mm256_i32gather_ps(a.data, index_1_p, sizeof(float));
                        __m256 pixel_value_2_p = _mm256_i32gather_ps(a.data, index_2_p, sizeof(float));
                        __m256 pixel_value_3_p = _mm256_i32gather_ps(a.data, index_3_p, sizeof(float));

                        sum_1_p = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_1_p), sum_1_p);
                        sum_2_p = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_2_p), sum_2_p);
                        sum_3_p = _mm256_add_ps(_mm256_mul_ps(res_p, pixel_value_3_p), sum_3_p);
                    }

                    for (i = (gv.length - deta) / 8 * 8; i < gv.length - deta; ++i) {
                        offset = i - ext;
//                        sum += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[channel];
                        sum_1 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[0];
                        sum_2 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[1];
                        sum_3 += gv.data[i] / gv.sum[ext - deta] * (float) get_pixel(a, x, y + offset)[2];
                    }

                    for (i = 0; i < 8; ++i) {
                        sum_1 += sum_1_p[i];
                        sum_2 += sum_2_p[i];
                        sum_3 += sum_3_p[i];
                    }
                    pc[0] = sum_1;
                    pc[1] = sum_2;
                    pc[2] = sum_3;
            }
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
    unsigned int dim, min_dim;

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