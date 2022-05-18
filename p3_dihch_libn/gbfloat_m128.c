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

int min(int a, int b) {
    return a < b ? a : b;
}

int max(int a, int b) {
    return a > b ? a : b;
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

#pragma omp parallel for schedule(dynamic) default(none) shared(gv, b, ext, a)
    for (int x = 0; x < a.dimX; x++) {
        for (int y = 0; y < a.dimY; y++) {

            float *pc = get_pixel(b, x, y);
            int deta = min(min(min(a.dimY - y - 1, y), min(a.dimX - x - 1, x)), gv.min_deta);

            __m128 sum = _mm_setzero_ps();

            int i = deta;
            int bound = gv.length - deta;

            for (; i < bound; ++i) {
                sum = _mm_fmadd_ps(_mm_load1_ps(gv.data + i), _mm_loadu_ps(get_pixel(a, x + i - ext, y)), sum);
            }
            pc[0] += sum[0] / gv.sum[ext - deta];
            pc[1] += sum[1] / gv.sum[ext - deta];
            pc[2] += sum[2] / gv.sum[ext - deta];
        }
    }
    return b;
}


Image gb_v(Image a, FVec gv) {
    Image b = img_sc(a);
    int ext = (int) gv.length / 2;

#pragma omp parallel for schedule(dynamic) default(none) shared(gv, b, ext, a)
    for (int y = 0; y < a.dimY; ++y) {
        for (int x = 0; x < a.dimX; ++x) {

            float *pc = get_pixel(b, x, y);
            int deta = min(min(min(a.dimY - y - 1, y), min(a.dimX - x - 1, x)), gv.min_deta);

            __m128 sum = _mm_setzero_ps();

            int i = deta;
            int bound = gv.length - deta;

            for (; i < bound; ++i) {
                sum = _mm_fmadd_ps(_mm_load1_ps(gv.data + i), _mm_loadu_ps(get_pixel(a, x, y + i - ext)), sum);
            }
            pc[0] = sum[0] / gv.sum[ext - deta];
            pc[1] = sum[1] / gv.sum[ext - deta];
            pc[2] = sum[2] / gv.sum[ext - deta];
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

