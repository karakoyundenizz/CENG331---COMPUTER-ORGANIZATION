/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"

/*
 * Please fill in the following team struct
 */
team_t team = {
    "biltes", /* Team name */

    "Deniz Karakoyun",             /* First member full name */
    "karakoyun.deniz@metu.edu.tr", /* First member email address */

    "", /* Second member full name (leave blank if none) */
    ""  /* Second member email addr (leave blank if none) */
};

/***************
 * BILATERAL FILTER KERNEL
 ***************/

/******************************************************
 * Your different versions of the bilateral filter kernel go here
 ******************************************************/

/* Global parameters for bilateral filter */
static int radius = 2;        /* Filter radius (window size is 2*radius+1) */
static double sigma_s = 2.0;  /* Spatial parameter */
static double sigma_r = 50.0; /* Range parameter */

/*
 * naive_bilateral - The naive baseline version of bilateral filter
 */
char naive_bilateral_descr[] = "naive_bilateral: Naive baseline implementation";
void naive_bilateral(int dim, pixel *src, pixel *dst)
{
    int i, j, u, v;
    int src_i, src_j;
    double sum_r, sum_g, sum_b;
    double weight_sum;
    double spatial_weight, range_weight, weight;
    double dist_sq, color_diff_sq;
    double two_sigma_s_sq = 2.0 * sigma_s * sigma_s;
    double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            sum_r = sum_g = sum_b = 0.0;
            weight_sum = 0.0;

            for (u = -radius; u <= radius; u++)
            {
                for (v = -radius; v <= radius; v++)
                {
                    src_i = i + u;
                    src_j = j + v;

                    /* Mirror padding at boundaries */
                    if (src_i < 0)
                        src_i = -src_i;
                    if (src_i >= dim)
                        src_i = 2 * dim - 1 - src_i;
                    if (src_j < 0)
                        src_j = -src_j;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    /* Spatial weight */
                    dist_sq = u * u + v * v;
                    spatial_weight = exp(-dist_sq / two_sigma_s_sq);

                    /* Range weight (using L2 norm of RGB difference) */
                    color_diff_sq = (src[RIDX(i, j, dim)].red -
                                     src[RIDX(src_i, src_j, dim)].red) *
                                        (src[RIDX(i, j, dim)].red -
                                         src[RIDX(src_i, src_j, dim)].red) +
                                    (src[RIDX(i, j, dim)].green -
                                     src[RIDX(src_i, src_j, dim)].green) *
                                        (src[RIDX(i, j, dim)].green -
                                         src[RIDX(src_i, src_j, dim)].green) +
                                    (src[RIDX(i, j, dim)].blue -
                                     src[RIDX(src_i, src_j, dim)].blue) *
                                        (src[RIDX(i, j, dim)].blue -
                                         src[RIDX(src_i, src_j, dim)].blue);
                    range_weight = exp(-color_diff_sq / two_sigma_r_sq);

                    weight = spatial_weight * range_weight;
                    weight_sum += weight;

                    sum_r += src[RIDX(src_i, src_j, dim)].red * weight;
                    sum_g += src[RIDX(src_i, src_j, dim)].green * weight;
                    sum_b += src[RIDX(src_i, src_j, dim)].blue * weight;
                }
            }

            dst[RIDX(i, j, dim)].red = (unsigned short)(sum_r / weight_sum);
            dst[RIDX(i, j, dim)].green = (unsigned short)(sum_g / weight_sum);
            dst[RIDX(i, j, dim)].blue = (unsigned short)(sum_b / weight_sum);
        }
    }
}

/*
 * bilateral - Your current working version of bilateral filter
 * IMPORTANT: This is the version you will be graded on
 */
char bilateral_descr[] = "bilateral: Current working version";
void bilateral(int dim, pixel *src, pixel *dst)
{
    int i, j, u, v;
    int src_i, src_j;
    double sum_r, sum_g, sum_b;
    double weight_sum;
    double spatial_weight, range_weight, weight;
    double dist_sq, color_diff_sq;
    double two_sigma_s_sq = 2.0 * sigma_s * sigma_s;
    double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            sum_r = sum_g = sum_b = 0.0;
            weight_sum = 0.0;

            for (u = -radius; u <= radius; u++)
            {
                // moved i+u here since it doesnt need to be executed for every v it is independent of v
                // after that ı decided to compare with 0 immediatly so that CC are already set no need for cmpq 0 etc
                src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                // part of dist_sq calculation moved here since again it is independent of v
                dist_sq = u * u;
                for (v = -radius; v <= radius; v++)
                {

                    /* Spatial weight moved here maybe it will be good since we lastly used dist_sq and reaching will be easy*/
                    dist_sq += v * v;
                    spatial_weight = exp(-dist_sq / two_sigma_s_sq);
                    src_j = j + v;

                    /* Mirror padding at boundaries */
                    if (src_j < 0)
                        src_j = -src_j;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    /* Range weight (using L2 norm of RGB difference) */
                    color_diff_sq = (src[RIDX(i, j, dim)].red -
                                     src[RIDX(src_i, src_j, dim)].red) *
                                        (src[RIDX(i, j, dim)].red -
                                         src[RIDX(src_i, src_j, dim)].red) +
                                    (src[RIDX(i, j, dim)].green -
                                     src[RIDX(src_i, src_j, dim)].green) *
                                        (src[RIDX(i, j, dim)].green -
                                         src[RIDX(src_i, src_j, dim)].green) +
                                    (src[RIDX(i, j, dim)].blue -
                                     src[RIDX(src_i, src_j, dim)].blue) *
                                        (src[RIDX(i, j, dim)].blue -
                                         src[RIDX(src_i, src_j, dim)].blue);
                    range_weight = exp(-color_diff_sq / two_sigma_r_sq);

                    weight = spatial_weight * range_weight;
                    weight_sum += weight;

                    sum_r += src[RIDX(src_i, src_j, dim)].red * weight;
                    sum_g += src[RIDX(src_i, src_j, dim)].green * weight;
                    sum_b += src[RIDX(src_i, src_j, dim)].blue * weight;
                }
            }

            dst[RIDX(i, j, dim)].red = (unsigned short)(sum_r / weight_sum);
            dst[RIDX(i, j, dim)].green = (unsigned short)(sum_g / weight_sum);
            dst[RIDX(i, j, dim)].blue = (unsigned short)(sum_b / weight_sum);
        }
    }
}

/*********************************************************************
 * register_bilateral_functions - Register all of your different versions
 *     of the bilateral filter kernel with the driver by calling the
 *     add_bilateral_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 *********************************************************************/

void register_bilateral_functions()
{
    add_bilateral_function(&naive_bilateral, naive_bilateral_descr);
    add_bilateral_function(&bilateral, bilateral_descr);
    /* ... Register additional test functions here */
}

/***************
 * BOX BLUR KERNEL
 **************/

/******************************************************
 * Your different versions of the box blur kernel go here
 ******************************************************/

/* Global parameter for box blur */
static int blur_radius = 2; /* Filter radius (window size is 2*radius+1) */

/*
 * naive_blur - The naive baseline version of box blur
 */
char naive_blur_descr[] = "naive_blur: Naive baseline implementation";
void naive_blur(int dim, pixel *src, pixel *dst)
{
    int i, j, u, v;
    int src_i, src_j;
    int sum_r, sum_g, sum_b;
    int window_size = 2 * blur_radius + 1;
    int window_area = window_size * window_size;

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            sum_r = sum_g = sum_b = 0;

            for (u = -blur_radius; u <= blur_radius; u++)
            {
                for (v = -blur_radius; v <= blur_radius; v++)
                {
                    src_i = i + u;
                    src_j = j + v;

                    /* Mirror padding at boundaries */
                    if (src_i < 0)
                        src_i = -src_i;
                    if (src_i >= dim)
                        src_i = 2 * dim - 1 - src_i;
                    if (src_j < 0)
                        src_j = -src_j;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    sum_r += src[RIDX(src_i, src_j, dim)].red;
                    sum_g += src[RIDX(src_i, src_j, dim)].green;
                    sum_b += src[RIDX(src_i, src_j, dim)].blue;
                }
            }

            dst[RIDX(i, j, dim)].red = (unsigned short)(sum_r / window_area);
            dst[RIDX(i, j, dim)].green = (unsigned short)(sum_g / window_area);
            dst[RIDX(i, j, dim)].blue = (unsigned short)(sum_b / window_area);
        }
    }
}

/*
 * blur - Your current working version of box blur.
 * IMPORTANT: This is the version you will be graded on
 */
char blur_descr[] = "blur: Current working version";
void blur(int dim, pixel *src, pixel *dst)
{
    naive_blur(dim, src, dst);
}

/*********************************************************************
 * register_blur_functions - Register all of your different versions
 *     of the blur kernel with the driver by calling the
 *     add_blur_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.
 *********************************************************************/

void register_blur_functions()
{
    add_blur_function(&blur, blur_descr);
    add_blur_function(&naive_blur, naive_blur_descr);
    /* ... Register additional test functions here */
}
