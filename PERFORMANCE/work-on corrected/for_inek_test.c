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

    "Deniz Karakoyun",     /* First member full name */
    "e258067@metu.edu.tr", /* First member email address */

    "Simon N. Nganga",     /* Second member full name (leave blank if none) */
    "e278084@metu.edu.tr", /* Second member email addr (leave blank if none) */

    "Can Okulmuş",        /* Third member full name (leave blank if none) */
    "e240026@metu.edu.tr" /* Third member email addr (leave blank if none) */
};

/***************
 * BILATERAL FILTER KERNEL
 ***************/

/******************************************************
 * Your different versions of the bilateral filter kernel go here
 ******************************************************/

/* Global parameters for bilateral filter */
static int radius = 2;           /* Filter radius (window size is 2*radius+1) */
static double sigma_s = 2.0;     /* Spatial parameter */
static double sigma_r = 50000.0; /* Range parameter */




char bilateral_descr0[] =
    "bilateral: more ı want more before if reversing";

void bilateral0(int dim, pixel *src, pixel *dst)
{
    // the constants I will use during computation
    const int r = radius;
    const int BLOCK = 32;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    // all declarations since it is ANSI C
    int i_blocked, j_blocked, i_in_block, j_in_block, i, j, u, v;

    // Precompute spatial weights
    double spatial_weights[2 * r + 1][2 * r + 1];
    for (u = -r; u <= r; u++)
    {
        for (v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (aggressive 4x unroll) =================
    for (i_blocked = r; i_blocked < block_limit; i_blocked += BLOCK)
    {
        int i_end = (i_blocked + BLOCK < block_limit) ? i_blocked + BLOCK : block_limit;
        for (j_blocked = r; j_blocked < block_limit; j_blocked += BLOCK)
        {

            int j_end = (j_blocked + BLOCK < block_limit) ? j_blocked + BLOCK : block_limit;

            for (i_in_block = i_blocked; i_in_block < i_end; i_in_block++)
            {
                for (j_in_block = j_blocked; j_in_block < j_end; j_in_block++)
                {
                    int center_idx = i_in_block * dim + j_in_block;
                    pixel center = src[center_idx];
                    double center_red = center.red, center_green = center.green, center_blue = center.blue;
                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;
                    double inv_wsum;

                    for (u = -r; u <= r; u++)
                    {
                        // taking the begining of the window so that we can use positive indexing
                        pixel *row = src + (i_in_block + u) * dim + (j_in_block - r);
                        double *sw_row = spatial_weights[u + r];

                        // positive indexing since pointer points to the most far away element
                        for (v = 0; v + 3 < window_width; v += 4)
                        {
                            pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                            double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                            double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                            double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                            double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                            double color_sq0 = difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0;
                            double color_sq1 = difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1;
                            double color_sq2 = difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2;
                            double color_sq3 = difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3;

                            double weight0 = sw_row[v] * exp(-color_sq0 * inv_two_sigma_r_sq);
                            double weight1 = sw_row[v + 1] * exp(-color_sq1 * inv_two_sigma_r_sq);
                            double weight2 = sw_row[v + 2] * exp(-color_sq2 * inv_two_sigma_r_sq);
                            double weight3 = sw_row[v + 3] * exp(-color_sq3 * inv_two_sigma_r_sq);

                            weight_sum += weight0 + weight1 + weight2 + weight3;
                            sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                            sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                            sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                        }

                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                            double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                            weight_sum += weight;
                            sum_r += p.red * weight;
                            sum_g += p.green * weight;
                            sum_b += p.blue * weight;
                        }
                    }

                    // multiplication is cheaper than division so first divide after use it as multiplicatiob
                    inv_wsum = 1.0 / weight_sum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                // mostly src_i will be >=0 so dont waste cycle
                if (src_i < 0)
                    src_i = -src_i;
                // int row = src_i * dim;

                for (v = -r; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                int v = 0;
                for (; v + 3 < window_width; v += 4)
                {
                    pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = 2 * dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = 2 * dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = 2 * dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = 2 * dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = 2 * dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = 2 * dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = 2 * dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = 2 * dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = 2 * dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = 2 * dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = 2 * dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = 2 * dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = src[src_i * dim + src_j0];
                    pixel src_i_u_j_v1 = src[src_i * dim + src_j1];
                    pixel src_i_u_j_v2 = src[src_i * dim + src_j2];
                    pixel src_i_u_j_v3 = src[src_i * dim + src_j3];

                    double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                    double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                    double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                    double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                    double weight0 = spatial_weights[u + r][v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = spatial_weights[u + r][v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = spatial_weights[u + r][v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = spatial_weights[u + r][v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                    sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                    sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}

























char bilateral_descr[] =
    "BASE BILATERAL MAKE CHANGES ON THIS ONE THIS IS THE BEST YET ANDD SUBMIT FUNC BLOCK 32";

void bilateral(int dim, pixel *src, pixel *dst)
{
    // the constants I will use during computation
    const int r = radius;
    const int BLOCK = 32;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    // all declarations since it is ANSI C
    int i_blocked, j_blocked, i_in_block, j_in_block, i, j, u, v;

    // Precompute spatial weights
    double spatial_weights[2 * r + 1][2 * r + 1];
    for (u = -r; u <= r; u++)
    {
        for (v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (aggressive 4x unroll) =================
    for (i_blocked = r; i_blocked < block_limit; i_blocked += BLOCK)
    {
        int i_end = (i_blocked + BLOCK < block_limit) ? i_blocked + BLOCK : block_limit;
        for (j_blocked = r; j_blocked < block_limit; j_blocked += BLOCK)
        {

            int j_end = (j_blocked + BLOCK < block_limit) ? j_blocked + BLOCK : block_limit;

            for (i_in_block = i_blocked; i_in_block < i_end; i_in_block++)
            {
                for (j_in_block = j_blocked; j_in_block < j_end; j_in_block++)
                {
                    int center_idx = i_in_block * dim + j_in_block;
                    pixel center = src[center_idx];
                    double center_red = center.red, center_green = center.green, center_blue = center.blue;
                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;
                    double inv_wsum;

                    for (u = -r; u <= r; u++)
                    {
                        // taking the begining of the window so that we can use positive indexing
                        pixel *row = src + (i_in_block + u) * dim + (j_in_block - r);
                        double *sw_row = spatial_weights[u + r];

                        // positive indexing since pointer points to the most far away element
                        for (v = 0; v + 3 < window_width; v += 4)
                        {
                            pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                            double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                            double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                            double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                            double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                            double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                            double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                            double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                            double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                            weight_sum += weight0 + weight1 + weight2 + weight3;
                            sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                            sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                            sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                        }

                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                            double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                            weight_sum += weight;
                            sum_r += p.red * weight;
                            sum_g += p.green * weight;
                            sum_b += p.blue * weight;
                        }
                    }

                    // multiplication is cheaper than division so first divide after use it as multiplicatiob
                    inv_wsum = 1.0 / weight_sum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                // mostly src_i will be >=0 so dont waste cycle
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (v = -r; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                v = 0;
                for (; v + 3 < window_width; v += 4)
                {

                    pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;
                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}











/*
                    BLOCK SIZE 16
                    16 16 16 16 16 16 16 16 16 16 16 16 16 16 16
                    */

char base_descr_16[] =
    "BASE BILATERAL MAKE CHANGES ON THIS ONE BEST FUNCTION VERSION BLOCK 16";

void base_bilateral_16(int dim, pixel *src, pixel *dst)
{
    // the constants I will use during computation
    const int r = radius;
    const int BLOCK = 16;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    // all declarations since it is ANSI C
    int i_blocked, j_blocked, i_in_block, j_in_block, i, j, u, v;

    // Precompute spatial weights
    double spatial_weights[2 * r + 1][2 * r + 1];
    for (u = -r; u <= r; u++)
    {
        for (v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (aggressive 4x unroll) =================
    for (i_blocked = r; i_blocked < block_limit; i_blocked += BLOCK)
    {
        int i_end = (i_blocked + BLOCK < block_limit) ? i_blocked + BLOCK : block_limit;
        for (j_blocked = r; j_blocked < block_limit; j_blocked += BLOCK)
        {

            int j_end = (j_blocked + BLOCK < block_limit) ? j_blocked + BLOCK : block_limit;

            for (i_in_block = i_blocked; i_in_block < i_end; i_in_block++)
            {
                for (j_in_block = j_blocked; j_in_block < j_end; j_in_block++)
                {
                    int center_idx = i_in_block * dim + j_in_block;
                    pixel center = src[center_idx];
                    double center_red = center.red, center_green = center.green, center_blue = center.blue;
                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;
                    double inv_wsum;

                    for (u = -r; u <= r; u++)
                    {
                        // taking the begining of the window so that we can use positive indexing
                        pixel *row = src + (i_in_block + u) * dim + (j_in_block - r);
                        double *sw_row = spatial_weights[u + r];

                        // positive indexing since pointer points to the most far away element
                        for (v = 0; v + 3 < window_width; v += 4)
                        {
                            pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                            double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                            double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                            double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                            double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                            double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                            double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                            double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                            double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                            weight_sum += weight0 + weight1 + weight2 + weight3;
                            sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                            sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                            sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                        }

                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                            double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                            weight_sum += weight;
                            sum_r += p.red * weight;
                            sum_g += p.green * weight;
                            sum_b += p.blue * weight;
                        }
                    }

                    // multiplication is cheaper than division so first divide after use it as multiplicatiob
                    inv_wsum = 1.0 / weight_sum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                // mostly src_i will be >=0 so dont waste cycle
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (v = -r; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                v = 0;
                for (; v + 3 < window_width; v += 4)
                {

                    pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;
                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}



















/*
                    BLOCK SIZE 64 64 64 64 64 64 64 64
                    64 64 64 64 64 64 64 64 64 64
                    */

char base_descr_64[] =
    "BASE BILATERAL MAKE CHANGES ON THIS ONE BEST FUNCTION VERSION BLOCK 64";

void base_bilateral_64(int dim, pixel *src, pixel *dst)
{
    // the constants I will use during computation
    const int r = radius;
    const int BLOCK = 64;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    // all declarations since it is ANSI C
    int i_blocked, j_blocked, i_in_block, j_in_block, i, j, u, v;

    // Precompute spatial weights
    double spatial_weights[2 * r + 1][2 * r + 1];
    for (u = -r; u <= r; u++)
    {
        for (v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (aggressive 4x unroll) =================
    for (i_blocked = r; i_blocked < block_limit; i_blocked += BLOCK)
    {
        int i_end = (i_blocked + BLOCK < block_limit) ? i_blocked + BLOCK : block_limit;
        for (j_blocked = r; j_blocked < block_limit; j_blocked += BLOCK)
        {

            int j_end = (j_blocked + BLOCK < block_limit) ? j_blocked + BLOCK : block_limit;

            for (i_in_block = i_blocked; i_in_block < i_end; i_in_block++)
            {
                for (j_in_block = j_blocked; j_in_block < j_end; j_in_block++)
                {
                    int center_idx = i_in_block * dim + j_in_block;
                    pixel center = src[center_idx];
                    double center_red = center.red, center_green = center.green, center_blue = center.blue;
                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;
                    double inv_wsum;

                    for (u = -r; u <= r; u++)
                    {
                        // taking the begining of the window so that we can use positive indexing
                        pixel *row = src + (i_in_block + u) * dim + (j_in_block - r);
                        double *sw_row = spatial_weights[u + r];

                        // positive indexing since pointer points to the most far away element
                        for (v = 0; v + 3 < window_width; v += 4)
                        {
                            pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                            double difference_red0 = center_red - src_i_u_j_v0.red, difference_green0 = center_green - src_i_u_j_v0.green, difference_blue0 = center_blue - src_i_u_j_v0.blue;
                            double difference_red1 = center_red - src_i_u_j_v1.red, difference_green1 = center_green - src_i_u_j_v1.green, difference_blue1 = center_blue - src_i_u_j_v1.blue;
                            double difference_red2 = center_red - src_i_u_j_v2.red, difference_green2 = center_green - src_i_u_j_v2.green, difference_blue2 = center_blue - src_i_u_j_v2.blue;
                            double difference_red3 = center_red - src_i_u_j_v3.red, difference_green3 = center_green - src_i_u_j_v3.green, difference_blue3 = center_blue - src_i_u_j_v3.blue;

                            double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                            double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                            double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                            double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                            weight_sum += weight0 + weight1 + weight2 + weight3;
                            sum_r += src_i_u_j_v0.red * weight0 + src_i_u_j_v1.red * weight1 + src_i_u_j_v2.red * weight2 + src_i_u_j_v3.red * weight3;
                            sum_g += src_i_u_j_v0.green * weight0 + src_i_u_j_v1.green * weight1 + src_i_u_j_v2.green * weight2 + src_i_u_j_v3.green * weight3;
                            sum_b += src_i_u_j_v0.blue * weight0 + src_i_u_j_v1.blue * weight1 + src_i_u_j_v2.blue * weight2 + src_i_u_j_v3.blue * weight3;
                        }

                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                            double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                            weight_sum += weight;
                            sum_r += p.red * weight;
                            sum_g += p.green * weight;
                            sum_b += p.blue * weight;
                        }
                    }

                    // multiplication is cheaper than division so first divide after use it as multiplicatiob
                    inv_wsum = 1.0 / weight_sum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                // mostly src_i will be >=0 so dont waste cycle
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (v = -r; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                v = 0;
                for (; v + 3 < window_width; v += 4)
                {

                    pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;
                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}

















char bilateral_optimized_descr[] =
    "bilateral: optimized with precomputation and better memory access";

void bilateral_optimized(int dim, pixel *src, pixel *dst)
{
    const int r = radius;
    const int BLOCK = 32;
    const double two_sigma_s_sq = 2.0 * sigma_s * sigma_s;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    // all declarations since it is ANSI C
    int i_blocked, j_blocked, i_in_block, j_in_block, i, j, u, v;
    const int block_limit = dim - r;
    const int window_width = r + r + 1;
    // Precompute spatial weights
    double spatial_weights[2 * r + 1][2 * r + 1];
    for (int u = -r; u <= r; u++)
    {
        for (int v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / two_sigma_s_sq);
        }
    }

    // ================= INTERIOR (optimized) =================
    for (int ii = r; ii < dim - r; ii += BLOCK)
    {
        int i_end = (ii + BLOCK < dim - r) ? ii + BLOCK : dim - r;
        for (int jj = r; jj < dim - r; jj += BLOCK)
        {
            int j_end = (jj + BLOCK < dim - r) ? jj + BLOCK : dim - r;

            for (int i = ii; i < i_end; i++)
            {
                for (int j = jj; j < j_end; j++)
                {
                    int center_idx = i * dim + j;
                    pixel center = src[center_idx];
                    double cr = center.red;
                    double cg = center.green;
                    double cb = center.blue;

                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
                    double wsum = 0.0;

                    // Process neighborhood with better memory access pattern
                    for (int u = -r; u <= r; u++)
                    {
                        int row_idx = (i + u) * dim + j - r;
                        pixel *row = src + row_idx;
                        double *sw_row = spatial_weights[u + r];

                        // Manual unroll by 4 for better ILP
                        int v = 0;
                        int window_width = 2 * r + 1;

                        for (; v + 3 < window_width; v += 4)
                        {
                            // Process 4 pixels at once
                            pixel p0 = row[v];
                            pixel p1 = row[v + 1];
                            pixel p2 = row[v + 2];
                            pixel p3 = row[v + 3];

                            double dr0 = cr - p0.red;
                            double dg0 = cg - p0.green;
                            double db0 = cb - p0.blue;
                            double color_diff_sq0 = dr0 * dr0 + dg0 * dg0 + db0 * db0;
                            double w0 = sw_row[v] * exp(-color_diff_sq0 * inv_two_sigma_r_sq);

                            double dr1 = cr - p1.red;
                            double dg1 = cg - p1.green;
                            double db1 = cb - p1.blue;
                            double color_diff_sq1 = dr1 * dr1 + dg1 * dg1 + db1 * db1;
                            double w1 = sw_row[v + 1] * exp(-color_diff_sq1 * inv_two_sigma_r_sq);

                            double dr2 = cr - p2.red;
                            double dg2 = cg - p2.green;
                            double db2 = cb - p2.blue;
                            double color_diff_sq2 = dr2 * dr2 + dg2 * dg2 + db2 * db2;
                            double w2 = sw_row[v + 2] * exp(-color_diff_sq2 * inv_two_sigma_r_sq);

                            double dr3 = cr - p3.red;
                            double dg3 = cg - p3.green;
                            double db3 = cb - p3.blue;
                            double color_diff_sq3 = dr3 * dr3 + dg3 * dg3 + db3 * db3;
                            double w3 = sw_row[v + 3] * exp(-color_diff_sq3 * inv_two_sigma_r_sq);

                            wsum += w0 + w1 + w2 + w3;
                            sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                            sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                            sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                        }

                        // Handle remaining pixels
                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double dr = cr - p.red;
                            double dg = cg - p.green;
                            double db = cb - p.blue;
                            double color_diff_sq = dr * dr + dg * dg + db * db;
                            double w = sw_row[v] * exp(-color_diff_sq * inv_two_sigma_r_sq);

                            wsum += w;
                            sum_r += p.red * w;
                            sum_g += p.green * w;
                            sum_b += p.blue * w;
                        }
                    }

                    double inv_wsum = 1.0 / wsum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

// ================= TOP-LEFT CORNER - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                // mostly src_i will be >=0 so dont waste cycle
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (v = -r; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (i = 0; i < r; i++)
    {
        for (j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                v = 0;
                for (; v + 3 < window_width; v += 4)
                {

                    pixel src_i_u_j_v0 = row[v], src_i_u_j_v1 = row[v + 1], src_i_u_j_v2 = row[v + 2], src_i_u_j_v3 = row[v + 3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = sw_row[v] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;
                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 < 0)
                        src_j0 = -src_j0;
                    if (src_j1 < 0)
                        src_j1 = -src_j1;
                    if (src_j2 < 0)
                        src_j2 = -src_j2;
                    if (src_j3 < 0)
                        src_j3 = -src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================

    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            double center_red = center.red, center_green = center.green, center_blue = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = dim + dim - 1 - src_i;

                int v = -r;
                pixel *row = src + src_i * dim;
                double *sw_row = spatial_weights[u + r];
                for (; v + 3 <= r; v += 4)
                {
                    int src_j0 = j + v, src_j1 = j + v + 1, src_j2 = j + v + 2, src_j3 = j + v + 3;
                    if (src_j0 >= dim)
                        src_j0 = dim + dim - 1 - src_j0;
                    if (src_j1 >= dim)
                        src_j1 = dim + dim - 1 - src_j1;
                    if (src_j2 >= dim)
                        src_j2 = dim + dim - 1 - src_j2;
                    if (src_j3 >= dim)
                        src_j3 = dim + dim - 1 - src_j3;

                    pixel src_i_u_j_v0 = row[src_j0], src_i_u_j_v1 = row[src_j1], src_i_u_j_v2 = row[src_j2], src_i_u_j_v3 = row[src_j3];

                    unsigned short src_i_u_j_v0_red = src_i_u_j_v0.red, src_i_u_j_v0_green = src_i_u_j_v0.green, src_i_u_j_v0_blue = src_i_u_j_v0.blue;
                    unsigned short src_i_u_j_v1_red = src_i_u_j_v1.red, src_i_u_j_v1_green = src_i_u_j_v1.green, src_i_u_j_v1_blue = src_i_u_j_v1.blue;
                    unsigned short src_i_u_j_v2_red = src_i_u_j_v2.red, src_i_u_j_v2_green = src_i_u_j_v2.green, src_i_u_j_v2_blue = src_i_u_j_v2.blue;
                    unsigned short src_i_u_j_v3_red = src_i_u_j_v3.red, src_i_u_j_v3_green = src_i_u_j_v3.green, src_i_u_j_v3_blue = src_i_u_j_v3.blue;

                    double difference_red0 = center_red - src_i_u_j_v0_red, difference_green0 = center_green - src_i_u_j_v0_green, difference_blue0 = center_blue - src_i_u_j_v0_blue;
                    double difference_red1 = center_red - src_i_u_j_v1_red, difference_green1 = center_green - src_i_u_j_v1_green, difference_blue1 = center_blue - src_i_u_j_v1_blue;
                    double difference_red2 = center_red - src_i_u_j_v2_red, difference_green2 = center_green - src_i_u_j_v2_green, difference_blue2 = center_blue - src_i_u_j_v2_blue;
                    double difference_red3 = center_red - src_i_u_j_v3_red, difference_green3 = center_green - src_i_u_j_v3_green, difference_blue3 = center_blue - src_i_u_j_v3_blue;

                    double weight0 = sw_row[v + r] * exp(-(difference_red0 * difference_red0 + difference_green0 * difference_green0 + difference_blue0 * difference_blue0) * inv_two_sigma_r_sq);
                    double weight1 = sw_row[v + r + 1] * exp(-(difference_red1 * difference_red1 + difference_green1 * difference_green1 + difference_blue1 * difference_blue1) * inv_two_sigma_r_sq);
                    double weight2 = sw_row[v + r + 2] * exp(-(difference_red2 * difference_red2 + difference_green2 * difference_green2 + difference_blue2 * difference_blue2) * inv_two_sigma_r_sq);
                    double weight3 = sw_row[v + r + 3] * exp(-(difference_red3 * difference_red3 + difference_green3 * difference_green3 + difference_blue3 * difference_blue3) * inv_two_sigma_r_sq);

                    weight_sum += weight0 + weight1 + weight2 + weight3;
                    sum_r += src_i_u_j_v0_red * weight0 + src_i_u_j_v1_red * weight1 + src_i_u_j_v2_red * weight2 + src_i_u_j_v3_red * weight3;
                    sum_g += src_i_u_j_v0_green * weight0 + src_i_u_j_v1_green * weight1 + src_i_u_j_v2_green * weight2 + src_i_u_j_v3_green * weight3;
                    sum_b += src_i_u_j_v0_blue * weight0 + src_i_u_j_v1_blue * weight1 + src_i_u_j_v2_blue * weight2 + src_i_u_j_v3_blue * weight3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = dim + dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double difference_red = center_red - p.red, difference_green = center_green - p.green, difference_blue = center_blue - p.blue;
                    double weight = spatial_weights[u + r][v + r] * exp(-(difference_red * difference_red + difference_green * difference_green + difference_blue * difference_blue) * inv_two_sigma_r_sq);
                    weight_sum += weight;
                    sum_r += p.red * weight;
                    sum_g += p.green * weight;
                    sum_b += p.blue * weight;
                }
            }

            double inv_wsum = 1.0 / weight_sum;
            dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
            dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
            dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}

char bilateral_split_borders_descr[] =
    "bilateral: split borders, no if checks, aggressive optimization";

void bilateral_split_borders(int dim, pixel *src, pixel *dst)
{
    const int r = radius;
    const int BLOCK = 32;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;

    // Precompute spatial weights
    double spatial_weights[21][21];
    for (int u = -r; u <= r; u++)
    {
        for (int v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (fully optimized, no boundary checks) =================
    for (int ii = r; ii < block_limit; ii += BLOCK)
    {
        for (int jj = r; jj < block_limit; jj += BLOCK)
        {
            int i_end = (ii + BLOCK < block_limit) ? ii + BLOCK : block_limit;
            int j_end = (jj + BLOCK < block_limit) ? jj + BLOCK : block_limit;

            for (int i = ii; i < i_end; i++)
            {
                for (int j = jj; j < j_end; j++)
                {
                    int center_idx = i * dim + j;
                    pixel center = src[center_idx];
                    double cr = center.red, cg = center.green, cb = center.blue;

                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

                    for (int u = -r; u <= r; u++)
                    {
                        pixel *row = src + (i + u) * dim + (j - r);
                        double *sw_row = spatial_weights[u + r];

                        int v = 0;
                        for (; v + 3 < 2 * r + 1; v += 4)
                        {
                            pixel p0 = row[v], p1 = row[v + 1], p2 = row[v + 2], p3 = row[v + 3];

                            double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                            double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                            double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                            double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                            double w0 = sw_row[v] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                            double w1 = sw_row[v + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                            double w2 = sw_row[v + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                            double w3 = sw_row[v + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                            wsum += w0 + w1 + w2 + w3;
                            sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                            sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                            sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                        }

                        for (; v < 2 * r + 1; v++)
                        {
                            pixel p = row[v];
                            double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                            double w = sw_row[v] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                            wsum += w;
                            sum_r += p.red * w;
                            sum_g += p.green * w;
                            sum_b += p.blue * w;
                        }
                    }

                    double inv_wsum = 1.0 / wsum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER (i < r, j < r) =================
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE (i < r, r <= j < block_limit) =================
    for (int i = 0; i < r; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v; // No boundary check needed

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER (i < r, j >= block_limit) =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE (r <= i < block_limit, j < r) =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u; // No boundary check needed

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE (r <= i < block_limit, j >= block_limit) =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u; // No boundary check needed

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v; // No boundary check needed

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}

char bilateral_max_performance_descr0[] =
    "bilateral: maximum performance - 4x unroll everywhere, optimized for Xeon E-2224";

void bilateral_max_performance(int dim, pixel *src, pixel *dst)
{
    const int r = radius;
    const int BLOCK = 32;
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    // Precompute spatial weights
    double spatial_weights[21][21];
    for (int u = -r; u <= r; u++)
    {
        for (int v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (aggressive 4x unroll) =================
    for (int ii = r; ii < block_limit; ii += BLOCK)
    {
        for (int jj = r; jj < block_limit; jj += BLOCK)
        {
            int i_end = (ii + BLOCK < block_limit) ? ii + BLOCK : block_limit;
            int j_end = (jj + BLOCK < block_limit) ? jj + BLOCK : block_limit;

            for (int i = ii; i < i_end; i++)
            {
                for (int j = jj; j < j_end; j++)
                {
                    int center_idx = i * dim + j;
                    pixel center = src[center_idx];
                    double cr = center.red, cg = center.green, cb = center.blue;

                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

                    for (int u = -r; u <= r; u++)
                    {
                        pixel *row = src + (i + u) * dim + (j - r);
                        double *sw_row = spatial_weights[u + r];

                        int v = 0;
                        for (; v + 3 < window_width; v += 4)
                        {
                            pixel p0 = row[v], p1 = row[v + 1], p2 = row[v + 2], p3 = row[v + 3];

                            double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                            double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                            double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                            double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                            double color_sq0 = dr0 * dr0 + dg0 * dg0 + db0 * db0;
                            double color_sq1 = dr1 * dr1 + dg1 * dg1 + db1 * db1;
                            double color_sq2 = dr2 * dr2 + dg2 * dg2 + db2 * db2;
                            double color_sq3 = dr3 * dr3 + dg3 * dg3 + db3 * db3;

                            double w0 = sw_row[v] * exp(-color_sq0 * inv_two_sigma_r_sq);
                            double w1 = sw_row[v + 1] * exp(-color_sq1 * inv_two_sigma_r_sq);
                            double w2 = sw_row[v + 2] * exp(-color_sq2 * inv_two_sigma_r_sq);
                            double w3 = sw_row[v + 3] * exp(-color_sq3 * inv_two_sigma_r_sq);

                            wsum += w0 + w1 + w2 + w3;
                            sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                            sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                            sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                        }

                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                            double w = sw_row[v] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                            wsum += w;
                            sum_r += p.red * w;
                            sum_g += p.green * w;
                            sum_b += p.blue * w;
                        }
                    }

                    double inv_wsum = 1.0 / wsum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 < 0)
                        sj0 = -sj0;
                    if (sj1 < 0)
                        sj1 = -sj1;
                    if (sj2 < 0)
                        sj2 = -sj2;
                    if (sj3 < 0)
                        sj3 = -sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                int v = 0;
                for (; v + 3 < window_width; v += 4)
                {
                    pixel p0 = row[v], p1 = row[v + 1], p2 = row[v + 2], p3 = row[v + 3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = sw_row[v] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = sw_row[v + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = sw_row[v + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = sw_row[v + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = sw_row[v] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 >= dim)
                        sj0 = 2 * dim - 1 - sj0;
                    if (sj1 >= dim)
                        sj1 = 2 * dim - 1 - sj1;
                    if (sj2 >= dim)
                        sj2 = 2 * dim - 1 - sj2;
                    if (sj3 >= dim)
                        sj3 = 2 * dim - 1 - sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 < 0)
                        sj0 = -sj0;
                    if (sj1 < 0)
                        sj1 = -sj1;
                    if (sj2 < 0)
                        sj2 = -sj2;
                    if (sj3 < 0)
                        sj3 = -sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 >= dim)
                        sj0 = 2 * dim - 1 - sj0;
                    if (sj1 >= dim)
                        sj1 = 2 * dim - 1 - sj1;
                    if (sj2 >= dim)
                        sj2 = 2 * dim - 1 - sj2;
                    if (sj3 >= dim)
                        sj3 = 2 * dim - 1 - sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v; // No boundary check needed

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}

char bilateral_xeon_optimized_descr[] =
    "bilateral: highly optimized for Xeon E-2224 - 8x unroll + prefetch + cache opt";

void bilateral_xeon_optimized(int dim, pixel *src, pixel *dst)
{
    const int r = radius;
    const int BLOCK = 32; // Larger blocks for better cache utilization
    const double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;
    const double inv_two_sigma_r_sq = 1.0 / two_sigma_r_sq;
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    // Precompute spatial weights - align for better cache performance
    double spatial_weights[21][21] __attribute__((aligned(64)));
    for (int u = -r; u <= r; u++)
    {
        for (int v = -r; v <= r; v++)
        {
            spatial_weights[u + r][v + r] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    // ================= INTERIOR (aggressive 8x unroll with prefetch) =================
    for (int ii = r; ii < block_limit; ii += BLOCK)
    {
        for (int jj = r; jj < block_limit; jj += BLOCK)
        {
            int i_end = (ii + BLOCK < block_limit) ? ii + BLOCK : block_limit;
            int j_end = (jj + BLOCK < block_limit) ? jj + BLOCK : block_limit;

            for (int i = ii; i < i_end; i++)
            {
                for (int j = jj; j < j_end; j++)
                {
                    int center_idx = i * dim + j;
                    pixel center = src[center_idx];
                    double cr = center.red, cg = center.green, cb = center.blue;

                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

                    for (int u = -r; u <= r; u++)
                    {
                        pixel *row = src + (i + u) * dim + (j - r);
                        double *sw_row = spatial_weights[u + r];

                        // Prefetch next row for better cache utilization
                        if (u < r)
                        {
                            __builtin_prefetch(src + (i + u + 1) * dim + (j - r), 0, 3);
                        }

                        int v = 0;
                        // 8x unroll for maximum throughput
                        for (; v + 7 < window_width; v += 8)
                        {
                            pixel p0 = row[v], p1 = row[v + 1], p2 = row[v + 2], p3 = row[v + 3];
                            pixel p4 = row[v + 4], p5 = row[v + 5], p6 = row[v + 6], p7 = row[v + 7];

                            // Compute color differences
                            double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                            double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                            double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                            double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;
                            double dr4 = cr - p4.red, dg4 = cg - p4.green, db4 = cb - p4.blue;
                            double dr5 = cr - p5.red, dg5 = cg - p5.green, db5 = cb - p5.blue;
                            double dr6 = cr - p6.red, dg6 = cg - p6.green, db6 = cb - p6.blue;
                            double dr7 = cr - p7.red, dg7 = cg - p7.green, db7 = cb - p7.blue;

                            // Compute color distance squared
                            double color_sq0 = dr0 * dr0 + dg0 * dg0 + db0 * db0;
                            double color_sq1 = dr1 * dr1 + dg1 * dg1 + db1 * db1;
                            double color_sq2 = dr2 * dr2 + dg2 * dg2 + db2 * db2;
                            double color_sq3 = dr3 * dr3 + dg3 * dg3 + db3 * db3;
                            double color_sq4 = dr4 * dr4 + dg4 * dg4 + db4 * db4;
                            double color_sq5 = dr5 * dr5 + dg5 * dg5 + db5 * db5;
                            double color_sq6 = dr6 * dr6 + dg6 * dg6 + db6 * db6;
                            double color_sq7 = dr7 * dr7 + dg7 * dg7 + db7 * db7;

                            // Compute weights
                            double w0 = sw_row[v] * exp(-color_sq0 * inv_two_sigma_r_sq);
                            double w1 = sw_row[v + 1] * exp(-color_sq1 * inv_two_sigma_r_sq);
                            double w2 = sw_row[v + 2] * exp(-color_sq2 * inv_two_sigma_r_sq);
                            double w3 = sw_row[v + 3] * exp(-color_sq3 * inv_two_sigma_r_sq);
                            double w4 = sw_row[v + 4] * exp(-color_sq4 * inv_two_sigma_r_sq);
                            double w5 = sw_row[v + 5] * exp(-color_sq5 * inv_two_sigma_r_sq);
                            double w6 = sw_row[v + 6] * exp(-color_sq6 * inv_two_sigma_r_sq);
                            double w7 = sw_row[v + 7] * exp(-color_sq7 * inv_two_sigma_r_sq);

                            // Accumulate (use FMA-friendly pattern)
                            wsum += (w0 + w1) + (w2 + w3) + (w4 + w5) + (w6 + w7);
                            sum_r += (p0.red * w0 + p1.red * w1) + (p2.red * w2 + p3.red * w3) +
                                     (p4.red * w4 + p5.red * w5) + (p6.red * w6 + p7.red * w7);
                            sum_g += (p0.green * w0 + p1.green * w1) + (p2.green * w2 + p3.green * w3) +
                                     (p4.green * w4 + p5.green * w5) + (p6.green * w6 + p7.green * w7);
                            sum_b += (p0.blue * w0 + p1.blue * w1) + (p2.blue * w2 + p3.blue * w3) +
                                     (p4.blue * w4 + p5.blue * w5) + (p6.blue * w6 + p7.blue * w7);
                        }

                        // Handle remaining elements
                        for (; v < window_width; v++)
                        {
                            pixel p = row[v];
                            double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                            double w = sw_row[v] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                            wsum += w;
                            sum_r += p.red * w;
                            sum_g += p.green * w;
                            sum_b += p.blue * w;
                        }
                    }

                    double inv_wsum = 1.0 / wsum;
                    dst[center_idx].red = (unsigned short)(sum_r * inv_wsum);
                    dst[center_idx].green = (unsigned short)(sum_g * inv_wsum);
                    dst[center_idx].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    // ================= TOP-LEFT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 < 0)
                        sj0 = -sj0;
                    if (sj1 < 0)
                        sj1 = -sj1;
                    if (sj2 < 0)
                        sj2 = -sj2;
                    if (sj3 < 0)
                        sj3 = -sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP EDGE - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                pixel *row = src + src_i * dim + j - r;
                double *sw_row = spatial_weights[u + r];

                int v = 0;
                for (; v + 3 < window_width; v += 4)
                {
                    pixel p0 = row[v], p1 = row[v + 1], p2 = row[v + 2], p3 = row[v + 3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = sw_row[v] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = sw_row[v + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = sw_row[v + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = sw_row[v + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v < window_width; v++)
                {
                    pixel p = row[v];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = sw_row[v] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= TOP-RIGHT CORNER - 4x unrolled =================
    for (int i = 0; i < r; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 >= dim)
                        sj0 = 2 * dim - 1 - sj0;
                    if (sj1 >= dim)
                        sj1 = 2 * dim - 1 - sj1;
                    if (sj2 >= dim)
                        sj2 = 2 * dim - 1 - sj2;
                    if (sj3 >= dim)
                        sj3 = 2 * dim - 1 - sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= LEFT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 < 0)
                        sj0 = -sj0;
                    if (sj1 < 0)
                        sj1 = -sj1;
                    if (sj2 < 0)
                        sj2 = -sj2;
                    if (sj3 < 0)
                        sj3 = -sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= RIGHT EDGE - 4x unrolled =================
    for (int i = r; i < block_limit; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;

                int v = -r;
                for (; v + 3 <= r; v += 4)
                {
                    int sj0 = j + v, sj1 = j + v + 1, sj2 = j + v + 2, sj3 = j + v + 3;
                    if (sj0 >= dim)
                        sj0 = 2 * dim - 1 - sj0;
                    if (sj1 >= dim)
                        sj1 = 2 * dim - 1 - sj1;
                    if (sj2 >= dim)
                        sj2 = 2 * dim - 1 - sj2;
                    if (sj3 >= dim)
                        sj3 = 2 * dim - 1 - sj3;

                    pixel p0 = src[src_i * dim + sj0];
                    pixel p1 = src[src_i * dim + sj1];
                    pixel p2 = src[src_i * dim + sj2];
                    pixel p3 = src[src_i * dim + sj3];

                    double dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                    double dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                    double dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                    double dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                    double w0 = spatial_weights[u + r][v + r] * exp(-(dr0 * dr0 + dg0 * dg0 + db0 * db0) * inv_two_sigma_r_sq);
                    double w1 = spatial_weights[u + r][v + r + 1] * exp(-(dr1 * dr1 + dg1 * dg1 + db1 * db1) * inv_two_sigma_r_sq);
                    double w2 = spatial_weights[u + r][v + r + 2] * exp(-(dr2 * dr2 + dg2 * dg2 + db2 * db2) * inv_two_sigma_r_sq);
                    double w3 = spatial_weights[u + r][v + r + 3] * exp(-(dr3 * dr3 + dg3 * dg3 + db3 * db3) * inv_two_sigma_r_sq);

                    wsum += w0 + w1 + w2 + w3;
                    sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                    sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                    sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                }

                for (; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-LEFT CORNER (i >= block_limit, j < r) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = 0; j < r; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM EDGE (i >= block_limit, r <= j < block_limit) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = r; j < block_limit; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v; // No boundary check needed

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }

    // ================= BOTTOM-RIGHT CORNER (i >= block_limit, j >= block_limit) =================
    for (int i = block_limit; i < dim; i++)
    {
        for (int j = block_limit; j < dim; j++)
        {
            int index = i * dim + j;
            pixel center = src[index];
            double cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, wsum = 0.0;

            for (int u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (int v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    double dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    double w = spatial_weights[u + r][v + r] * exp(-(dr * dr + dg * dg + db * db) * inv_two_sigma_r_sq);
                    wsum += w;
                    sum_r += p.red * w;
                    sum_g += p.green * w;
                    sum_b += p.blue * w;
                }
            }

            double inv_wsum = 1.0 / wsum;
            dst[index].red = (unsigned short)(sum_r * inv_wsum);
            dst[index].green = (unsigned short)(sum_g * inv_wsum);
            dst[index].blue = (unsigned short)(sum_b * inv_wsum);
        }
    }
}










void register_bilateral_functions()
{
    add_bilateral_function(&bilateral0, bilateral_descr0); // second best till now 
    add_bilateral_function(&base_bilateral_16, base_descr_16);
    add_bilateral_function(&bilateral, bilateral_descr);// best till now 
    add_bilateral_function(&base_bilateral_64, base_descr_64);
    add_bilateral_function(&bilateral_optimized, bilateral_optimized_descr);
    add_bilateral_function(&bilateral_split_borders, bilateral_split_borders_descr);
    add_bilateral_function(&bilateral_max_performance, bilateral_max_performance_descr0);
    add_bilateral_function(&bilateral_xeon_optimized, bilateral_xeon_optimized_descr);
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
char blur_descr[] = "blur: Optmized function which works with 32*32 Cache Blocking";

void blur(int dim, pixel *src, pixel *dst)
{
    int i, j, ii, jj;
    int BLOCK = 32;

    for (ii = 0; ii < dim; ii += BLOCK)
    {
        for (jj = 0; jj < dim; jj += BLOCK)
        {

            int i_limit = (ii + BLOCK < dim) ? (ii + BLOCK) : dim;
            int j_limit = (jj + BLOCK < dim) ? (jj + BLOCK) : dim;

            for (i = ii; i < i_limit; i++)
            {

                // We precompute row indices
                int r0 = (i - 2) * dim;
                int r1 = (i - 1) * dim;
                int r2 = i * dim;
                int r3 = (i + 1) * dim;
                int r4 = (i + 2) * dim;
                // We put the row starts in an array for the loop to easily access them
                int rows[] = {r0, r1, r2, r3, r4};

                for (j = jj; j < j_limit; j++)
                {

                    int sum_r = 0, sum_g = 0, sum_b = 0;

                    // FAST PATH: Check if we are safe from edges
                    if (i >= 2 && i < dim - 2 && j >= 2 && j < dim - 2)
                    {

                        // Loop over the 5 rows
                        for (int k = 0; k < 5; k++)
                        {
                            // Point to start of row + current column - 2
                            pixel *p = src + rows[k] + (j - 2);

                            // Loop over the 5 columns
                            for (int col = 0; col < 5; col++)
                            {
                                sum_r += p[col].red;
                                sum_g += p[col].green;
                                sum_b += p[col].blue;
                            }
                        }

                        int dst_idx = r2 + j;
                        dst[dst_idx].red = (unsigned short)(sum_r / 25);
                        dst[dst_idx].green = (unsigned short)(sum_g / 25);
                        dst[dst_idx].blue = (unsigned short)(sum_b / 25);
                    }
                    else
                    {
                        // Same approach like our naive function for the Edges
                        for (int u = -2; u <= 2; u++)
                        {
                            for (int v = -2; v <= 2; v++)
                            {
                                int src_i = i + u;
                                int src_j = j + v;
                                if (src_i < 0)
                                    src_i = -src_i;
                                else if (src_i >= dim)
                                    src_i = 2 * dim - 1 - src_i;
                                if (src_j < 0)
                                    src_j = -src_j;
                                else if (src_j >= dim)
                                    src_j = 2 * dim - 1 - src_j;

                                int idx = src_i * dim + src_j;
                                sum_r += src[idx].red;
                                sum_g += src[idx].green;
                                sum_b += src[idx].blue;
                            }
                        }
                        int dst_idx = i * dim + j;
                        dst[dst_idx].red = (unsigned short)(sum_r / 25);
                        dst[dst_idx].green = (unsigned short)(sum_g / 25);
                        dst[dst_idx].blue = (unsigned short)(sum_b / 25);
                    }
                }
            }
        }
    }
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
