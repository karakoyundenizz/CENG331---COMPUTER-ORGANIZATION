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





char base_descr[] =
    "BASE BILATERAL MAKE CHANGES ON THIS ONE";

void base_bilateral(int dim, pixel *src, pixel *dst)
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