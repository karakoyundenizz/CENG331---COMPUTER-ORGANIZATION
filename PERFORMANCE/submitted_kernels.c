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
                    color_diff_sq = ((double)src[RIDX(i, j, dim)].red -
                                     (double)src[RIDX(src_i, src_j, dim)].red) *
                                        ((double)src[RIDX(i, j, dim)].red -
                                         (double)src[RIDX(src_i, src_j, dim)].red) +
                                    ((double)src[RIDX(i, j, dim)].green -
                                     (double)src[RIDX(src_i, src_j, dim)].green) *
                                        ((double)src[RIDX(i, j, dim)].green -
                                         (double)src[RIDX(src_i, src_j, dim)].green) +
                                    ((double)src[RIDX(i, j, dim)].blue -
                                     (double)src[RIDX(src_i, src_j, dim)].blue) *
                                        ((double)src[RIDX(i, j, dim)].blue -
                                         (double)src[RIDX(src_i, src_j, dim)].blue);

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

char bilateral_descr[] = "bilateral: Optimized Blocking (64x64) + Full 5x Unrolling";

void bilateral(int dim, pixel *src, pixel *dst)
{
    const int r = radius;
    const int BLOCK = 64; /* Try larger blocks */
    const double inv_two_sigma_r_sq = 1.0 / (2.0 * sigma_r * sigma_r);
    const int block_limit = dim - r;
    const int window_width = 2 * r + 1;

    int i_blocked, j_blocked, i_in_block, j_in_block, i, j, u, v;

    /* Precompute spatial weights - flatten for better cache access */
    double spatial_weights_flat[25]; /* Max 5x5 window */
    int sw_idx = 0;
    for (u = -r; u <= r; u++)
    {
        for (v = -r; v <= r; v++)
        {
            spatial_weights_flat[sw_idx++] =
                exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));
        }
    }

    /* ================= INTERIOR - Optimized with minimal branches ================= */
    for (i_blocked = r; i_blocked < block_limit; i_blocked += BLOCK)
    {
        int i_end = (i_blocked + BLOCK < block_limit) ? i_blocked + BLOCK : block_limit;
        for (j_blocked = r; j_blocked < block_limit; j_blocked += BLOCK)
        {
            int j_end = (j_blocked + BLOCK < block_limit) ? j_blocked + BLOCK : block_limit;

            for (i_in_block = i_blocked; i_in_block < i_end; i_in_block++)
            {
                pixel *dst_row = dst + i_in_block * dim;
                pixel *center_row = src + i_in_block * dim;

                for (j_in_block = j_blocked; j_in_block < j_end; j_in_block++)
                {
                    pixel center = center_row[j_in_block];
                    int cr = center.red, cg = center.green, cb = center.blue;
                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

                    int sw_base = 0;
                    for (u = -r; u <= r; u++)
                    {
                        pixel *row = src + (i_in_block + u) * dim + (j_in_block - r);

                        /* Process entire row without inner loop variable */
                        v = 0;

                        /* Unroll by 5 (window is 5 wide) */
                        if (window_width == 5)
                        {
                            pixel p0 = row[0], p1 = row[1], p2 = row[2], p3 = row[3], p4 = row[4];

                            int dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                            int dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                            int dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                            int dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;
                            int dr4 = cr - p4.red, dg4 = cg - p4.green, db4 = cb - p4.blue;

                            long long cd0 = (long long)dr0 * dr0 + (long long)dg0 * dg0 + (long long)db0 * db0;
                            long long cd1 = (long long)dr1 * dr1 + (long long)dg1 * dg1 + (long long)db1 * db1;
                            long long cd2 = (long long)dr2 * dr2 + (long long)dg2 * dg2 + (long long)db2 * db2;
                            long long cd3 = (long long)dr3 * dr3 + (long long)dg3 * dg3 + (long long)db3 * db3;
                            long long cd4 = (long long)dr4 * dr4 + (long long)dg4 * dg4 + (long long)db4 * db4;

                            double w0 = spatial_weights_flat[sw_base] * exp(-(double)cd0 * inv_two_sigma_r_sq);
                            double w1 = spatial_weights_flat[sw_base + 1] * exp(-(double)cd1 * inv_two_sigma_r_sq);
                            double w2 = spatial_weights_flat[sw_base + 2] * exp(-(double)cd2 * inv_two_sigma_r_sq);
                            double w3 = spatial_weights_flat[sw_base + 3] * exp(-(double)cd3 * inv_two_sigma_r_sq);
                            double w4 = spatial_weights_flat[sw_base + 4] * exp(-(double)cd4 * inv_two_sigma_r_sq);

                            weight_sum += w0 + w1 + w2 + w3 + w4;
                            sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3 + p4.red * w4;
                            sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3 + p4.green * w4;
                            sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3 + p4.blue * w4;
                        }
                        else
                        {
                            /* General case - 4-way unroll */
                            for (v = 0; v + 3 < window_width; v += 4)
                            {
                                pixel p0 = row[v], p1 = row[v + 1], p2 = row[v + 2], p3 = row[v + 3];

                                int dr0 = cr - p0.red, dg0 = cg - p0.green, db0 = cb - p0.blue;
                                int dr1 = cr - p1.red, dg1 = cg - p1.green, db1 = cb - p1.blue;
                                int dr2 = cr - p2.red, dg2 = cg - p2.green, db2 = cb - p2.blue;
                                int dr3 = cr - p3.red, dg3 = cg - p3.green, db3 = cb - p3.blue;

                                long long cd0 = (long long)dr0 * dr0 + (long long)dg0 * dg0 + (long long)db0 * db0;
                                long long cd1 = (long long)dr1 * dr1 + (long long)dg1 * dg1 + (long long)db1 * db1;
                                long long cd2 = (long long)dr2 * dr2 + (long long)dg2 * dg2 + (long long)db2 * db2;
                                long long cd3 = (long long)dr3 * dr3 + (long long)dg3 * dg3 + (long long)db3 * db3;

                                double w0 = spatial_weights_flat[sw_base + v] * exp(-(double)cd0 * inv_two_sigma_r_sq);
                                double w1 = spatial_weights_flat[sw_base + v + 1] * exp(-(double)cd1 * inv_two_sigma_r_sq);
                                double w2 = spatial_weights_flat[sw_base + v + 2] * exp(-(double)cd2 * inv_two_sigma_r_sq);
                                double w3 = spatial_weights_flat[sw_base + v + 3] * exp(-(double)cd3 * inv_two_sigma_r_sq);

                                weight_sum += w0 + w1 + w2 + w3;
                                sum_r += p0.red * w0 + p1.red * w1 + p2.red * w2 + p3.red * w3;
                                sum_g += p0.green * w0 + p1.green * w1 + p2.green * w2 + p3.green * w3;
                                sum_b += p0.blue * w0 + p1.blue * w1 + p2.blue * w2 + p3.blue * w3;
                            }

                            for (; v < window_width; v++)
                            {
                                pixel p = row[v];
                                int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                                long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                                double weight = spatial_weights_flat[sw_base + v] * exp(-(double)cd * inv_two_sigma_r_sq);
                                weight_sum += weight;
                                sum_r += p.red * weight;
                                sum_g += p.green * weight;
                                sum_b += p.blue * weight;
                            }
                        }

                        sw_base += window_width;
                    }

                    double inv_wsum = 1.0 / weight_sum;
                    dst_row[j_in_block].red = (unsigned short)(sum_r * inv_wsum);
                    dst_row[j_in_block].green = (unsigned short)(sum_g * inv_wsum);
                    dst_row[j_in_block].blue = (unsigned short)(sum_b * inv_wsum);
                }
            }
        }
    }

    /* ================= BOUNDARY REGIONS (keep simple) ================= */
    double spatial_weights[2 * r + 1][2 * r + 1];
    for (u = -r; u <= r; u++)
        for (v = -r; v <= r; v++)
            spatial_weights[u + r][v + r] = exp(-(double)(u * u + v * v) / (2.0 * sigma_s * sigma_s));

    /* TOP-LEFT CORNER */
    for (i = 0; i < r; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* TOP EDGE */
    for (i = 0; i < r; i++)
    {
        for (j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* TOP-RIGHT CORNER */
    for (i = 0; i < r; i++)
    {
        for (j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i < 0)
                    src_i = -src_i;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* LEFT EDGE */
    for (i = r; i < block_limit; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* RIGHT EDGE */
    for (i = r; i < block_limit; i++)
    {
        for (j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* BOTTOM-LEFT CORNER */
    for (i = block_limit; i < dim; i++)
    {
        for (j = 0; j < r; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j < 0)
                        src_j = -src_j;

                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* BOTTOM EDGE */
    for (i = block_limit; i < dim; i++)
    {
        for (j = r; j < block_limit; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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

    /* BOTTOM-RIGHT CORNER */
    for (i = block_limit; i < dim; i++)
    {
        for (j = block_limit; j < dim; j++)
        {
            int center_idx = i * dim + j;
            pixel center = src[center_idx];
            int cr = center.red, cg = center.green, cb = center.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0, weight_sum = 0.0;

            for (u = -r; u <= r; u++)
            {
                int src_i = i + u;
                if (src_i >= dim)
                    src_i = 2 * dim - 1 - src_i;

                for (v = -r; v <= r; v++)
                {
                    int src_j = j + v;
                    if (src_j >= dim)
                        src_j = 2 * dim - 1 - src_j;

                    pixel p = src[src_i * dim + src_j];
                    int dr = cr - p.red, dg = cg - p.green, db = cb - p.blue;
                    long long cd = (long long)dr * dr + (long long)dg * dg + (long long)db * db;
                    double weight = spatial_weights[u + r][v + r] * exp(-(double)cd * inv_two_sigma_r_sq);
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
