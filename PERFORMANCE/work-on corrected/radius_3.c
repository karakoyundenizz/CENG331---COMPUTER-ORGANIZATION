// WHOLE FUNCTION ASSUMES RADIUS IS 2

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
static int radius = 3;           /* Filter radius (window size is 2*radius+1) */
static double sigma_s = 2.0;     /* Spatial parameter */
static double sigma_r = 50000.0; /* Range parameter */

#define ACC(pix, i, j)                                                                                                                                                              \
    do                                                                                                                                                                              \
    {                                                                                                                                                                               \
        double calculated_red = current_pixel_red - (pix).red;                                                                                                                      \
        double calculated_green = current_pixel_green - (pix).green;                                                                                                                \
        double calculated_blue = current_pixel_blue - (pix).blue;                                                                                                                   \
        double weight = spatial_weights[i][j] * exp(-(calculated_red * calculated_red + calculated_green * calculated_green + calculated_blue * calculated_blue) / two_sigma_r_sq); \
        weight_sum += weight;                                                                                                                                                       \
        sum_r += (pix).red * weight;                                                                                                                                                \
        sum_g += (pix).green * weight;                                                                                                                                              \
        sum_b += (pix).blue * weight;                                                                                                                                               \
    } while (0)

char fully_unrolled_descr_3[] = "bilateral: unrolled every fucking loop maxing the divided version with radius 3  ";

void fully_unrolled_bilateral_3(int dim, pixel *src, pixel *dst)
{
    int BLOCK = 32;
    int i_blocked, j_blocked, u, v, i, j;
    int i_in_block, j_in_block, block_limit;
    // precompute the spatial weights so that dont need to caclulate for every value in inner iteration
    double spatial_weights[2 * radius + 1][2 * radius + 1];
    double two_sigma_s_sq = 2.0 * sigma_s * sigma_s;
    double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;

    for (u = -radius; u <= radius; u++)
        for (v = -radius; v <= radius; v++)
            spatial_weights[u + radius][v + radius] =
                exp(-(double)(u * u + v * v) / two_sigma_s_sq);

    /*
    WHEN WE REMOVE THE BLOCKS WITH RADIUS DISTINCE TO THE EDGES SO INDEXES DIFFER LIKE
    (RADIUS-(DIMENSION-RADIUS), RADIUS-(DIMENSION-RADIUS))
    SO NO NEED FOR CHECKING EDGES IN EVERY ITERATION
    */
    block_limit = dim - radius;
    for (i_blocked = radius; i_blocked < block_limit; i_blocked += BLOCK)
    {
        for (j_blocked = radius; j_blocked < block_limit; j_blocked += BLOCK)
        {
            int i_end = (i_blocked + BLOCK < block_limit) ? i_blocked + BLOCK : block_limit;
            int j_end = (j_blocked + BLOCK < block_limit) ? j_blocked + BLOCK : block_limit;

            for (i_in_block = i_blocked; i_in_block < i_end; i_in_block++)
            {
                // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=3 SO DO IT WITH UNROLLING FACTOR K=3
                // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
                pixel *row0 = src + (i_in_block - 3) * dim;
                pixel *row1 = src + (i_in_block - 2) * dim;
                pixel *row2 = src + (i_in_block - 1) * dim;
                pixel *row3 = src + i_in_block * dim; // middle pointer so we can clculate just by using it
                pixel *row4 = src + (i_in_block + 1) * dim;
                pixel *row5 = src + (i_in_block + 2) * dim;
                pixel *row6 = src + (i_in_block + 3) * dim;

                for (j_in_block = j_blocked; j_in_block < j_end; j_in_block++)
                {
                    // PICK THE MIDDLE ROW TO CALCULATE OTHERS
                    pixel current_pixel = row3[j_in_block];
                    double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
                    double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
                    double weight_sum = 0.0;
                    pixel *destination = dst + i_in_block * dim + j_in_block;

                    /* fully unrolled 6x6 */
                    ACC(row0[j_in_block - 3], 0, 0);
                    ACC(row0[j_in_block - 2], 0, 1);
                    ACC(row0[j_in_block - 1], 0, 2);
                    ACC(row0[j_in_block], 0, 3);
                    ACC(row0[j_in_block + 1], 0, 4);
                    ACC(row0[j_in_block + 2], 0, 5);
                    ACC(row0[j_in_block + 3], 0, 6);

                    ACC(row1[j_in_block - 3], 1, 0);
                    ACC(row1[j_in_block - 2], 1, 1);
                    ACC(row1[j_in_block - 1], 1, 2);
                    ACC(row1[j_in_block], 1, 3);
                    ACC(row1[j_in_block + 1], 1, 4);
                    ACC(row1[j_in_block + 2], 1, 5);
                    ACC(row1[j_in_block + 3], 1, 6);

                    ACC(row2[j_in_block - 3], 2, 0);
                    ACC(row2[j_in_block - 2], 2, 1);
                    ACC(row2[j_in_block - 1], 2, 2);
                    ACC(row2[j_in_block], 2, 3);
                    ACC(row2[j_in_block + 1], 2, 4);
                    ACC(row2[j_in_block + 2], 2, 5);
                    ACC(row2[j_in_block + 3], 2, 6);

                    ACC(row3[j_in_block - 3], 3, 0);
                    ACC(row3[j_in_block - 2], 3, 1);
                    ACC(row3[j_in_block - 1], 3, 2);
                    ACC(row3[j_in_block], 3, 3);
                    ACC(row3[j_in_block + 1], 3, 4);
                    ACC(row3[j_in_block + 2], 3, 5);
                    ACC(row3[j_in_block + 3], 3, 6);

                    ACC(row4[j_in_block - 3], 4, 0);
                    ACC(row4[j_in_block - 2], 4, 1);
                    ACC(row4[j_in_block - 1], 4, 2);
                    ACC(row4[j_in_block], 4, 3);
                    ACC(row4[j_in_block + 1], 4, 4);
                    ACC(row4[j_in_block + 2], 4, 5);
                    ACC(row4[j_in_block + 3], 4, 6);

                    ACC(row5[j_in_block - 3], 5, 0);
                    ACC(row5[j_in_block - 2], 5, 1);
                    ACC(row5[j_in_block - 1], 5, 2);
                    ACC(row5[j_in_block], 5, 3);
                    ACC(row5[j_in_block + 1], 5, 4);
                    ACC(row5[j_in_block + 2], 5, 5);
                    ACC(row5[j_in_block + 3], 5, 6);

                    ACC(row6[j_in_block - 3], 6, 0);
                    ACC(row6[j_in_block - 2], 6, 1);
                    ACC(row6[j_in_block - 1], 6, 2);
                    ACC(row6[j_in_block], 6, 3);
                    ACC(row6[j_in_block + 1], 6, 4);
                    ACC(row6[j_in_block + 2], 6, 5);
                    ACC(row6[j_in_block + 3], 6, 6);

                    destination->red = (unsigned short)(sum_r / weight_sum);
                    destination->green = (unsigned short)(sum_g / weight_sum);
                    destination->blue = (unsigned short)(sum_b / weight_sum);
                }
            }
        }
    }

    /*
    4 EDGE WITH RADIUS DIFFERENCE RADIUS SO INDEXES DIFFER LIKE
    (0-RADIUS,0-RADIUS) SO NO LOOSE FOR IF CONTROL IN EVERY LOOP
    MAY SPLIT FOR EVERY EDGE LIKE BOTTOM, LEFT, TOP, RIGHT
    SO THAT WE USE NO IF'S JUST FOR LOOPS IN SEQUENCE
    ---------- LOOP INDEX CONTROLS ARE NOT THAT STRONG DEFINITELY CHECK THAT BEFORE SUBMIT ------------------
    */

    /* TRYING TO PROCESS THE TOP LEFT PART the part row < radius and column <radius
    so row index < radius since i+radi
    */

    for (i = 0; i < radius; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=3 SO DO IT WITH UNROLLING FACTOR K=3
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (- i + 3) * dim;
        pixel *row1 = src + (-i +2) * dim;
        pixel *row2;
        // only for 2 it is positive so no need for changing sign but later maybe ı can implement directly without an if
        // since it only need to go to else once for 2 others need to go if part so since always taken it makes sense to put here
        if(i != 2){
            row2 = src + (-i + 1) * dim;
        }
        else{
            row2 = src + (i - 1) * dim;
        }
        pixel *row3 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row4 = src + (i + 1) * dim;
        pixel *row5 = src + (i + 2) * dim;
        pixel *row6 = src + (i + 3) * dim;

        for (j = 0; j < radius; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row3[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 6x6 */
                    ACC(row0[j - 3], 0, 0);
                    ACC(row0[j - 2], 0, 1);
                    ACC(row0[j - 1], 0, 2);
                    ACC(row0[j], 0, 3);
                    ACC(row0[j + 1], 0, 4);
                    ACC(row0[j + 2], 0, 5);
                    ACC(row0[j + 3], 0, 6);

                    ACC(row1[j - 3], 1, 0);
                    ACC(row1[j - 2], 1, 1);
                    ACC(row1[j - 1], 1, 2);
                    ACC(row1[j], 1, 3);
                    ACC(row1[j + 1], 1, 4);
                    ACC(row1[j + 2], 1, 5);
                    ACC(row1[j + 3], 1, 6);

                    ACC(row2[j - 3], 2, 0);
                    ACC(row2[j - 2], 2, 1);
                    ACC(row2[j - 1], 2, 2);
                    ACC(row2[j], 2, 3);
                    ACC(row2[j + 1], 2, 4);
                    ACC(row2[j + 2], 2, 5);
                    ACC(row2[j + 3], 2, 6);

                    ACC(row3[j - 3], 3, 0);
                    ACC(row3[j - 2], 3, 1);
                    ACC(row3[j - 1], 3, 2);
                    ACC(row3[j], 3, 3);
                    ACC(row3[j + 1], 3, 4);
                    ACC(row3[j + 2], 3, 5);
                    ACC(row3[j + 3], 3, 6);

                    ACC(row4[j - 3], 4, 0);
                    ACC(row4[j - 2], 4, 1);
                    ACC(row4[j - 1], 4, 2);
                    ACC(row4[j], 4, 3);
                    ACC(row4[j + 1], 4, 4);
                    ACC(row4[j + 2], 4, 5);
                    ACC(row4[j + 3], 4, 6);

                    ACC(row5[j - 3], 5, 0);
                    ACC(row5[j - 2], 5, 1);
                    ACC(row5[j - 1], 5, 2);
                    ACC(row5[j], 5, 3);
                    ACC(row5[j + 1], 5, 4);
                    ACC(row5[j + 2], 5, 5);
                    ACC(row5[j + 3], 5, 6);

                    ACC(row6[j - 3], 6, 0);
                    ACC(row6[j - 2], 6, 1);
                    ACC(row6[j - 1], 6, 2);
                    ACC(row6[j], 6, 3);
                    ACC(row6[j + 1], 6, 4);
                    ACC(row6[j + 2], 6, 5);
                    ACC(row6[j + 3], 6, 6);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE TOP MIDDLE PART the part row < radius and column >radius and column < (dimension-radius)
    so row index < radius since i+radi
    */

    for (i = 0; i < radius; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (-i + 2) * dim;
        pixel *row1 = src + (-i + 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 2) * dim;

        for (j = radius; j < block_limit; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 2], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 2], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 2], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 2], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 2], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE TOP RIGHT PART the part row < radius and dimension > column > (dimension-radius)
    so row index < radius since i+radi
    */
    for (i = 0; i < radius; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (-i + 2) * dim;
        pixel *row1 = src + (-i + 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 2) * dim;

        for (j = block_limit; j < dim - 1; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }

        // to calculate the second column where both +1 and +2 exceeds the dimension
        for (j = dim - 1; j < dim; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j], 0, 3);
            ACC(row0[j - 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j], 1, 3);
            ACC(row1[j - 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j], 2, 3);
            ACC(row2[j - 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j], 3, 3);
            ACC(row3[j - 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j], 4, 3);
            ACC(row4[j - 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE LEFTMOST MIDDLE PART the part row > radius and column < radius
    so row index < radius since i+radi
    */
    for (i = radius; i < block_limit; i++)
    {
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 2) * dim;

        for (j = 0; j < radius; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[-j + 2], 0, 0);
            ACC(row0[-j + 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 2], 0, 4);

            ACC(row1[-j + 2], 1, 0);
            ACC(row1[-j + 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 2], 1, 4);

            ACC(row2[-j + 2], 2, 0);
            ACC(row2[-j + 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 2], 2, 4);

            ACC(row3[-j + 2], 3, 0);
            ACC(row3[-j + 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 2], 3, 4);

            ACC(row4[-j + 2], 4, 0);
            ACC(row4[-j + 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 2], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE RIGHTMOST MIDDLE PART the part (dimension-radius) > row > radius and column > (dimension-radius)
    so row index < radius since i+radi
    */
    for (i = radius; i < block_limit; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 2) * dim;

        for (j = block_limit; j < dim - 1; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }

        // to calculate the second column where both +1 and +2 exceeds the dimension
        for (j = dim - 1; j < dim; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j], 0, 3);
            ACC(row0[j - 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j], 1, 3);
            ACC(row1[j - 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j], 2, 3);
            ACC(row2[j - 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j], 3, 3);
            ACC(row3[j - 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j], 4, 3);
            ACC(row4[j - 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE BOTTOM LEFT PART the part dimension > row > (dimension-radius) and column <radius
        so row index < radius since i+radi
        */
    for (i = block_limit; i < dim - 1; i++)
    {
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 1) * dim;

        for (j = 0; j < radius; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[-j + 2], 0, 0);
            ACC(row0[-j + 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 2], 0, 4);

            ACC(row1[-j + 2], 1, 0);
            ACC(row1[-j + 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 2], 1, 4);

            ACC(row2[-j + 2], 2, 0);
            ACC(row2[-j + 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 2], 2, 4);

            ACC(row3[-j + 2], 3, 0);
            ACC(row3[-j + 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 2], 3, 4);

            ACC(row4[-j + 2], 4, 0);
            ACC(row4[-j + 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 2], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    for (i = dim - 1; i < dim; i++)
    {
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i)*dim;
        pixel *row4 = src + (i - 1) * dim;

        for (j = 0; j < radius; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[-j + 2], 0, 0);
            ACC(row0[-j + 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 2], 0, 4);

            ACC(row1[-j + 2], 1, 0);
            ACC(row1[-j + 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 2], 1, 4);

            ACC(row2[-j + 2], 2, 0);
            ACC(row2[-j + 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 2], 2, 4);

            ACC(row3[-j + 2], 3, 0);
            ACC(row3[-j + 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 2], 3, 4);

            ACC(row4[-j + 2], 4, 0);
            ACC(row4[-j + 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 2], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE BOTTOM MIDDLE PART the part row < radius and column >radius and column < (dimension-radius)
    so row index < radius since i+radi
    */

    for (i = block_limit; i < dim - 1; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 1) * dim;

        for (j = radius; j < block_limit; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 2], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 2], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 2], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 2], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 2], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    for (i = dim - 1; i < dim; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i)*dim;
        pixel *row4 = src + (i - 1) * dim;

        for (j = radius; j < block_limit; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 2], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 2], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 2], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 2], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 2], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /* TRYING TO PROCESS THE BOTTOM RIGHT PART the part row < radius and dimension > column > (dimension-radius)
    so row index < radius since i+radi
    */
    for (i = block_limit; i < dim - 1; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i + 1) * dim;
        pixel *row4 = src + (i + 1) * dim;

        for (j = block_limit; j < dim - 1; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }

        // to calculate the second column where both +1 and +2 exceeds the dimension
        for (j = dim - 1; j < dim; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j], 0, 3);
            ACC(row0[j - 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j], 1, 3);
            ACC(row1[j - 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j], 2, 3);
            ACC(row2[j - 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j], 3, 3);
            ACC(row3[j - 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j], 4, 3);
            ACC(row4[j - 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }

    /*
    THERE IS 2 POSSIBLITY FIRST ROW EXCEEDS THE LIMIT ONLY ON +2 BUT SECOND ROW ALSO EXCEDDS ON +1,+2 SO HANDLE THAT SEPERATELY
    BOTTOM RIGHT AT DIMENSION - 1 TH ROW
    */

    for (i = dim - 1; i < dim; i++)
    {
        // can use block limit since it is same (dimension-radius)

        // DOING UNROLLING BUT IT ASSUMES THAT RADIUS=2 SO DO IT WITH UNROLLING FACTOR K=2
        // since it has 2 multilply unit and latency 4 it can process 8 multiplication at once
        pixel *row0 = src + (i - 2) * dim;
        pixel *row1 = src + (i - 1) * dim;
        pixel *row2 = src + i * dim; // middle pointer so we can clculate just by using it
        pixel *row3 = src + (i)*dim;
        pixel *row4 = src + (i - 1) * dim;

        // COLUMNS ALSO SURPASS THE DIMENSIONS SO HANDLE 2 CASES WHERE EACH
        for (j = block_limit; j < dim - 1; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j + 1], 0, 3);
            ACC(row0[j + 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j + 1], 1, 3);
            ACC(row1[j + 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j + 1], 2, 3);
            ACC(row2[j + 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j + 1], 3, 3);
            ACC(row3[j + 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j + 1], 4, 3);
            ACC(row4[j + 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }

        // to calculate the second column where both +1 and +2 exceeds the dimension
        for (j = dim - 1; j < dim; j++)
        {
            // PICK THE MIDDLE ROW TO CALCULATE OTHERS
            pixel current_pixel = row2[j];
            double current_pixel_red = current_pixel.red, current_pixel_green = current_pixel.green, current_pixel_blue = current_pixel.blue;
            double sum_r = 0.0, sum_g = 0.0, sum_b = 0.0;
            double weight_sum = 0.0;
            pixel *destination = dst + i * dim + j;

            /* fully unrolled 5x5 */
            ACC(row0[j - 2], 0, 0);
            ACC(row0[j - 1], 0, 1);
            ACC(row0[j], 0, 2);
            ACC(row0[j], 0, 3);
            ACC(row0[j - 1], 0, 4);

            ACC(row1[j - 2], 1, 0);
            ACC(row1[j - 1], 1, 1);
            ACC(row1[j], 1, 2);
            ACC(row1[j], 1, 3);
            ACC(row1[j - 1], 1, 4);

            ACC(row2[j - 2], 2, 0);
            ACC(row2[j - 1], 2, 1);
            ACC(row2[j], 2, 2);
            ACC(row2[j], 2, 3);
            ACC(row2[j - 1], 2, 4);

            ACC(row3[j - 2], 3, 0);
            ACC(row3[j - 1], 3, 1);
            ACC(row3[j], 3, 2);
            ACC(row3[j], 3, 3);
            ACC(row3[j - 1], 3, 4);

            ACC(row4[j - 2], 4, 0);
            ACC(row4[j - 1], 4, 1);
            ACC(row4[j], 4, 2);
            ACC(row4[j], 4, 3);
            ACC(row4[j - 1], 4, 4);

            destination->red = (unsigned short)(sum_r / weight_sum);
            destination->green = (unsigned short)(sum_g / weight_sum);
            destination->blue = (unsigned short)(sum_b / weight_sum);
        }
    }
}