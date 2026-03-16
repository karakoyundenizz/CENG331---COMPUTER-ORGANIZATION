/*******************************************************************
 * 
 * driver.c - Driver program for CS:APP Performance Lab
 * 
 * In kernels.c, students generate an arbitrary number of bilateral and
 * invert test functions, which they then register with the driver
 * program using the add_bilateral_function() and add_invert_function()
 * functions.
 * 
 * The driver program runs and measures the registered test functions
 * and reports their performance.
 * 
 * Copyright (c) 2002, R. Bryant and D. O'Hallaron, All rights
 * reserved.  May not be used, modified, or copied without permission.
 *
 ********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <math.h>

/* Platform-specific includes */
#if defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
    #include <io.h>
    #define getopt(argc, argv, optstring) _getopt(argc, argv, optstring)
    extern char *optarg;
    extern int optind;
#else
    #include <unistd.h>
    #include <sys/time.h>
#endif
#include "fcyc.h"
#include "defs.h"
#include "config.h"

/* Team structure that identifies the students */
extern team_t team; 

/* Keep track of a number of different test functions */
#define MAX_BENCHMARKS 100
#define DIM_CNT 5

/* Misc constants */
#define BSIZE 32     /* cache block size in bytes */     
#define MAX_DIM 1280 /* 1024 + 256 */
#define ODD_DIM 96   /* not a power of 2 */

/* fast versions of min and max */
#define min(a,b) (a < b ? a : b)
#define max(a,b) (a > b ? a : b)

/* This struct characterizes the results for one benchmark test */
typedef struct {
    lab_test_func tfunct; /* The test function */
    double cpes[DIM_CNT]; /* One CPE result for each dimension */
    char *description;    /* ASCII description of the test function */
    unsigned short valid; /* The function is tested if this is non zero */
} bench_t;

/* The range of image dimensions that we will be testing */
static int test_dim_bilateral[] = {64, 128, 256, 512, 1024};
static int test_dim_blur[] = {64, 128, 256, 512, 1024};

/* Baseline CPEs (see config.h) */
static double bilateral_baseline_cpes[] = {B64, B128, B256, B512, B1024};
static double blur_baseline_cpes[] = {BL64, BL128, BL256, BL512, BL1024};

/* These hold the results for all benchmarks */
static bench_t benchmarks_bilateral[MAX_BENCHMARKS];
static bench_t benchmarks_blur[MAX_BENCHMARKS];

/* These give the sizes of the above lists */
static int bilateral_benchmark_count = 0;
static int blur_benchmark_count = 0;

/* 
 * An image is a dimxdim matrix of pixels stored in a 1D array.  The
 * data array holds three images (the input original, a copy of the original, 
 * and the output result array. There is also an additional BSIZE bytes
 * of padding for alignment to cache block boundaries.
 */
static pixel data[(3*MAX_DIM*MAX_DIM) + (BSIZE/sizeof(pixel))];

/* Various image pointers */
static pixel *orig = NULL;         /* original image */
static pixel *copy_of_orig = NULL; /* copy of original for checking result */
static pixel *result = NULL;       /* result image */

/* Keep track of the best bilateral and blur score for grading */
double bilateral_maxmean = 0.0;
char *bilateral_maxmean_desc = NULL;

double blur_maxmean = 0.0;
char *blur_maxmean_desc = NULL;


/******************** Functions begin *************************/

void add_blur_function(lab_test_func f, char *description) 
{
    benchmarks_blur[blur_benchmark_count].tfunct = f;
    benchmarks_blur[blur_benchmark_count].description = description;
    benchmarks_blur[blur_benchmark_count].valid = 0;  
    blur_benchmark_count++;
}


void add_bilateral_function(lab_test_func f, char *description) 
{
    benchmarks_bilateral[bilateral_benchmark_count].tfunct = f;
    benchmarks_bilateral[bilateral_benchmark_count].description = description;
    benchmarks_bilateral[bilateral_benchmark_count].valid = 0;
    bilateral_benchmark_count++;
}

/* 
 * random_in_interval - Returns random integer in interval [low, high) 
 */
static int random_in_interval(int low, int high) 
{
    int size = high - low;
    return (rand()% size) + low;
}

/*
 * create - creates a dimxdim image aligned to a BSIZE byte boundary
 */
static void create(int dim)
{
    int i, j;

    /* Align the images to BSIZE byte boundaries */
    orig = data;
    while ((unsigned long long)orig % BSIZE) {
        char *tmp = (char *) orig;
        tmp++;
        orig = (pixel *) tmp;
    }
    result = orig + dim*dim;
    copy_of_orig = result + dim*dim;

    for (i = 0; i < dim; i++) {
		for (j = 0; j < dim; j++) {
			/* Original image initialized to random colors */
			orig[RIDX(i,j,dim)].red = random_in_interval(0, 65536);
			orig[RIDX(i,j,dim)].green = random_in_interval(0, 65536);
			orig[RIDX(i,j,dim)].blue = random_in_interval(0, 65536);

			/* Copy of original image for checking result */
			copy_of_orig[RIDX(i,j,dim)].red = orig[RIDX(i,j,dim)].red;
			copy_of_orig[RIDX(i,j,dim)].green = orig[RIDX(i,j,dim)].green;
			copy_of_orig[RIDX(i,j,dim)].blue = orig[RIDX(i,j,dim)].blue;

			/* Result image initialized to all black */
			result[RIDX(i,j,dim)].red = 0;
			result[RIDX(i,j,dim)].green = 0;
			result[RIDX(i,j,dim)].blue = 0;
		}
    }

    return;
}


/* 
 * compare_pixels - Returns 1 if the two arguments don't have same RGB
 *    values, 0 o.w.  
 */
static int compare_pixels(pixel p1, pixel p2) 
{
    return 
	(p1.red != p2.red) || 
	(p1.green != p2.green) || 
	(p1.blue != p2.blue);
}


/* Make sure the orig array is unchanged */
static int check_orig(int dim) 
{
    int i, j;

    for (i = 0; i < dim; i++) 
	for (j = 0; j < dim; j++) 
	    if (compare_pixels(orig[RIDX(i,j,dim)], copy_of_orig[RIDX(i,j,dim)])) {
		printf("\n");
		printf("Error: Original image has been changed!\n");
		return 1;
	    }

    return 0;
}

/* 
 * check_bilateral - Make sure the bilateral filter actually works. 
 * The orig array should not have been tampered with! 
 */
static int check_bilateral(int dim) 
{
    int err = 0;
    int i, j;
    int badi = 0;
    int badj = 0;
    pixel expected, actual;
    int u, v;
    int src_i, src_j;
    double sum_r, sum_g, sum_b;
    double weight_sum;
    double spatial_weight, range_weight, weight;
    double dist_sq, color_diff_sq;
    int radius = 2;
    double sigma_s = 2.0;
    double sigma_r = 50000.0;  /* Updated to match implementation */
    double two_sigma_s_sq = 2.0 * sigma_s * sigma_s;
    double two_sigma_r_sq = 2.0 * sigma_r * sigma_r;

    /* return 1 if the original image has been changed */
    if (check_orig(dim)) 
	return 1; 

    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++) {
	    sum_r = sum_g = sum_b = 0.0;
	    weight_sum = 0.0;
	    
	    for (u = -radius; u <= radius; u++) {
		for (v = -radius; v <= radius; v++) {
		    src_i = i + u;
		    src_j = j + v;
		    
		    /* Mirror padding at boundaries */
		    if(src_i < 0) src_i = -src_i;
		    if(src_i >= dim) src_i = 2*dim - 1 - src_i;
		    if(src_j < 0) src_j = -src_j;
		    if(src_j >= dim) src_j = 2*dim - 1 - src_j;
		    
		    /* Spatial weight */
		    dist_sq = u*u + v*v;
		    spatial_weight = exp(-dist_sq / two_sigma_s_sq);
		    
		    /* Range weight (using L2 norm of RGB difference) */
		    /* Cast to double to prevent unsigned underflow */
		    color_diff_sq = ((double)orig[RIDX(i, j, dim)].red - 
				     (double)orig[RIDX(src_i, src_j, dim)].red) *
				    ((double)orig[RIDX(i, j, dim)].red - 
				     (double)orig[RIDX(src_i, src_j, dim)].red) +
				    ((double)orig[RIDX(i, j, dim)].green - 
				     (double)orig[RIDX(src_i, src_j, dim)].green) *
				    ((double)orig[RIDX(i, j, dim)].green - 
				     (double)orig[RIDX(src_i, src_j, dim)].green) +
				    ((double)orig[RIDX(i, j, dim)].blue - 
				     (double)orig[RIDX(src_i, src_j, dim)].blue) *
				    ((double)orig[RIDX(i, j, dim)].blue - 
				     (double)orig[RIDX(src_i, src_j, dim)].blue);
		    range_weight = exp(-color_diff_sq / two_sigma_r_sq);
		    
		    weight = spatial_weight * range_weight;
		    weight_sum += weight;
		    
		    sum_r += orig[RIDX(src_i, src_j, dim)].red * weight;
		    sum_g += orig[RIDX(src_i, src_j, dim)].green * weight;
		    sum_b += orig[RIDX(src_i, src_j, dim)].blue * weight;
		}
	    }
	    
	    expected.red = (unsigned short)(sum_r / weight_sum);
	    expected.green = (unsigned short)(sum_g / weight_sum);
	    expected.blue = (unsigned short)(sum_b / weight_sum);
	    
	    if (compare_pixels(result[RIDX(i, j, dim)], expected)) {
		err++;
		badi = i;
		badj = j;
		actual = result[RIDX(i, j, dim)];
	    }
	}
    }

    if (err) {
	printf("\n");
	printf("ERROR: Dimension=%d, %d errors\n", dim, err);    
	printf("E.g., The following two pixels should have equal value:\n");
	printf("Expected[%d][%d].{red,green,blue} = {%d,%d,%d}\n",
	       badi, badj, expected.red, expected.green, expected.blue);
	printf("Got[%d][%d].{red,green,blue} = {%d,%d,%d}\n",
	       badi, badj, actual.red, actual.green, actual.blue);
    }

    return err;
}

/*
 * check_blur - Make sure the box blur function actually works.  The
 * orig array should not have been tampered with!  
 */
static int check_blur(int dim) {
    int err = 0;
    int i, j, u, v;
    int src_i, src_j;
    int badi = 0;
    int badj = 0;
    pixel expected, actual;
    int radius = 2;
    int window_size = 2 * radius + 1;
    int sum_r, sum_g, sum_b;

    /* return 1 if original image has been changed */
    if (check_orig(dim)) 
	return 1; 

    for (i = 0; i < dim; i++) {
	for (j = 0; j < dim; j++) {
	    sum_r = sum_g = sum_b = 0;
	    
	    for (u = -radius; u <= radius; u++) {
		for (v = -radius; v <= radius; v++) {
		    src_i = i + u;
		    src_j = j + v;
		    
		    /* Mirror padding at boundaries */
		    if(src_i < 0) src_i = -src_i;
		    if(src_i >= dim) src_i = 2*dim - 1 - src_i;
		    if(src_j < 0) src_j = -src_j;
		    if(src_j >= dim) src_j = 2*dim - 1 - src_j;
		    
		    sum_r += orig[RIDX(src_i, src_j, dim)].red;
		    sum_g += orig[RIDX(src_i, src_j, dim)].green;
		    sum_b += orig[RIDX(src_i, src_j, dim)].blue;
		}
	    }
	    
	    expected.red = (unsigned short)(sum_r / (window_size * window_size));
	    expected.green = (unsigned short)(sum_g / (window_size * window_size));
	    expected.blue = (unsigned short)(sum_b / (window_size * window_size));
	    
	    if (compare_pixels(result[RIDX(i,j,dim)], expected)) {
		err++;
		badi = i;
		badj = j;
		actual = result[RIDX(i,j,dim)];
	    }
	}
    }

    if (err) {
	printf("\n");
	printf("ERROR: Dimension=%d, %d errors\n", dim, err);    
	printf("E.g., \n");
	printf("You have dst[%d][%d].{red,green,blue} = {%d,%d,%d}\n",
	       badi, badj, actual.red, actual.green, actual.blue);
	printf("It should be dst[%d][%d].{red,green,blue} = {%d,%d,%d}\n",
	       badi, badj, expected.red, expected.green, expected.blue);
    }

    return err;
}


void func_wrapper(void *arglist[]) 
{
    pixel *src, *dst;
    int mydim;
    lab_test_func f;

    f = (lab_test_func) arglist[0];
    mydim = *((int *) arglist[1]);
    src = (pixel *) arglist[2];
    dst = (pixel *) arglist[3];

    (*f)(mydim, src, dst);

    return;
}

void run_bilateral_benchmark(int idx, int dim) 
{
    benchmarks_bilateral[idx].tfunct(dim, orig, result);
}

void test_bilateral(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_bilateral[bench_index].description;
  
    for (test_num = 0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check for odd dimension */
	create(ODD_DIM);
	run_bilateral_benchmark(bench_index, ODD_DIM);
	if (check_bilateral(ODD_DIM)) {
	    printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		   benchmarks_bilateral[bench_index].description, ODD_DIM);
	    return;
	}

	/* Create a test image of the required dimension */
	dim = test_dim_bilateral[test_num];
	create(dim);
#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_bilateral[bench_index].description);
#endif

	/* Check that the code works */
	run_bilateral_benchmark(bench_index, dim);
	if (check_bilateral(dim)) {
	    printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		   benchmarks_bilateral[bench_index].description, dim);
	    return;
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_bilateral[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;

	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_bilateral[bench_index].cpes[test_num] = cpe;
	}
    }

    /* 
     * Print results as a table 
     */
    printf("Bilateral: Version = %s:\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_bilateral[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_bilateral[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", bilateral_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute Speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_bilateral[bench_index].cpes[i] > 0.0) {
		ratio = bilateral_baseline_cpes[i]/
		    benchmarks_bilateral[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}

	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > bilateral_maxmean) {
	    bilateral_maxmean = mean;
	    bilateral_maxmean_desc = benchmarks_bilateral[bench_index].description;
	}
    }


#ifdef DEBUG
    fflush(stdout);
#endif
    return;  
}

void run_blur_benchmark(int idx, int dim) 
{
    benchmarks_blur[idx].tfunct(dim, orig, result);
}

void test_blur(int bench_index) 
{
    int i;
    int test_num;
    char *description = benchmarks_blur[bench_index].description;
  
    for(test_num=0; test_num < DIM_CNT; test_num++) {
	int dim;

	/* Check correctness for odd (non power of two dimensions */
	create(ODD_DIM);
	run_blur_benchmark(bench_index, ODD_DIM);
	if (check_blur(ODD_DIM)) {
	    printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		   benchmarks_blur[bench_index].description, ODD_DIM);
	    return;
	}

	/* Create a test image of the required dimension */
	dim = test_dim_blur[test_num];
	create(dim);

#ifdef DEBUG
	printf("DEBUG: Running benchmark \"%s\"\n", benchmarks_blur[bench_index].description);
#endif
	/* Check that the code works */
	run_blur_benchmark(bench_index, dim);
	if (check_blur(dim)) {
	    printf("Benchmark \"%s\" failed correctness check for dimension %d.\n",
		   benchmarks_blur[bench_index].description, dim);
	    return;
	}

	/* Measure CPE */
	{
	    double num_cycles, cpe;
	    int tmpdim = dim;
	    void *arglist[4];
	    double dimension = (double) dim;
	    double work = dimension*dimension;
#ifdef DEBUG
	    printf("DEBUG: dimension=%.1f\n",dimension);
	    printf("DEBUG: work=%.1f\n",work);
#endif
	    arglist[0] = (void *) benchmarks_blur[bench_index].tfunct;
	    arglist[1] = (void *) &tmpdim;
	    arglist[2] = (void *) orig;
	    arglist[3] = (void *) result;
        
	    create(dim);
	    num_cycles = fcyc_v((test_funct_v)&func_wrapper, arglist); 
	    cpe = num_cycles/work;
	    benchmarks_blur[bench_index].cpes[test_num] = cpe;
	}
    }

    /* Print results as a table */
    printf("Blur: Version = %s:\n", description);
    printf("Dim\t");
    for (i = 0; i < DIM_CNT; i++)
	printf("\t%d", test_dim_blur[i]);
    printf("\tMean\n");
  
    printf("Your CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", benchmarks_blur[bench_index].cpes[i]);
    }
    printf("\n");

    printf("Baseline CPEs");
    for (i = 0; i < DIM_CNT; i++) {
	printf("\t%.1f", blur_baseline_cpes[i]);
    }
    printf("\n");

    /* Compute speedup */
    {
	double prod, ratio, mean;
	prod = 1.0; /* Geometric mean */
	printf("Speedup\t");
	for (i = 0; i < DIM_CNT; i++) {
	    if (benchmarks_blur[bench_index].cpes[i] > 0.0) {
		ratio = blur_baseline_cpes[i]/
		    benchmarks_blur[bench_index].cpes[i];
	    }
	    else {
		printf("Fatal Error: Non-positive CPE value...\n");
		exit(EXIT_FAILURE);
	    }
	    prod *= ratio;
	    printf("\t%.1f", ratio);
	}
	/* Geometric mean */
	mean = pow(prod, 1.0/(double) DIM_CNT);
	printf("\t%.1f", mean);
	printf("\n\n");
	if (mean > blur_maxmean) {
	    blur_maxmean = mean;
	    blur_maxmean_desc = benchmarks_blur[bench_index].description;
	}
    }

    return;  
}


void usage(char *progname) 
{
    fprintf(stderr, "Usage: %s [-hqg] [-f <func_file>] [-d <dump_file>]\n", progname);    
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -h         Print this message\n");
    fprintf(stderr, "  -q         Quit after dumping (use with -d )\n");
    fprintf(stderr, "  -g         Autograder mode: checks only bilateral() and invert()\n");
    fprintf(stderr, "  -f <file>  Get test function names from dump file <file>\n");
    fprintf(stderr, "  -d <file>  Emit a dump file <file> for later use with -f\n");
    exit(EXIT_FAILURE);
}



int main(int argc, char *argv[])
{
    int i;
    int quit_after_dump = 0;
    int skip_teamname_check = 0;
    int autograder = 0;
    int seed = 1729;
    char c = '0';
    char *bench_func_file = NULL;
    char *func_dump_file = NULL;

    /* register all the defined functions */
    register_bilateral_functions();
    register_blur_functions();

    /* parse command line args */
    while ((c = getopt(argc, argv, "tgqf:d:s:h")) != -1)
	switch (c) {

	case 't': /* skip team name check (hidden flag) */
	    skip_teamname_check = 1;
	    break;

	case 's': /* seed for random number generator (hidden flag) */
	    seed = atoi(optarg);
	    break;

	case 'g': /* autograder mode (checks only bilateral() and blur()) */
	    autograder = 1;
	    break;

	case 'q':
	    quit_after_dump = 1;
	    break;

	case 'f': /* get names of benchmark functions from this file */
	    bench_func_file = strdup(optarg);
	    break;

	case 'd': /* dump names of benchmark functions to this file */
	    func_dump_file = strdup(optarg);
	    {
		int i;
		FILE *fp = fopen(func_dump_file, "w");	

		if (fp == NULL) {
		    printf("Can't open file %s\n",func_dump_file);
		    exit(-5);
		}

		for(i = 0; i < bilateral_benchmark_count; i++) {
		    fprintf(fp, "B:%s\n", benchmarks_bilateral[i].description); 
		}
		for(i = 0; i < blur_benchmark_count; i++) {
		    fprintf(fp, "L:%s\n", benchmarks_blur[i].description); 
		}
		fclose(fp);
	    }
	    break;

	case 'h': /* print help message */
	    usage(argv[0]);

	default: /* unrecognized argument */
	    usage(argv[0]);
	}

    if (quit_after_dump) 
	exit(EXIT_SUCCESS);


    /* Print team info */
    if (!skip_teamname_check) {
	if (strcmp("bovik", team.team) == 0) {
	    printf("%s: Please fill in the team struct in kernels.c.\n", argv[0]);
	    exit(1);
	}
	printf("Teamname: %s\n", team.team);
	printf("Member 1: %s\n", team.name1);
	printf("Email 1: %s\n", team.email1);
	if ((team.name2 && *team.name2) || (team.email2 && *team.email2)) {
	    printf("Member 2: %s\n", team.name2);
	    printf("Email 2: %s\n", team.email2);
	}
	if ((team.name3 && *team.name3) || (team.email3 && *team.email3)) {
	    printf("Member 3: %s\n", team.name3);
	    printf("Email 3: %s\n", team.email3);
	}
        printf("\n");
    }

    srand(seed);

    /* 
     * If we are running in autograder mode, we will only test
     * the bilateral() and blur() functions.
     */
    if (autograder) {
	bilateral_benchmark_count = 1;
	blur_benchmark_count = 1;

	benchmarks_bilateral[0].tfunct = bilateral;
	benchmarks_bilateral[0].description = "bilateral() function";
	benchmarks_bilateral[0].valid = 1;

	benchmarks_blur[0].tfunct = blur;
	benchmarks_blur[0].description = "blur() function";
	benchmarks_blur[0].valid = 1;
    }

    /* 
     * If the user specified a file name using -f, then use
     * the file to determine the versions of bilateral and blur to test
     */
    else if (bench_func_file != NULL) {
	char flag;
	char func_line[256];
	FILE *fp = fopen(bench_func_file, "r");

	if (fp == NULL) {
	    printf("Can't open file %s\n",bench_func_file);
	    exit(-5);
	}
    
	while(func_line == fgets(func_line, 256, fp)) {
	    char *func_name = func_line;
	    char **strptr = &func_name;
	    char *token = strsep(strptr, ":");
	    flag = token[0];
	    func_name = strsep(strptr, "\n");
#ifdef DEBUG
	    printf("Function Description is %s\n",func_name);
#endif

	    if (flag == 'B') {
		for(i=0; i<bilateral_benchmark_count; i++) {
		    if (strcmp(benchmarks_bilateral[i].description, func_name) == 0)
			benchmarks_bilateral[i].valid = 1;
		}
	    }
	    else if (flag == 'L') {
		for(i=0; i<blur_benchmark_count; i++) {
		    if (strcmp(benchmarks_blur[i].description, func_name) == 0)
			benchmarks_blur[i].valid = 1;
		}
	    }      
	}

	fclose(fp);
    }

    /* 
     * If the user didn't specify a dump file using -f, then 
     * test all of the functions
     */
    else { /* set all valid flags to 1 */
	for (i = 0; i < bilateral_benchmark_count; i++)
	    benchmarks_bilateral[i].valid = 1;
	for (i = 0; i < blur_benchmark_count; i++)
	    benchmarks_blur[i].valid = 1;
    }

    /* Set measurement (fcyc) parameters */
    set_fcyc_cache_size(1 << 14); /* 16 KB cache size */
    set_fcyc_clear_cache(1); /* clear the cache before each measurement */
    set_fcyc_compensate(1); /* try to compensate for timer overhead */
 
    for (i = 0; i < bilateral_benchmark_count; i++) {
	if (benchmarks_bilateral[i].valid)
	    test_bilateral(i);
    
}
    for (i = 0; i < blur_benchmark_count; i++) {
	if (benchmarks_blur[i].valid)
	    test_blur(i);
    }


    if (autograder) {
	printf("\nbestscores:%.1f:%.1f:\n", bilateral_maxmean, blur_maxmean);
    }
    else {
	printf("Summary of Your Best Scores:\n");
	printf("  Bilateral: %3.1f (%s)\n", bilateral_maxmean, bilateral_maxmean_desc);
	printf("  Blur: %3.1f (%s)\n", blur_maxmean, blur_maxmean_desc);
    }

    return 0;
}














