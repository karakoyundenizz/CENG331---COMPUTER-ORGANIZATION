/*
 * driver.h - Various definitions for the Performance Lab.
 * 
 * DO NOT MODIFY ANYTHING IN THIS FILE
 */
#ifndef _DEFS_H_
#define _DEFS_H_

#include <stdlib.h>

#define RIDX(i,j,n) ((i)*(n)+(j))

typedef struct {
  char *team;
  char *name1, *email1;
  char *name2, *email2;
  char *name3, *email3;
} team_t;

extern team_t team;

typedef struct {
   unsigned short red;
   unsigned short green;
   unsigned short blue;
} pixel;

typedef void (*lab_test_func) (int, pixel*, pixel*);

void blur(int, pixel *, pixel *);
void bilateral(int, pixel *, pixel *);

void register_bilateral_functions(void);
void register_blur_functions(void);
void add_blur_function(lab_test_func, char*);
void add_bilateral_function(lab_test_func, char*);

#endif /* _DEFS_H_ */

