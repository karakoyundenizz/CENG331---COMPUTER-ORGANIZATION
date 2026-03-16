/*********************************************************
 * config.h - Configuration data for the driver.c program.
 *********************************************************/
#ifndef _CONFIG_H_
#define _CONFIG_H_

/*
 * CPEs for the baseline (naive) version of the bilateral filter function that
 * was handed out to the students. Bd is the measured CPE for a dxd
 * image. Run the driver.c program on your system to get these
 * numbers.
 */
#define B64    1442.9
#define B128   1440.9
#define B256   1447.5
#define B512   1448.6
#define B1024  1451.9

/*
 * CPEs for the baseline (naive) version of the box blur function that
 * was handed out to the students. Bd is the measured CPE for a dxd
 * image. Run the driver.c program on your system to get these
 * numbers.
 */
#define BL64    103.6
#define BL128   103.5
#define BL256   103.5
#define BL512   103.8
#define BL1024  104.5


#endif /* _CONFIG_H_ */
