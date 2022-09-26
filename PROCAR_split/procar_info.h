/*
 * procar_info.h
 *
 *  Created on: 7 juin 2019
 *      Author: gosteau
 */

//#ifndef PROCAR_INFO_H_
//#define PROCAR_INFO_H_
#define N 20000

typedef struct {
	int n_kpt;
	int n_ions;
	int n_bands;
	int n_char;
	int ispin;
	int soc;
	int space;
	char orbs[N];

} INFO_PROCAR;

INFO_PROCAR get_info(char * path);
void print_info(INFO_PROCAR info);

//#endif /* PROCAR_INFO_H_ */
