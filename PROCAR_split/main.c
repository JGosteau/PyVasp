/*
 * main.c
 *
 *  Created on: 7 juin 2019
 *      Author: gosteau
 */
#include"procar_info.h"
#include"split.h"
#include<stdio.h>
/*
int main(){
	INFO_PROCAR info;
	info.n_kpt = 8000;
	info.n_ions = 5;
	info.n_bands = 200;
	info.ispin = 1;
	info.soc = 1;
	info.space = 0;
	//print_info(info);
	//printf("%s\n", argv[1]);
	//char path[200] = "/home/gosteau/Documents/These/1ere année/STO/OUTPUT_dir_ncl/3.9/BS_SOC/";
	char path[200] = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/constrain/a375_c454/BS/";
	//snprintf(path, 200,"%s", argv[1]);
	printf("%s\n", path);
	info = get_info(path);
	print_info(info);
	split_procar(path, info);
	printf("Done\n");
}
*/

int main(int argc, char *argv[]){
	INFO_PROCAR info;
	info.n_kpt = 8000;
	info.n_ions = 5;
	info.n_bands = 200;
	info.ispin = 1;
	info.soc = 1;
	info.space = 0;
	//print_info(info);
	//printf("%s\n", argv[1]);
	//char path[200] = "/home/gosteau/Documents/These/1ere année/STO/OUTPUT_dir_ncl/3.9/BS_SOC/";
	char path[200] = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/constrain/a375_c454/BS/";
	snprintf(path, 200,"%s", argv[1]);
	printf("%s\n", path);
	info = get_info(path);
	print_info(info);
	split_procar(path, info);
	printf("Done\n");
}

