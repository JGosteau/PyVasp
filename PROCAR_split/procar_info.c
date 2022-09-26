/*
 * procar_info.c
 *
 *  Created on: 6 juin 2019
 *      Author: gosteau
 */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
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


int cfileexists(const char * filename){
    /* try to open file to read */
    FILE *file;
    if (file = fopen(filename, "r")){
        fclose(file);
        return 1;
    }
    return 0;
}

void print_info(INFO_PROCAR info){
	printf("n_kpt : %d\nn_ions : %d\nn_bands : %d\nn_char : %d\nispin : %d\nLSORBIT : %d\nsoc_space : %d\n",
			info.n_kpt, info.n_ions, info.n_bands, info.n_char, info.ispin, info.soc, info.space);
}

int get_ISPIN(char * path, INFO_PROCAR info){
	FILE *fp;
	char file[N];
	char line[N];
	int count = 0;
	int ispin = -1;
	char ISPIN[5];
	ispin = 1;
	if (info.soc != 1) {
		snprintf(file,N,"%sOUTCAR", path);
		if (cfileexists(file) == 1){
			fp = fopen(file,"r");
			while(count != 1 && fgets(line,N, fp)){
				sscanf(line, " %s = %d", ISPIN,&ispin);
				if(strcmp("ISPIN",ISPIN) == 0){
					//printf("%d , %d | %s" ,strcmp("ISPIN",ISPIN),ispin,line);
					count = 1;
				}

			}
		}
		else{
			snprintf(file,N,"%sPROCAR", path);
			fp = fopen(file,"r");
			while(fgets(line,N,fp)){
				count += 1;
			}
			fclose(fp);
			printf("%d || %d\n", count, info.n_bands*(info.n_ions+1)*info.n_kpt*2);
			if(count > info.n_bands*(info.n_ions+1)*info.n_kpt*2){
				ispin = 2;
			}
		}
		/*
		snprintf(file,N,"%sOUTCAR", path);
		fp = fopen(file,"r");
		while(count != 1 && fgets(line,N, fp)){
			sscanf(line, " %s = %d", ISPIN,&ispin);
			if(strcmp("ISPIN",ISPIN) == 0){
				//printf("%d , %d | %s" ,strcmp("ISPIN",ISPIN),ispin,line);
				count = 1;
			}
		}
		*/
	}
	return ispin;
}

INFO_PROCAR get_info(char * path){
	FILE *fp;
	INFO_PROCAR info;
	char file[N];
	char line[N];
	char * buff;
	int buff_int;
	int count = 0;
	snprintf(file,N,"%sPROCAR", path);
	fp = fopen(file,"r");
	while(count < 2){
		buff = fgets(line,N, fp);
		//printf("%s", line);
		if(line[0] == '#'){
			sscanf(line, "# of k-points: %d # of bands:  %d # of ions: %d", &info.n_kpt, &info.n_bands, &info.n_ions);
		}
		if(line[0] == 'b'){
			count ++;
		}
		if(line[0] == 'i'){
			/*
			printf("%s\n",line);
			while(i < strlen(line)){
				while(*(line+i) == ' ') {
					//printf("%c\n", *(line+i));
					i++;
				}
				sscanf(line+i, "%s ", tmp);
				printf(" %s ",tmp);
				i += strlen(tmp);
			}
			printf("\n");
			*/
			snprintf(info.orbs, N, "%s",line);
			buff = fgets(line,N,fp);
			//printf("gets : \n%s", line);
			info.n_char = strlen(line);
			buff_int = fread(line,info.n_char,info.n_ions,fp);
			//printf("read : \n%s", line);
			buff = fgets(line,N,fp);
			//printf("line : %s END %d\n", line, strlen(line));
			if(strlen(line) == 1){
				info.space = 1;
				buff = fgets(line,N,fp);
			}
			else {
				info.space = 0;
			}
			if(strlen(line) != 2){
				info.soc = 1;
				info.ispin = 1;
			}
			else {
				info.soc = 0;
				info.ispin = get_ISPIN(path, info);
			}
			count ++;

		}


	}
	fclose(fp);
	return info;
}
