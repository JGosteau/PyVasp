/*
 * split.c
 *
 *  Created on: 7 juin 2019
 *      Author: gosteau
 */


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include <locale.h>
#include"procar_info.h"
#define N 20000

void split_procar(char * path, INFO_PROCAR info){
	FILE * fp, * rfp;
	//setlocale(LC_ALL,"en_GB");
	//INFO_PROCAR info;
	//info = get_info(path);
	char * buff;
	int buffsn;
	int kpt;
	int write_kpt = 0;
	//char kx[10],ky[10],kz[10],weight[10];
	/*
	char *kx,*ky,*kz,*weight;
	kx = (char *) malloc(sizeof(char)*11);
	ky = (char *) malloc(sizeof(char)*11);
	kz = (char *) malloc(sizeof(char)*11);
	weight = (char *) malloc(sizeof(char)*11);
	*/
	double *kx,*ky,*kz,*weight;
	kx = malloc(sizeof(double));
	ky = malloc(sizeof(double));
	kz = malloc(sizeof(double));
	weight = malloc(sizeof(double));
	char file[N];
	char line[N];
	//char *line;
	//line = (char *) malloc(sizeof(char)*(info.n_char*(info.n_ions +1)));
	//line = (char *) malloc(sizeof(char)*N;
	char name[N];
	//char sp[3] = "";
	char *sp;
	sp = malloc(sizeof(char)*10);
	char * soc[3];
	soc[0] = malloc(sizeof(char)*10);

	buffsn = snprintf(soc[0], 10,"_sx");
	soc[1] = malloc(sizeof(char)*10);
	buffsn = snprintf(soc[1], 10,"_sy");
	soc[2] = malloc(sizeof(char)*10);
	buffsn = snprintf(soc[2], 10,"_sz");
	int count = 0, nband;
	//char Energy[14], occ[11];
	char *Energy, *occ;
	Energy = (char *) malloc(sizeof(char)*14);
	occ = (char *) malloc(sizeof(char)*11);

	buffsn = snprintf(file,N,"%sPROCAR", path);
	fp = fopen(file,"r");
	if (info.ispin == 1){
		buffsn = snprintf(name,N,"%sBANDS/band_info_tmp", path);
		rfp = fopen(name,"w");
		fclose(rfp);
	}
	else {
		buffsn = snprintf(name,N,"%sBANDS/band_info_tmp_up", path);
		rfp = fopen(name,"w");
		fclose(rfp);
		buffsn = snprintf(name,N,"%sBANDS/band_info_tmp_dn", path);
		rfp = fopen(name,"w");
		fclose(rfp);
	}

	buffsn = snprintf(name,N,"%sBANDS/kpts_tmp", path);
	rfp = fopen(name,"w");
	fclose(rfp);

	for(int i = 1; i < info.n_bands+1; i++){
		if(info.soc == 0){
			if(info.ispin == 1){
				buffsn = snprintf(name,N,"%sBANDS/BAND_%d", path,i);
				rfp = fopen(name,"w");
				fclose(rfp);
			}
			else{
				buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,i, "_up");
				rfp = fopen(name,"w");
				fclose(rfp);
				buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,i, "_dn");
				rfp = fopen(name,"w");
				fclose(rfp);
			}
		}
		else{
			buffsn = snprintf(name,N,"%sBANDS/BAND_%d", path,i);
			rfp = fopen(name,"w");
			fclose(rfp);
			buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,i, "_sx");
			rfp = fopen(name,"w");
			fclose(rfp);
			buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,i, "_sy");
			rfp = fopen(name,"w");
			fclose(rfp);
			buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,i, "_sz");
			rfp = fopen(name,"w");
			fclose(rfp);

		}
	}
	while( fgets(line,N,fp) ){
		if( line[0] == '#' && info.ispin == 2 ){
			if(count == 0){
				snprintf(sp, 10,"_up");
				write_kpt = 0;
				printf("%s\n", sp);
				count ++;
			}
			else {
				snprintf(sp, 10,"_dn");
				write_kpt = 1;
				printf("%s\n", sp);
			}
		}
		else {
			if( info.ispin == 1){
				snprintf(sp, 10,"");
			}
		}
		if(line[1] == 'k' && write_kpt == 0){
			//sscanf(line," k-point     %d :   %11s%11s%11s     weight =%11s", &kpt,kx, ky,kz,weight);
			sscanf(line," k-point     %d :   %lf%lf%lf     weight =%lf", &kpt,kx, ky,kz,weight);
			buffsn = snprintf(name,N,"%sBANDS/kpts_tmp", path);
			rfp = fopen(name,"a");
			//printf("%s %s %s %s\n", kx, ky,kz, weight);
			//buffsn = fprintf(rfp,"%s %s %s %s\n",kx, ky,kz,weight);
			buffsn = fprintf(rfp,"%1.8lf %1.8lf %1.8lf %1.8lf\n",*kx, *ky,*kz,*weight);
			fclose(rfp);
		}
		if(line[0] == 'b'){
			sscanf(line,"band %d # energy%14s # occ. %11s", &nband,Energy,occ);
			buffsn = snprintf(name,N,"%sBANDS/band_info_tmp%s", path, sp);


			rfp = fopen(name,"a");
			buffsn = fprintf(rfp,"%s %s\n",Energy,occ);
			fclose(rfp);
		}
		if(line[0] == 'i'){
			buffsn = fread(line, sizeof(char)*info.n_char, info.n_ions+1, fp);
			buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,nband, sp);
            rfp = fopen(name,"a");
            buffsn = fprintf(rfp,"%s", line);
            fclose(rfp);
            if(info.soc == 1){
            	for(int i =0; i<3;i++){
					if(info.space == 1){
						buff = fgets(line,N,fp);
					}
					buffsn = fread(line, sizeof(char)*info.n_char, info.n_ions+1, fp);
					buffsn = snprintf(name,N,"%sBANDS/BAND_%d%s", path,nband, soc[i]);
					rfp = fopen(name,"a");
					buffsn = fprintf(rfp,"%s", line);
					fclose(rfp);
            	}
            }
		}
	}
	fclose(fp);
}


