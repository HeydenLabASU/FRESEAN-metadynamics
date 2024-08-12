#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/dataTypes.h"
#include "../include/alloc.h"
#include "../include/io.h"
#include "../include/fatal.h"

int printTitle() {
        printf("THE PURPOSE OF THIS PROGRAM IS TO WRITE\n");
	printf("A SELECTED EIGVECTOR IN ASCII\n");
        printf("Version 1.0: Sept. 22, 2021\n");
        printf("Author:\n");
        printf(" Dr. Matthias Heyden\n");
        printf(" School of Molecular Sciences\n");
        printf(" Arizona State University\n");
        printf(" Tempe, AZ, USA\n");
        printf(" e-mail: mheyden1@asu.edu\n");
        printf("\n");
        return 0;
}

int getLineFromCOM(FILE *in,char *buffer,int max)
{
        if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        while(strncmp(buffer,"#",1)==0)
        {
                if(fgets(buffer,max,in)==NULL) fatal("input file incomplete\n");
        }
        return 0;
}

int printKeys() {
	printf("fnEigVec\nfreqSel\nmodeStart\nmodeEnd\nfnOut\n");
        return 0;
}

int getInput(char *fnCOM,char *fnEigVec,int *freqSel, float *timestep, int *modeStart,int *modeEnd, int *extractMode, char *fnOut) {
        FILE *io;
        char buffer[300];
        int format;
        int i=1;
        float tmp;

        saveOpenRead(&io,fnCOM);

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnEigVec);
        printf("%2d -> read %20s : %s\n",i,"fnEigVec",fnEigVec);i++;

	getLineFromCOM(io,buffer,300);
	sscanf(buffer,"%d",extractMode);
        printf("%2d -> read %20s : %d\n",i,"extractMode",extractMode[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",freqSel);
        printf("%2d -> read %20s : %d\n",i,"freqSel",freqSel[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%f",timestep);
        printf("%2d -> read %20s : %f\n",i,"timestep",timestep[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",modeStart);
        printf("%2d -> read %20s : %d\n",i,"modeStart",modeStart[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%d",modeEnd);
        printf("%2d -> read %20s : %d\n",i,"modeEnd",modeEnd[0]);i++;

	getLineFromCOM(io,buffer,300);
        sscanf(buffer,"%s",fnOut);
        printf("%2d -> read %20s : %s\n",i,"fnOut",fnOut);i++;

        fclose(io);
        return 0;
}

typedef struct {
	double *m;
	int nu,nv;
} t_mat;

int getMat(t_mat mat,int u,int v,double *muv) {
	int i;
	i=u*mat.nv+v;
	muv[0]=mat.m[i];
	return 0;
}

int setMat(t_mat mat,int u,int v,double muv) {
	int i;
	i=u*mat.nv+v;
        mat.m[i]=muv;
        return 0;
}

int main(int argc, char *argv[])
{
	FILE *io;
	char fnCOM[100],fnOut[100];
	char fnOut2[100];
	int i,j;
/*MODULAR*/
	char fnEigVec[100];
	int nCorr,n1,n2;
	double **eigVec;
	int nProj;
	double **proj;
	double tmp1,tmp2,tmp3;
	float floatRead;
	int freqSel;
	int modeStart,modeEnd;
	int bc1,bc2;
	float df;
	float timestep;
	int extractMode;
/*MODULAR*/

	printTitle();

/*PARTIALLY MODULAR*/
	/*PARSE COMMAND LINE INPUT HERE*/
	if(argc!=2) {
                printKeys();
                fatal("no input file specified\n");
        }
        strcpy(fnCOM,argv[1]);
        getInput(fnCOM,fnEigVec, &freqSel, &timestep, &modeStart, &modeEnd, &extractMode, fnOut);
/*PARTIALLY MODULAR*/

/*MODULAR*/
	nProj=(modeEnd-modeStart)+1;
	io=fopen(fnEigVec,"rb");
	
	// Read in the first two integers from the binary file.
	// What is bc1? I guess bc1=bc2 is to check a condition 	
	fread(&bc1,sizeof(int),1,io);

	// nCorr determines the max value of tau to use in the velocity cross correlation matrix
	// nCorr is unitless (number of frames)
	// To determine the time tau, need to multiply by the time step
	fread(&nCorr,sizeof(int),1,io);

	// Determine the frequency resolution
	df = 33/((2*nCorr-1)*timestep);

	// If extractMode is set to 0, then we will take in freqSel as an integer corresponding to the desired frequency
	// Since the frequency resolution is most likely to be a floating point value
	// We will compare freqSel to the frequency resolution to determine the closest floating point
	// the closest floating point frequency
	//
	// Then we will convert to an index and proceeed as normal
	if (extractMode == 0){
		printf("Frequency resolution is %.5f\n", df);
		freqSel = round(freqSel/df)+1;
		printf("Converting for frequency %.5f cm^-1 or index %d\n", df*(freqSel-1), freqSel);
	}

	// We can only look at nCorr number of frequencies
	if(nCorr<freqSel) {
		printf("not enough matrices in %s: freqSel %d nCorr %d\n",fnEigVec,freqSel,nCorr);
		exit(1);
	}

	//Read in the next three integers in the header
	//n1 = n2 = 3*DOF
	//bc2 = bc1
	fread(&n1,sizeof(int),1,io);
	fread(&n2,sizeof(int),1,io);
	fread(&bc2,sizeof(int),1,io);
	if(bc1!=bc2) {
		printf("format error in %s\n",fnEigVec);
		exit(1);
	}
	if(n1!=n2) {
		printf("incompatible dimensions in %s: nCorr %d n1 %d n2 %d\n",fnEigVec,nCorr,n1,n2);
		exit(1);
	}


	eigVec=(double**)save_malloc(nProj*sizeof(double*));
		for(i=0;i<nProj;i++) {
			eigVec[i]=(double*)save_malloc(n1*sizeof(double));
		}
	for(i=1;i<freqSel;i++) {
		fread(&bc1,sizeof(int),1,io);
		for(j=0;j<n1*n2;j++) {
			fread(&tmp1,sizeof(double),1,io);
		}
		fread(&bc2,sizeof(int),1,io);
		if(bc1!=bc2) {
			printf("format error in %s\n",fnEigVec);
			exit(1);
		}
	}
	fread(&bc1,sizeof(int),1,io);
	for(i=1;i<modeStart;i++) {
		for(j=0;j<n1;j++) {
			fread(&tmp1,sizeof(double),1,io);
		}
	}
	for(i=0;i<nProj;i++) {
		for(j=0;j<n1;j++) {
			fread(&tmp1,sizeof(double),1,io);
				eigVec[i][j]=tmp1;
			}
	}
	fclose(io);
	sprintf(fnOut2,"evec_freq%d_mode%d-%d_%s",freqSel,modeStart,modeEnd,fnOut);
	io=fopen(fnOut2,"w");
	for(i=0;i<nProj;i++) {
		fprintf(io,"%d\nFrequency %d Eigenvector %d\n",n1/3,freqSel,modeStart+i);
		for(j=0;j<n1;j++) {
			if(j%3==0) fprintf(io,"X  ");
			fprintf(io," %15e",eigVec[i][j]);
			if((j+1)%3==0) fprintf(io,"\n");
		}
	}


	fclose(io);
	return 0;
}

