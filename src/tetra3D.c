#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/*#include <malloc.h>*/
#include "../include/dataTypes.h"
#include "../include/geo.h"
#include "../include/fatal.h"
#include "../include/alloc.h"
#include "../include/tetra3D.h"

int allocTetra(t_atom *atoms,int nAtoms,t_mol *mols,t_tetra *tetra)
{
	int i,j,curmol,donCatch,donor;
	char err[300];

	allocInts(&tetra[0].coordList,nAtoms);

	for(i=0;i<nAtoms;i++)
	{
		atoms[i].HBinfo[0]=0;
		atoms[i].HBinfo[1]=0;
	}
	for(i=0;i<nAtoms;i++)
	{
		if(strncmp(atoms[i].atomName,"H",1)==0 && atoms[i].charge>=0.2)
		{
			atoms[i].HBinfo[0]=1;
			curmol=atoms[i].resNr;
			donCatch=0;
			j=0;
			while(donCatch==0 && j<mols[curmol].nBonds)
			{
				if(mols[curmol].bonds[j][0]==i/*&& atoms[mols[curmol].bonds[j][1]].charge<0.0*/)
				{
					donCatch=1;
					donor=mols[curmol].bonds[j][1];
				} else if(mols[curmol].bonds[j][1]==i/*&& atoms[mols[curmol].bonds[j][0]].charge<0.0*/)
                                {
                                        donCatch=1;
                                        donor=mols[curmol].bonds[j][0];
                                }
				j++;
			}
			if(donCatch!=1) 
			{
				sprintf(err,"DID NOT FIND DONOR ATOM FOR POLARIZED HYDROGEN %d WITH CHARGE %f!\n -> allocTetra\n -> tetra.c\n",i,atoms[i].charge);
				fatal(err);
			} else {
				atoms[i].HBinfo[1]=donor;
				/*this is detrimental in an actual HB analysis*/
				/*because atoms can be donors and acceptors simultaneously,*/
				/* but here it is ok*/
				atoms[donor].HBinfo[0]=2;
			}
		}
		if(atoms[i].charge<=-0.2)
		{
			atoms[i].HBinfo[0]=3;  
		}
	}

	return 0;
}

int gridPos(t_vecR crd,t_vecR ori,real dX,real dY,real dZ,int *pos,t_box box,int *ng) {
        t_vecR link;
        int check;
        real tmp1;
        int tmp2;

        /*grid ori is the origin of the grid (one of its corners, not the center)*/
        //linkPBC(ori,crd,box,&link);
        vecRSub(ori,crd,&link);

        if(link.x<0.0) link.x+=1e-5;
        if(link.y<0.0) link.y+=1e-5;
        if(link.z<0.0) link.z+=1e-5;
        if(link.x<0.0 || link.y<0.0 || link.z<0.0) {
                check=0; /*outside of grid volume*/
        } else {
                tmp1=link.x/dX;
                tmp2=(int)floor(tmp1);
                if(tmp2>=ng[0]) { tmp1-=1e-5; tmp2=(int)floor(tmp1); }
                pos[0]=tmp2;
                if(pos[0]>=ng[0]) {
                        check=0; /*outside of grid volume*/
                } else {
                        tmp1=link.y/dY;
                        tmp2=(int)floor(tmp1);
                        if(tmp2>=ng[1]) { tmp1-=1e-5; tmp2=(int)floor(tmp1); }
                        pos[1]=tmp2;
                        if(pos[1]>=ng[1]) {
                                check=0; /*outside of grid volume*/
                        } else {
                                tmp1=link.z/dZ;
                                tmp2=(int)floor(tmp1);
                                if(tmp2>=ng[2]) { tmp1-=1e-5; tmp2=(int)floor(tmp1); }
                                pos[2]=tmp2;
                                if(pos[2]>=ng[2]) {
                                        check=0; /*outside of grid volume*/
                                } else {
                                        check=1; /*INSIDE THE GRID VOLUME*/
                                }
                        }
                }
        }
        return check;
}

int doTetra(t_atom *atoms,int nAtoms,t_box box,t_grp centerGrp,t_grp coordGrp,t_tetra *tetra,t_vecR gridOri,int *ng,real dgrid,double *vox,int *vox2,int *vox3)
{
	FILE *err;
	int i,j,k;
	int m,n;
	real dSq;
	int closeList[1000];
	real closeDistSq[1000];
	int nCloseList,nCloseListCnt;
	int four[4];
	real fourDSq[4];
	real minDSq;
	int minIdx;
	real order,angle,tmp;
	int pos[3];
	int pos2,check;

	j=0;
	for(i=0;i<coordGrp.nAtoms;i++)
        {
		m=coordGrp.atoms[i];
		if(atoms[m].HBinfo[0]>1) {
			tetra[0].coordList[j]=m;
			j++;
		}
	}
	tetra[0].nCoordList=j;

	nCloseListCnt=0;
	for(i=0;i<centerGrp.nAtoms;i++) {
		m=centerGrp.atoms[i];
		check=gridPos(atoms[m].crd,gridOri,dgrid,dgrid,dgrid,pos,box,ng);
		if(check==1) {
			pos2=pos[2]+ng[1]*(pos[1]+ng[0]*pos[0]);
			vox2[pos2]++;
			vox3[pos2]++;
			k=0;
			for(j=0;j<tetra[0].nCoordList;j++) {
				n=tetra[0].coordList[j];
				distSqPBC(atoms[m].crd,atoms[n].crd,box,&dSq);
				if(dSq<=16.0 && m!=n) {
					closeList[k]=n;
					closeDistSq[k]=dSq;
					k++;
				}
			}
			nCloseList=k;
			nCloseListCnt+=k;
			/*if(nCloseList<4) {
				err=fopen("tetra_environment.gro","w");
				fprintf(err,"check environment of atom %d\n",m);
				fprintf(err,"%d\n",nAtoms);
				for(j=0;j<nAtoms;j++) {
					fprintf(err,"%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n",
						atoms[j].resNr+1,
						atoms[j].resName,
						atoms[j].atomName,
						(atoms[j].atomNr+1)%100000,
						0.1*atoms[j].crd.x,
						0.1*atoms[j].crd.y,
						0.1*atoms[j].crd.z
					);
				}
				fprintf(err,"%10.5f%10.5f%10.5f\n",0.0,0.0,0.0);
				fclose(err);
				fatal("encountered center atom with less than 4 neighbors within cutoff\n");
			}*/
	
			/*we skip the test above and just ignore atoms without 4 coordinating atoms withing the cutoff*/
	                /*for that case, we set the order parameter to -1 and remove it from the average (see 'else' statement below)*/
			if(nCloseList>=4) {	
				for(k=0;k<4;k++) {
					minDSq=9999.9;
					minIdx=0;
					for(j=0;j<nCloseList;j++) {
						if(closeDistSq[j]<minDSq) {
							minDSq=closeDistSq[j];
							minIdx=j;
						}
					}
					four[k]=closeList[minIdx];
					fourDSq[k]=closeDistSq[minIdx];
					closeDistSq[minIdx]=9999.9;
				}
				order=0.0;
				for(j=0;j<3;j++) { 
					for(k=j+1;k<4;k++) {
						getAnglePBC(atoms[four[j]].crd,atoms[four[k]].crd,atoms[m].crd,box,&angle);
						tmp=cos(angle)+1.0/3.0;
						order+=tmp*tmp;
					}
				}
				order=1-3.0/8.0*order;
				if(order<-1.0 || order>1.0) {
					for(j=0;j<4;j++) {
						printf("atom: %s\n",atoms[four[j]].atomName);
					}
					for(j=0;j<4;j++) {
						distPBC(atoms[m].crd,atoms[four[j]].crd,box,&angle);
						printf("dist: %f\n",angle);
					}
					for(j=0;j<3;j++) { for(k=j+1;k<4;k++) {
						getAnglePBC(atoms[four[j]].crd,atoms[four[k]].crd,atoms[m].crd,box,&angle);
						printf("angle: %f\n",angle/(2.0*3.14159)*360.0);
					}}
					printf("order: %f\n",order);
					fatal("order parameter out of range\n");
				}
				vox[pos2]+=order;
			} else {
				vox3[pos2]--;
			}
		}
	}
	/*printf("analyzed: %d centers with %d potential coordinators (average %5.2f within cutoff)\n",
		centerGrp.nAtoms,
		tetra[0].nCoordList,
		((float)nCloseListCnt)/((float)centerGrp.nAtoms));fflush(stdout);*/

	return 0;
}

