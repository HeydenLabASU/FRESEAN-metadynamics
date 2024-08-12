#include <math.h>
#include <stdio.h>
#include "../include/dataTypes.h"
#include "../include/geo.h"
#include "../include/fatal.h"
#include "../include/alloc.h"
#include "../include/msf.h"

int allocMSF(t_msf *msf,int nTimes)
{
	msf[0].nTimes=nTimes;
	return 0;
}

int initMSF1(t_atom *atoms,t_grp grp,t_msf *msf)
{
	int i;

	allocInts(&msf[0].grp.atoms,grp.nAtoms);
	allocVecRs(&msf[0].ref,grp.nAtoms);
	allocVecRs(&msf[0].last,grp.nAtoms);
	msf[0].grp.nAtoms=grp.nAtoms;

	for(i=0;i<grp.nAtoms;i++)
	{
		msf[0].grp.atoms[i]=grp.atoms[i];
		msf[0].ref[i].x=0.0; msf[0].ref[i].y=0.0; msf[0].ref[i].z=0.0;
		vecRCpy(atoms[grp.atoms[i]].crd,&msf[0].last[i]);
	}

	msf[0].msf=0.0;

	return 0;
}

int initMSF2(FILE *ref,t_atom *atoms,t_grp grp,t_msf *msf)
{
	int i;
	char dum;
	float x,y,z;
	char buffer[200];
	char err[200];

	fgets(buffer,200,ref);
	sscanf(buffer,"%d",&i);
	fgets(buffer,200,ref);
	if(i!=grp.nAtoms) {
		sprintf(err,"MSF error: number of atoms read from reference file (%d) does not match current group (%d)!!!\n -> initMSF2\n -> msf.c\n",i,grp.nAtoms);
		fatal(err);
	}
	msf[0].grp.atoms=(int*)save_realloc(msf[0].grp.atoms,grp.nAtoms*sizeof(int));
	msf[0].ref=(t_vecR*)save_realloc(msf[0].ref,grp.nAtoms*sizeof(t_vecR));
	msf[0].last=(t_vecR*)save_realloc(msf[0].last,grp.nAtoms*sizeof(t_vecR));
	msf[0].grp.nAtoms=grp.nAtoms;

	for(i=0;i<grp.nAtoms;i++)
	{
		fgets(buffer,200,ref);
		sscanf(buffer,"%c %f %f %f",&dum,&x,&y,&z);
		msf[0].ref[i].x=x; msf[0].ref[i].y=y; msf[0].ref[i].z=z;
		msf[0].grp.atoms[i]=grp.atoms[i];
		vecRCpy(atoms[grp.atoms[i]].crd,&msf[0].last[i]);
	}
	
	msf[0].msf=0.0;

	return 0;
}

int initMSFtotal(t_msfTotal *msfTotal)
{
	int i;

	msfTotal[0].msf=0.0;
	msfTotal[0].cnt=0;

	return 0;
}

int reinitMSF1(FILE *ref,t_atom *atoms,t_grp grp,t_msf *msf)
{
	int i;

	fprintf(ref,"%d\n\n",msf[0].grp.nAtoms);
	for(i=0;i<msf[0].grp.nAtoms;i++)
	{
		fprintf(ref,"%c%10.5f%10.5f%10.5f\n",atoms[msf[0].grp.atoms[i]].atomName[0],msf[0].ref[i].x/msf[0].nTimes,msf[0].ref[i].y/msf[0].nTimes,msf[0].ref[i].z/msf[0].nTimes);
	}

	msf[0].grp.atoms=(int*)save_realloc(msf[0].grp.atoms,grp.nAtoms*sizeof(int));
	msf[0].ref=(t_vecR*)save_realloc(msf[0].ref,grp.nAtoms*sizeof(t_vecR));
	msf[0].last=(t_vecR*)save_realloc(msf[0].last,grp.nAtoms*sizeof(t_vecR));
	msf[0].grp.nAtoms=grp.nAtoms;

	for(i=0;i<grp.nAtoms;i++)
	{
		msf[0].grp.atoms[i]=grp.atoms[i];
		msf[0].ref[i].x=0.0; msf[0].ref[i].y=0.0; msf[0].ref[i].z=0.0;
		vecRCpy(atoms[grp.atoms[i]].crd,&msf[0].last[i]);
	}
	return 0;
}

int reinitMSF2(FILE *ref,t_atom *atoms,t_grp grp,t_msf *msf,t_msfTotal *msfTotal)
{
	int i;
	char dum;
	float x,y,z;
	char buffer[200];
	char err[300];

	if(msf[0].grp.nAtoms!=0)
	{
               	msfTotal[0].msf+=msf[0].msf/msf[0].nTimes;
		msfTotal[0].cnt++;
	}

        if(fgets(buffer,200,ref)!=NULL)
	{
		sscanf(buffer,"%d",&i);
       		if(i!=grp.nAtoms) {
			sprintf(err,"MSF error: number of atoms read from reference file (%d) does not match current group (%d)!!!\n -> reinitMSF2\n -> msf.c\n",i,grp.nAtoms);
			fatal(err);
		}
		fgets(buffer,200,ref);

       		msf[0].grp.atoms=(int*)save_realloc(msf[0].grp.atoms,grp.nAtoms*sizeof(int));
        	msf[0].ref=(t_vecR*)save_realloc(msf[0].ref,grp.nAtoms*sizeof(t_vecR));
        	msf[0].last=(t_vecR*)save_realloc(msf[0].last,grp.nAtoms*sizeof(t_vecR));
        	msf[0].grp.nAtoms=grp.nAtoms;

        	for(i=0;i<grp.nAtoms;i++)
        	{
			fgets(buffer,200,ref);
			sscanf(buffer,"%c %f %f %f",&dum,&x,&y,&z);
			msf[0].ref[i].x=x; msf[0].ref[i].y=y; msf[0].ref[i].z=z;
        	        msf[0].grp.atoms[i]=grp.atoms[i];
        	        vecRCpy(atoms[grp.atoms[i]].crd,&msf[0].last[i]);
        	}
	}

        msf[0].msf=0.0;

	return 0;
}

int doMSF1(t_atom *atoms,t_box box,t_msf *msf)
{
        int i;
        t_vecR tmp,tmp2;

        for(i=0;i<msf[0].grp.nAtoms;i++)
        {
                linkPBC(msf[0].last[i],atoms[msf[0].grp.atoms[i]].crd,box,&tmp);
                vecRAdd(msf[0].last[i],tmp,&tmp2);
                vecRCpy(tmp2,&msf[0].last[i]);

                vecRAdd(msf[0].last[i],msf[0].ref[i],&tmp);
                vecRCpy(tmp,&msf[0].ref[i]);
        }

        return 0;
}

int doMSF2(t_atom *atoms,t_box box,t_msf *msf)
{
	int i;
	t_vecR tmp,tmp2;
	real single;

	for(i=0;i<msf[0].grp.nAtoms;i++)
	{
		linkPBC(msf[0].last[i],atoms[msf[0].grp.atoms[i]].crd,box,&tmp);
		vecRAdd(msf[0].last[i],tmp,&tmp2);
		vecRCpy(tmp2,&msf[0].last[i]);

		vecRSub(msf[0].last[i],msf[0].ref[i],&tmp);
		vecRNormSq(tmp,&single);

		msf[0].msf+=single/msf[0].grp.nAtoms;
	}

	return 0;
}

