#include <math.h>
#include <stdio.h>
#include <string.h>
#include "../include/dataTypes.h"
#include "../include/geo.h"
#include "../include/fatal.h"
#include "../include/alloc.h"
#include "../include/msd.h"
#include "../include/geo.h"

int allocMSD(t_msd *msd,int nTimes)
{
	allocReals(&msd[0].msd,nTimes);
	msd[0].nTimes=nTimes;
	return 0;
}

int initMSD(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_msd *msd)
{
	int i;

	if(strcmp(mode,"atoms")==0)
	{
		allocInts(&msd[0].grp.atoms,grp.nAtoms);
		allocVecRs(&msd[0].ref,grp.nAtoms);
		allocVecRs(&msd[0].last,grp.nAtoms);
		msd[0].grp.nAtoms=grp.nAtoms;

		for(i=0;i<grp.nAtoms;i++)
		{
			msd[0].grp.atoms[i]=grp.atoms[i];
			vecRCpy(atoms[grp.atoms[i]].crd,&msd[0].ref[i]);
			vecRCpy(atoms[grp.atoms[i]].crd,&msd[0].last[i]);
		}
	} else if(strcmp(mode,"molCOM")==0)
	{
		allocInts(&msd[0].grp.mols,grp.nMols);
		allocVecRs(&msd[0].ref,grp.nMols);
                allocVecRs(&msd[0].last,grp.nMols);
		msd[0].grp.nMols=grp.nMols;

		for(i=0;i<grp.nMols;i++)
		{
			msd[0].grp.mols[i]=grp.mols[i];
			vecRCpy(mols[grp.mols[i]].COM,&msd[0].ref[i]);
			vecRCpy(mols[grp.mols[i]].COM,&msd[0].last[i]);
		}
	} else fatal("UNKNOWM MSD MODE: POSSIBLE SELECTIONS 'atoms' or 'molCOM'\n -> initMSD\n -> msd.c\n");

	return 0;
}

int initMSDtotal(t_msdTotal *msdTotal,int nTimes)
{
	int i;

	allocReals(&msdTotal[0].time,nTimes);
	allocReals(&msdTotal[0].msd,nTimes);
	for(i=0;i<nTimes;i++) msdTotal[0].msd[i]=0.0;
	msdTotal[0].nTimes=nTimes;
	msdTotal[0].cnt=0;

	return 0;
}

int reinitMSD(char *mode,t_atom *atoms,t_mol *mols,t_grp grp,t_msd *msd,t_msdTotal *msdTotal)
{
	int i;

	if(strcmp(mode,"atoms")==0)
	{
		if(msd[0].grp.nAtoms!=0)
		{
			for(i=0;i<msdTotal[0].nTimes;i++)
        		{
                		msdTotal[0].msd[i]+=msd[0].msd[i]/msd[0].grp.nAtoms;
        		}		
			msdTotal[0].cnt++;
		}

		msd[0].grp.atoms=(int*)save_realloc(msd[0].grp.atoms,grp.nAtoms*sizeof(int));
		msd[0].ref=(t_vecR*)save_realloc(msd[0].ref,grp.nAtoms*sizeof(t_vecR));
		msd[0].last=(t_vecR*)save_realloc(msd[0].last,grp.nAtoms*sizeof(t_vecR));
                msd[0].grp.nAtoms=grp.nAtoms;

                for(i=0;i<grp.nAtoms;i++)
                {
                        msd[0].grp.atoms[i]=grp.atoms[i];
                        vecRCpy(atoms[grp.atoms[i]].crd,&msd[0].ref[i]);
                        vecRCpy(atoms[grp.atoms[i]].crd,&msd[0].last[i]);
                }
	} else if(strcmp(mode,"molCOM")==0)
	{
		if(msd[0].grp.nMols!=0)
		{
			for(i=0;i<msdTotal[0].nTimes;i++)
                	{
                	        msdTotal[0].msd[i]+=msd[0].msd[i]/msd[0].grp.nMols;
                	}
			msdTotal[0].cnt++;
		}

		msd[0].grp.mols=(int*)save_realloc(msd[0].grp.mols,grp.nMols*sizeof(int));
                msd[0].ref=(t_vecR*)save_realloc(msd[0].ref,grp.nMols*sizeof(t_vecR));
                msd[0].last=(t_vecR*)save_realloc(msd[0].last,grp.nMols*sizeof(t_vecR));
                msd[0].grp.nMols=grp.nMols;

                for(i=0;i<grp.nMols;i++)
                {
                        msd[0].grp.mols[i]=grp.mols[i];
                        vecRCpy(mols[grp.mols[i]].COM,&msd[0].ref[i]);
                        vecRCpy(mols[grp.mols[i]].COM,&msd[0].last[i]);
                }
	} else fatal("UNKNOWM MSD MODE: POSSIBLE SELECTIONS 'atoms' or 'molCOM'\n -> reinitMSD\n -> msd.c\n");

	return 0;
}

int doMSD(char *mode,t_atom *atoms,t_mol *mols,t_box box,t_msd *msd,int step)
{
	int i;
	t_vecR tmp,tmp2;
	real single;

	msd[0].msd[step]=0.0;
	
	if(strcmp(mode,"atoms")==0)
        {
		for(i=0;i<msd[0].grp.nAtoms;i++)
		{
			/*vecRSub(msd[0].last[i],atoms[msd[0].grp.atoms[i]].crd,&tmp);
			transXYZToPBCbasis(box,tmp,&tmp2);
			if(tmp2.x>0.5 || tmp2.x<-0.5) tmp2.x-=floor(tmp2.x+0.5);
			if(tmp2.y>0.5 || tmp2.y<-0.5) tmp2.y-=floor(tmp2.y+0.5);
			if(tmp2.z>0.5 || tmp2.z<-0.5) tmp2.z-=floor(tmp2.z+0.5);
			transPBCToXYZbasis(box,tmp2,&tmp);*/
			linkPBC(msd[0].last[i],atoms[msd[0].grp.atoms[i]].crd,box,&tmp);
			vecRAdd(msd[0].last[i],tmp,&tmp2);
			vecRCpy(tmp2,&msd[0].last[i]);

			vecRSub(msd[0].last[i],msd[0].ref[i],&tmp);
			vecRNormSq(tmp,&single);

			msd[0].msd[step]+=single;
		}
	} else if(strcmp(mode,"molCOM")==0)
	{
		for(i=0;i<msd[0].grp.nMols;i++)
                { 
                        /*vecRSub(msd[0].last[i],mols[msd[0].grp.mols[i]].COM,&tmp);
                        transXYZToPBCbasis(box,tmp,&tmp2);
			if(tmp2.x>0.5 || tmp2.x<-0.5) tmp2.x-=floor(tmp2.x+0.5);
                        if(tmp2.y>0.5 || tmp2.y<-0.5) tmp2.y-=floor(tmp2.y+0.5);
                        if(tmp2.z>0.5 || tmp2.z<-0.5) tmp2.z-=floor(tmp2.z+0.5);
                        transPBCToXYZbasis(box,tmp2,&tmp);*/
			linkPBC(msd[0].last[i],mols[msd[0].grp.mols[i]].COM,box,&tmp);
                        vecRAdd(msd[0].last[i],tmp,&tmp2);
                        vecRCpy(tmp2,&msd[0].last[i]);

                        vecRSub(msd[0].last[i],msd[0].ref[i],&tmp);
                        vecRNormSq(tmp,&single);

                        msd[0].msd[step]+=single;
                }
	} else fatal("UNKNOWM MSD MODE: POSSIBLE SELECTIONS 'atoms' or 'molCOM'\n -> doMSD\n -> msd.c\n");

	return 0;
}

