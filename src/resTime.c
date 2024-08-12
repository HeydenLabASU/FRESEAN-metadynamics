#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/dataTypes.h"
#include "../include/geo.h"
#include "../include/fatal.h"
#include "../include/alloc.h"
#include "../include/qsort.h"
#include "../include/select.h"
#include "../include/grps.h"
#include "../include/resTime.h"

int allocResTime(t_resTime *resTime,int nTimes)
{
        allocReals(&resTime[0].resTime,nTimes);
        resTime[0].nTimes=nTimes;
        return 0;
}

int initResTime(t_grp grp,t_resTime *resTime)
{
        int i;

	allocInts(&resTime[0].grp.atoms,grp.nAtoms);
        allocInts(&resTime[0].grp.mols,grp.nMols);
	resTime[0].grp.nAtoms=grp.nAtoms;
        resTime[0].grp.nMols=grp.nMols;

	for(i=0;i<grp.nAtoms;i++)
        {
                resTime[0].grp.atoms[i]=grp.atoms[i];
        }
        for(i=0;i<grp.nMols;i++)
        {
                resTime[0].grp.mols[i]=grp.mols[i];
        }

        return 0;
}

int initResTimeTotal(t_resTimeTotal *resTimeTotal,int nTimes)
{
        int i;

        allocReals(&resTimeTotal[0].time,nTimes);
        allocReals(&resTimeTotal[0].resTime,nTimes);
        for(i=0;i<nTimes;i++) resTimeTotal[0].resTime[i]=0.0;
        resTimeTotal[0].nTimes=nTimes;
        resTimeTotal[0].cnt=0;

        return 0;
}

int reinitResTime(t_grp grp,t_resTime *resTime,t_resTimeTotal *resTimeTotal)
{
        int i;

        for(i=0;i<resTimeTotal[0].nTimes;i++)
        {
                resTimeTotal[0].resTime[i]+=resTime[0].resTime[i];
        }

        resTimeTotal[0].cnt++;

	resTime[0].grp.atoms=(int*)save_realloc(resTime[0].grp.atoms,grp.nAtoms*sizeof(int));
        resTime[0].grp.mols=(int*)save_realloc(resTime[0].grp.mols,grp.nMols*sizeof(int));
	resTime[0].grp.nAtoms=grp.nAtoms;
        resTime[0].grp.nMols=grp.nMols;

	for(i=0;i<grp.nAtoms;i++)
        {
                resTime[0].grp.atoms[i]=grp.atoms[i];
        }
        for(i=0;i<grp.nMols;i++)
        {
                resTime[0].grp.mols[i]=grp.mols[i];
        }

        return 0;
}

int doResTime(t_atom *atoms,t_box box,int nAtoms,t_mol *mols,int nMols,t_statGrpDef *statGrpDef,int nStat,t_dynGrpDef dynGrpDef,t_grp *grps,t_resTime *resTime,int step)
{
	int i,idx;
	char err[300];
	t_grp tmp;

	tmp.nAtoms=resTime[0].grp.nAtoms;
	tmp.nMols=resTime[0].grp.nMols;
	tmp.atoms=(int*)save_malloc(tmp.nAtoms*sizeof(int));
	tmp.mols=(int*)save_malloc(tmp.nMols*sizeof(int));

	if(dynGrpDef.preGrpIdx>=nStat)
		fatal("RESIDENCE TIMES ONLY SUPPORTED FOR DYNAMIC GROUPS WITH A STATIC PRESELECTION GROUP!\n -> doResTime\n -> resTime.c\n");

	for(i=0;i<nStat;i++)
	{
		getGrpCOM(atoms,&grps[i]);
	}

	if(dynGrpDef.doRef==1 && dynGrpDef.doDist==1)
	{
		if(strcmp(dynGrpDef.refMode,"COM")==0 && strcmp(dynGrpDef.selMode,"atoms")==0)
		{
			dynAtomsByDistToCrd(grps[dynGrpDef.refGrpIdx].COM,dynGrpDef,atoms,box,resTime[0].grp,&tmp);
		} else if(strcmp(dynGrpDef.refMode,"COM")==0 && strcmp(dynGrpDef.selMode,"mol_by_COM")==0)
		{ 
			dynMolsCOMByDistToCrd(grps[dynGrpDef.refGrpIdx].COM,dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
		} else if(strcmp(dynGrpDef.refMode,"COM")==0 && strcmp(dynGrpDef.selMode,"mol_by_atom")==0)
		{
			dynMolsAtomByDistToCrd(grps[dynGrpDef.refGrpIdx].COM,dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
		} else if(strcmp(dynGrpDef.refMode,"closest_atom")==0 && strcmp(dynGrpDef.selMode,"atoms")==0)
                {
			dynAtomsByDistToGrp(grps[dynGrpDef.refGrpIdx],dynGrpDef,atoms,box,resTime[0].grp,&tmp);
                } else if(strcmp(dynGrpDef.refMode,"closest_atom")==0 && strcmp(dynGrpDef.selMode,"mol_by_COM")==0)
                {
			dynMolsCOMByDistToGrp(grps[dynGrpDef.refGrpIdx],dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
                } else if(strcmp(dynGrpDef.refMode,"closest_atom")==0 && strcmp(dynGrpDef.selMode,"mol_by_atom")==0)
                {
			dynMolsAtomByDistToGrp(grps[dynGrpDef.refGrpIdx],dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
                }
	} else if(dynGrpDef.doCrd==1 && dynGrpDef.doDist==1)
	{
		if(strcmp(dynGrpDef.selMode,"atoms")==0)
		{
			dynAtomsByDistToCrd(dynGrpDef.crd,dynGrpDef,atoms,box,resTime[0].grp,&tmp);
		} else if(strcmp(dynGrpDef.selMode,"mol_by_COM")==0)
		{
			dynMolsCOMByDistToCrd(dynGrpDef.crd,dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
		} else if(strcmp(dynGrpDef.selMode,"mol_by_atom")==0)
		{
			dynMolsAtomByDistToCrd(dynGrpDef.crd,dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
		}
	} else if(dynGrpDef.doX==1 || dynGrpDef.doY==1 || dynGrpDef.doZ==1)
	{
		if(strcmp(dynGrpDef.selMode,"atoms")==0)
		{
			dynAtomsByXYZ(dynGrpDef,atoms,box,resTime[0].grp,&tmp);
		} else if(strcmp(dynGrpDef.selMode,"mol_by_COM")==0)
		{
			dynMolsCOMByXYZ(dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
		} else if(strcmp(dynGrpDef.selMode,"mol_by_atom")==0)
		{
			dynMolsAtomByXYZ(dynGrpDef,atoms,mols,box,resTime[0].grp,statGrpDef[dynGrpDef.preGrpIdx],&tmp);
		}
	} else {
		sprintf(err,"NO VALID SELECTION PARAMETERS FOR DYNAMIC GROUP %s!\n -> \n -> resTime.c\n",dynGrpDef.name);
		fatal(err);
	}

	if(strcmp(dynGrpDef.selMode,"atoms")==0)
	{
		if(resTime[0].grp.nAtoms!=0) resTime[0].resTime[step]=((real)tmp.nAtoms)/((real)resTime[0].grp.nAtoms);
		else resTime[0].resTime[step]=0.0;
	} else {
		if(resTime[0].grp.nMols!=0) resTime[0].resTime[step]=((real)tmp.nMols)/((real)resTime[0].grp.nMols);
		else resTime[0].resTime[step]=0.0;
	}

	free(tmp.atoms);
	free(tmp.mols);

	return 0;
}

