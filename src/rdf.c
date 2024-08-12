#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../include/dataTypes.h"
#include "../include/geo.h"
#include "../include/fatal.h"
#include "../include/alloc.h"
#include "../include/rdf.h"

int initRDF(int grp1Idx,int grp2Idx,real binSize,int nBins,t_rdf *rdf)
{
        int i;

        rdf[0].grp1Idx=grp1Idx;
        rdf[0].grp2Idx=grp2Idx;
        rdf[0].binSize=binSize;
        rdf[0].nBins=nBins;
        rdf[0].cnt=0;
        allocReals(&rdf[0].bins,nBins);
	allocReals(&rdf[0].bins2,nBins);
        for(i=0;i<nBins;i++) 
	{
		rdf[0].bins[i]=0.0;
		rdf[0].bins2[i]=0.0;
	}
        return 0;
}

int getBoxVol(t_box box,real *vol)
{
        real aLen,bLen,cLen;
        real cosAB,sinAB;

        vecRNorm(box.a,&aLen);
        vecRNorm(box.b,&bLen);
        vecRNorm(box.c,&cLen);

        vecRProd(box.a,box.b,&cosAB);
        cosAB=cosAB/(aLen*bLen);
        sinAB=1-cosAB*cosAB;
        if(sinAB>=0.0) sinAB=sqrt(sinAB);
        else sinAB=0.0;

        vol[0]=sinAB*aLen*bLen*cLen;

        return 0;
}

int doRDF(t_atom *atoms,int nAtoms,t_grp *grps,t_rdf *rdf,t_box box)
{
        real tmp,inc;
        int i,j;
        int put;

        getBoxVol(box,&tmp);

        rdf[0].averDens=((real)grps[rdf[0].grp2Idx].nAtoms)/tmp;
        inc=1.0/rdf[0].averDens;

        for(i=0;i<grps[rdf[0].grp1Idx].nAtoms;i++)
        {
                for(j=0;j<grps[rdf[0].grp2Idx].nAtoms;j++)
                {
                        distPBC(atoms[grps[rdf[0].grp1Idx].atoms[i]].crd,atoms[grps[rdf[0].grp2Idx].atoms[j]].crd,box,&tmp);
                        put=(int)floor(tmp/rdf[0].binSize);
                        if(put<rdf[0].nBins)
                        {
                                rdf[0].bins[put]+=inc;
				rdf[0].bins2[put]+=1.0;
                        }
                }
                rdf[0].cnt++;
        }

        return 0;
}

int finalizeRDF(t_rdf *rdf)
{
        int i;
        real tmp;

        for(i=0;i<rdf[0].nBins;i++)
        {
                tmp=(4.0/3.0)*3.14159*(pow((i+1)*rdf[0].binSize,3.0)-pow(i*rdf[0].binSize,3.0));
		if(rdf[0].cnt!=0) rdf[0].bins[i]=rdf[0].bins[i]/(((real)rdf[0].cnt)*tmp);
                else rdf[0].bins[i]=0.0;
        }
	if(rdf[0].cnt!=0) rdf[0].bins2[0]/=rdf[0].cnt;
	else rdf[0].bins2[0]=0.0;
	for(i=1;i<rdf[0].nBins;i++)
        {
		if(rdf[0].cnt!=0) rdf[0].bins2[i]=rdf[0].bins2[i-1]+rdf[0].bins2[i]/rdf[0].cnt;
		else rdf[0].bins2[i]=0.0;
	}

        return 0;
}

