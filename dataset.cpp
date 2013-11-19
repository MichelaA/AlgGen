/***************************************************************************
 *   Copyright (C) 2008 by miki   *
 *   michela.antonelli@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "dataset.h"

void dataset::inizialize(int dimChrom, int numBlock)
{
  Blocks=RandintDistinctOrd(0,numBlock-1,dimChrom);
  //int s;
/*  int k=dimChrom*0.4;
  for (int i=0;i<k;i++)
  {	  s=Randint(0,dimChrom-1);
	  Blocks[s]=-1;
  }*/
  this->dimChrom=dimChrom;

  nObj=2;
  obj=new double[2];

  obj[0]=0;
  obj[1]=56456465;


 /* for (int i=0;i<dimChrom;i++)
	  cout<<Blocks[i]<<' ';
  cout<<endl;*/


}

dataset::dataset(const dataset& dt)
{
	dimChrom=dt.dimChrom;
	Blocks=new int[dimChrom];
	for (int i=0;i<dimChrom;i++)
		Blocks[i]=dt.Blocks[i];
	obj=new double[dt.nObj];
	nObj=dt.nObj;
	for (int i=0;i<nObj;i++)
		obj[i]=dt.obj[i];

}

fstream& dataset::writeDT(fstream& fp,int num)
{

  fp<<num<<' ';
  /*for (int i=0;i<dimChrom;i++)
	fp<<Blocks[i]<<' ';*/
  for (int i=0;i<nObj;i++)
	fp<<obj[i]<<' ';
  fp<<endl;
  return fp;
}

int dataset::conta()
{	int aa=0;
	for (int i=0;i<dimChrom;i++)
		if (Blocks[i]==-1)
			aa++;

	return aa;


}
dataset& dataset::operator=(const dataset& dt)
{
	if(this!=&dt)
	{	Elimina();
		dimChrom=dt.dimChrom;
		Blocks=new int[dimChrom];
		for (int i=0;i<dimChrom;i++)
			Blocks[i]=dt.Blocks[i];
		obj=new double[dt.nObj];
		nObj=dt.nObj;
		for (int i=0;i<nObj;i++)
			obj[i]=dt.obj[i];
	}
	return *this;
}

bool dataset::operator==(const dataset& dt)
{
	if (dimChrom!=dt.dimChrom || nObj!=dt.nObj)
		return false;
	for (int i=0;i<dimChrom;i++)
		if (Blocks[i]!=dt.Blocks[i])
			return false;
	for (int i=0;i<nObj;i++)
		if (obj[i]!=dt.obj[i])
			return false;
	return true;
}


bool dataset::operator!=(const dataset& dt)
{
	if (dimChrom!=dt.dimChrom || nObj!=dt.nObj)
		return true;
	for (int i=0;i<dimChrom;i++)
		if (Blocks[i]!=dt.Blocks[i])
			return true;
	for (int i=0;i<nObj;i++)
		if (obj[i]!=dt.obj[i])
			return true;

	return false;
}

void dataset::Elimina()
{
	delete[] Blocks;
	delete[] obj;

}

dataset::~dataset()
{	Elimina();
}


void dataset::copyDataset(dataset& dt)
{
	dimChrom=dt.dimChrom;
	nObj=dt.nObj;
	for (int i=0;i<dimChrom;i++)
	  Blocks[i]=dt.Blocks[i];
	for (int i=0;i<nObj;i++)
	  obj[i]=dt.obj[i];
}


void cross(dataset& par1,dataset& par2,dataset& suc1,dataset& suc2,int Nblock)
{
	int num=Randint(1,par1.dimChrom-1);
	int i;
	//crossing the partition cromosome

	for (i=0;i<num;i++)
	{	suc1.Blocks[i]=par1.Blocks[i];
		suc2.Blocks[i]=par2.Blocks[i];
	}
	for (;i<par1.dimChrom;i++)
	{	suc1.Blocks[i]=par2.Blocks[i];
		suc2.Blocks[i]=par1.Blocks[i];
	}


	qsort(suc1.Blocks,par1.dimChrom, sizeof(int), &compareCre);
	qsort(suc2.Blocks,par1.dimChrom, sizeof(int), &compareCre);

	//check if there are duplicates in suc1 and suc2
	//bool trovato=false;
	for (int i=0;i<par1.dimChrom-1;i++)
		if (suc1.Blocks[i]==suc1.Blocks[i+1])
		{	suc1.Blocks[i]+=1;
			if (suc1.Blocks[i+1]==Nblock)
			{	for (int i=par1.dimChrom-1;i>0;i--)
				  suc1.Blocks[i]=suc1.Blocks[i-1];
				suc1.Blocks[0]=0;
				for (int i=0;i<par1.dimChrom-1 && suc1.Blocks[i]==suc1.Blocks[i+1];i++)
				  suc1.Blocks[i+1]+=1;
			}
		}


	for (int i=0;i<par1.dimChrom-1;i++)
		if (suc2.Blocks[i]==suc2.Blocks[i+1])
		{	suc2.Blocks[i]+=1;
			if (suc2.Blocks[i+1]==Nblock)
			{	for (int i=par1.dimChrom-1;i>0;i--)
				  suc2.Blocks[i]=suc2.Blocks[i-1];
				suc2.Blocks[0]=0;
				for (int i=0;i<par1.dimChrom-1 && suc2.Blocks[i]==suc2.Blocks[i+1];i++)
				  suc2.Blocks[i+1]+=1;
			}
		}




}


void dataset::Mutacion(int Nblock)
{

	int num=Randint(0,dimChrom-1);

	//double prob=Rand();

/*	if (prob<0.5)
	{	Blocks[num]=-1;
		return;
	}*/
	int nB, i;
	while(true)
	{	nB=Randint(0,Nblock-1);
		for ( i=0;i<dimChrom && nB!=Blocks[i];i++);
		if (i==dimChrom)
		{	Blocks[num]=nB;
			qsort(Blocks,dimChrom, sizeof(int), &compareCre);
			break;
		}
	}

}



int compareChromBlock(const  void * _a, const  void * _b)
{	dataset* elem1;
    dataset* elem2;

    elem1 =(dataset*) _a;
    elem2 = (dataset*) _b;
   if ( (*elem1).getObj(1)<(*elem2).getObj(1))
      return -1;
   else if ((*elem1).getObj(1)> (*elem2).getObj(1))
      return 1;
   else
      return 0;
}

