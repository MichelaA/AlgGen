/***************************************************************************
 *   Copyright (C) 2008 by miki   *
 *   michela.antonelli@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                               <usc[0                          *
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
#include "paes.h"

#include <fstream>

extern int numpieces,DIMBLOCK, NUMBLOCK, MINPART;
#include "ctime"


void paes::evaluatePop(double** inOutRed, int dimDt)
{	for (int i=0;i<arclengthV;i++)
		c[i].evaluateChrom(fis,inOutRed,dimDt);

}





paes::paes(int depth,int nsolution) //depth=4 nsolotion=2 fixed in the main program to this typical  values to implement the 2+2 Paes
{
	c=new chromosome[nsolution];

	for (int i=0;i<nsolution;i++)
		c[i].inizialize(fis);

	int numObj=c[0].getSizeObjAlg();

	gl_offset=new double[numObj];
	gl_range=new double[numObj];
	gl_largest=new double[numObj];
	for (int i=0;i<numObj;i++)
	{	gl_range[i]=0;
		gl_offset[i]=0;
	}

	RBcross=0;
	DBmut=0;
	RBmut=0;


	this->depth=depth;
	this->nsolution=nsolution;

	arclengthV=0; // current length of the virtual archive
	arclengthC=0; // current length of the cocrete archive

	iteration=Totaltrials/2;
	archive=popsize;

	//evaluating the current solutions vector
	for (int i=0;i<nsolution;i++)
		c[i].evaluateChrom(fis,inOutTr,numPatterTr);


	arcV=new chromosome[archive];

	for (int i=0;i<archive;i++)
		arcV[i].inizializewithoutrule();

	//Initializing the grid
	dimGrid=(int)pow((double)2,(numObj*depth))+1;  //grid dimension
	grid_pop=new int[dimGrid];
	for (int i=0;i<dimGrid;i++)
		grid_pop[i]=0;

	for (int i=0; i<nsolution;i++)
	{	update_grid(c[i],arcV,arclengthV);
		add_to_archive(c[i],arcV,arclengthV);
	}
	m=new chromosome[nsolution];
	for (int i=0;i<nsolution;i++)
		m[i].inizializewithoutrule();
}


paes::~paes(void)
{	delete[] grid_pop;
	delete[] c;
	delete[] arcV;
	delete[] m;
	delete[] gl_offset;
	delete[] gl_range;
	delete[] gl_largest;
}


void paes::evolvRed()
{
	double** bestInOut;
	sogaDT popDT(numPatterTr,DimChromR);

	chromosome* appo=new chromosome[archive];
	for (int i=0;i<archive;i++)
		appo[i].inizializewithoutrule();

	for (int i=0;i<arclengthV;i++) //valuto l'archivio considerando tutti i punti del dataset
	{	appo[i]=arcV[i];
		appo[i].evaluateChrom(fis,inOutTr,numPatterTr);
	}

	int numPoints=0;

	popDT.evaluate_pop(appo,arclengthV,fis,0);  //valuto la popolazione di dataset e scelgo il migliore
	dataset* bestPred(popDT.getDt(0));
	bestInOut=popDT.buildMatrix(*bestPred,numPoints); //costruisco la matrice dei dati ridotti con cui far evolvere la popolazione di frbs
	evaluatePop(bestInOut,numPoints);  //valuto l'archivio considerando il miglior dataset ridotto

	int num=iteration*nsolution/limVat;
	int iterSOGA = limVat*0.2/SIZE_POP_PR;
	int iterPaes=limVat/2;
	int salvait=num/48;
	if (salvait==0) salvait=1;
	int dimappo;
	for (int it=1;it<=num;it++)//num
	{	evolvPop(bestInOut,iterPaes,numPoints,it);	 //evolvo la popolazione di FRBS con bestPred e limVat valutazioni

	 	 //valuto l'archivio condierando tutti i punti del dataset
		dimappo=arclengthV;
		for (int i=0;i<arclengthV;i++)
		{	appo[i]=arcV[i];
			appo[i].evaluateChrom(fis,inOutTr,numPatterTr);
		}

		#ifdef _stampa
			if (it%10==0)
			{	save(appo,arcV,arclengthV,it*limVat);
				popDT.writePop(it*limVat);
			}
		#endif

		popDT.evol(appo, dimappo,fis,iterSOGA);	//evolvo il dataset

		bestPred=popDT.getDt(0);


		delete[] bestInOut;

		bestInOut=popDT.buildMatrix(*bestPred,numPoints);
		updateArc(bestInOut,numPoints,arcV[0].getSizeObjAlg());

	}

    //calcolo la fitness utilizzando tutto il dataset
	for (int i=0;i<arclengthV;i++)
	{	appo[i]=arcV[i];
		//if (CLASSIFICATION)
			appo[i].evaluateChrom(fis,inOutTr,numPatterTr);
		/*else
			appo[i].evaluateChromReg(fis,inOutTr,numPatterTr);*/
	}

	#ifdef _stampa
		write(appo,arclengthV,"PAES");
	#endif
}


void paes::updateArcDaChrom(chromosome* pop,int& dimpop,double** inOutRed, int numPatterRed)
{	// Resetting the grid of arcV

	int numObj=pop[0].getRealNumRule();
	for (int i=0;i<numObj;i++)
	{	gl_range[i]=0;
		gl_offset[i]=0;
	}
	for (int i=0;i<dimGrid;i++)
		grid_pop[i]=0;

	//Updating arcV
	chromosome* appo= new chromosome[dimpop];
	int dimAppo= dimpop;
	for (int i=0; i<dimAppo; i++)
	{	appo[i].inizializewithoutrule();
		appo[i]=pop[i];
	}
	dimpop=0;

	for (int i=0;i<dimAppo;i++)
	{	//appo[i].evaluateChrom(fis,inOutRed,numPatterRed);
		if (i==0 || compare_to_archive(appo[i],pop,dimpop)!=-1) 	//the new solution is no dominated with the archive. It is candidate to be added to the archive
		{	update_grid(appo[i],pop,dimpop);
			archive_soln(appo[i],pop,dimpop,inOutTr,numPatterTr);// the solution is added to the archive
        }
	}
	delete[] appo;
}

void paes::updateArc(double** inOutRed, int numPatterRed,int numObj)
{	// Resetting the grid of arcV

	for (int i=0;i<numObj;i++)
	{	gl_range[i]=0;
		gl_offset[i]=0;
	}
	for (int i=0;i<dimGrid;i++)
		grid_pop[i]=0;

	//Updating arcV
	chromosome* appo= new chromosome[arclengthV];
	int dimAppo= arclengthV;
	for (int i=0; i<dimAppo; i++)
	{	appo[i].inizializewithoutrule();
		appo[i]=arcV[i];
	}
	arclengthV=0;

	for (int i=0;i<dimAppo;i++)
	{	appo[i].evaluateChrom(fis,inOutRed,numPatterRed);
		if (i==0 || compare_to_archive(appo[i],arcV,arclengthV)!=-1) 	//the new solution is no dominated with the archive. It is candidate to be added to the archive
		{	update_grid(appo[i],arcV,arclengthV);
			archive_soln(appo[i],arcV,arclengthV,inOutTr,numPatterTr);// the solution is added to the archive
        }
	}
	delete[] appo;
}


/****************************************************************
*****************************************************************
Performs the genetic evolution of the population using
the predictor pred
*****************************************************************
****************************************************************/
/*void paes::evolvPop(double** inOutRed, int numEvol, int numPatterRed,int num)
{

	int* indPar;		// vector of the random indices for choosing current solutions
	double MutProbAtt;	//RB mutation probability (it depends if the RB crossover is applied or not (mutProbAtt=1) )
	int numVal;


	for (int it=1;it<=numEvol;it++)
	{
		numVal=it*nsolution;

		if (arclengthV>nsolution)
			indPar=RandintDistinct(0,arclengthV-1,nsolution);// estracting two random indeces for choosing current solutions
		else
		{	indPar=new int[nsolution];
			if (arclengthV==nsolution)
				for (int i=0; i<nsolution;i++)
					indPar[i]=i;
			else
				for (int i=0; i<nsolution;i++)
					indPar[i]=0;
		}

		for (int i=0; i<nsolution;i++)	// Copying the current solutions into the chromosome vector m to generate offspring
			m[i].copyChrom(arcV[indPar[i]]);
		for (int i=0; i<nsolution/2;i++)
			crossover(arcV[indPar[i*2]],arcV[indPar[i*2+1]],m[i*2],m[i*2+1],MutProbAtt); //(different crossovers are implemented in the function


		for (int z=0;z<nsolution;z++)
			mutation(m[z],MutProbAtt,mutProbPart);//(different mutations are implemented in the function

		for (int z=0;z<nsolution;z++) //Evaluating the new solutions and archive updating
		{	if (m[z].getNumRule()==0)
			{	printf("Error somewhere FRBS with zero rules\n");
				system("PAUSE");
			}
			m[z].evaluateChrom(fis,inOutRed,numPatterRed);

			if (compare_to_archive(m[z],arcV,arclengthV)!=-1) 	//the new solution is no dominated with the archive. It is candidate to be added to the archive
			{	update_grid(m[z],arcV,arclengthV);
				archive_soln(m[z],arcV,arclengthV,inOutRed,numPatterRed);// the solution is added to the archive

			}
        }
		delete[] indPar;
		qsort(arcV,arclengthV, sizeof(chromosome),&compareChrom);

	if (!(it%100))
		{ // for (int k=0;k<arclengthV;k++)
			stampaBest(arcV[0],it+(num-1)*numEvol,arclengthV);
			stampaBest(arcV[arclengthV-1],it+(num-1)*numEvol,arclengthV);
		}
		/*	if (!((it*nsolution)%50000))
			savestep(arcV,arclengthV,numVal,ProbVarPartion);
	}
}*/


/****************************************************************
*****************************************************************
Performs the genetic evolution
*****************************************************************
****************************************************************/
void paes::evolvPop(double** inOutRed, int numEvol, int numPatterRed,int num)
{
	int* indPar;// vector of the random indices for choosing current solutions
	double MutProbAtt;//RB mutation probability (it depends if the RB crossover is applied or not (mutProbAtt=1) )
	int numVal;

	if (numEvol==0)
		numEvol=iteration;

	for (int it=1;it<=numEvol;it++)
	{	numVal=it*nsolution;
		if (arclengthV>nsolution)
			indPar=RandintDistinct(0,arclengthV-1,nsolution);// estracting two random indeces for choosing current solutions
		else
		{	indPar=new int[nsolution];
			if (arclengthV==nsolution)
				for (int i=0; i<nsolution;i++)
					indPar[i]=i;
			else
				for (int i=0; i<nsolution;i++)
					indPar[i]=0;
		}

		for (int i=0; i<nsolution;i++)	// Copying the current solutions into the chromosome vector m to generate offspring
			m[i].copyChrom(arcV[indPar[i]]);

		for (int i=0; i<nsolution/2;i++)
			crossover(arcV[indPar[i*2]],arcV[indPar[i*2+1]],m[i*2],m[i*2+1],MutProbAtt); //(different crossovers are implemented in the function

		for (int z=0;z<nsolution;z++)
			mutation(m[z],MutProbAtt,probMutPart);//(different mutations are implemented in the function

		for (int z=0;z<nsolution;z++) //Evaluating the new solutions and archive updating
		{	if (m[z].getNumRule()==0)
			{	printf("Error somewhere FRBS with zero rules\n");

			}
			m[z].evaluateChrom(fis,inOutRed,numPatterRed);

			if (m[z].getRealNumRule()<minRule)
				m[z].copyChrom(arcV[indPar[z]]);

			if (compare_to_archive(m[z],arcV,arclengthV)!=-1) 	//the new solution is no dominated with the archive. It is candidate to be added to the archive
			{	update_grid(m[z],arcV,arclengthV);
				archive_soln(m[z],arcV,arclengthV,inOutRed,numPatterRed);// the solution is added to the archive
			}
        }
		delete[] indPar;
	    qsort(arcV,arclengthV, sizeof(chromosome),&compareChrom);


		#ifdef _stampa
			if (!(it%100))
				stampaBest(arcV[0],it+(num-1)*numEvol,arclengthV);

			if ((numPatterRed==numPatterTr) and !(numVal%5000))
				savestep(arcV, arclengthV,numVal,"NoD");

		#endif
	}// end principal loop

    #ifdef _stampa
		if (numPatterRed==numPatterTr)
		{	if (CLASSIFICATION)
				WriteResults("PAES",1-arcV[0].getObjTot(5),1-*arcV[0].ECM(fis,inOutTs,numPatterTs,true),arcV[0].getNumRule(),arcV[0].getComp());

			else
				WriteResults("PAES",arcV[0].getObjTot(5),*arcV[0].ECM(fis,inOutTs,numPatterTs,true),arcV[0].getNumRule(),arcV[0].getComp());

			saveresults(arcV,arclengthV);    //saving the results in a file
		}
		#endif
}

/*****************************************************************
******************************************************************
Compares a solution to every member of the archive.
Returns -1 if dominated by any member, 1 if dominates any member,
and 0 otherwise
	- m chromosome to be compared with the solutions in the archive
******************************************************************
******************************************************************/
int paes::compare_to_archive(chromosome& m, chromosome* arc,int arclength)
{
  int i=0;
  int result=0;
double* appo;

  while((i<arclength)&&(result!=1)&&(result!=-1))
  {		appo= arc[i].getObjAlg();
		result = compare_min(m.getObjAlg(), arc[i].getObjAlg(), arc[i].getSizeObjAlg());
      	i++;
   }
  return (result);
}



/*****************************************************************
******************************************************************
Update (if needed) the grid structures  and the location of the
archive solutions. Calculate the location of the chromosome s
	- s chromosome candidate to be insert into the archive
******************************************************************
******************************************************************/
void paes::update_grid(chromosome& s, chromosome* arc, int arclength)
{
  // recalculate ranges for grid in the light of a new solution s
	int square;
	int numObj=s.getSizeObjAlg();

	double* offset=new double[numObj];
	double* largest=new double[numObj];

	for (int i = 0; i < numObj; i++)
	{	offset[i] =  arc[0].getObjAlg(i);//minimun value for each objective
 		largest[i] = arc[0].getObjAlg(i);
	}
	for (int i = 0; i < numObj; i++)
  	{	for (int j = 1; j < arclength; j++)
		{	if (arc[j].getObjAlg(i) < offset[i])
				offset[i] = arc[j].getObjAlg(i);
			if (arc[j].getObjAlg(i) > largest[i])
				largest[i] = arc[j].getObjAlg(i);
        }
 	}
  	for (int b = 0; b < numObj; b++)
  	{	if (s.getObjAlg(b) < offset[b])
			offset[b] = s.getObjAlg(b);
    	if (s.getObjAlg(b) > largest[b])
			largest[b] = s.getObjAlg(b);
  	}

	if (find_loc(s.getObjAlg(),numObj)== dimGrid-1) // almeno un obiettivo e' fuori range
    	for (int a = 0; a < numObj; a++)
		{	gl_largest[a] = largest[a]+(0.2*largest[a]);
			gl_offset[a] = offset[a]-(0.2*offset[a]);
			gl_range[a] = gl_largest[a] - gl_offset[a];
		}

	for (int a = 0; a < dimGrid-1 ; a++)
		grid_pop[a] = 0;
	for (int a = 0; a < arclength; a++)
	{	square = find_loc(arc[a].getObjAlg(),numObj);
		arc[a].setGrid_loc(square);
		grid_pop[square]++;
	}
	square = find_loc(s.getObjAlg(),numObj);
	s.setGrid_loc(square);
	grid_pop[dimGrid-1] = -5;

	delete[] offset;
	delete[] largest;
}

/**************************************************************************
***************************************************************************
Given a solution s, add it to the archive if
   a) the archive is empty
   b) the archive is not full and s is not dominated or equal to anything
	 currently in the archive
   c) s dominates anything in the archive
   d) the archive is full but s is nondominated and is in a no more crowded
	square than at least one solution
 in addition, maintain the archive such that all solutions are nondominated.
****************************************************************************
****************************************************************************/
void paes::archive_soln(chromosome& s, chromosome* arc, int& arclength,double** matrPoints,int numPatt)
{
  	int i, repl=0, yes = 0, most, result, join = 0, old_arclength, set = 0;
  	int* tag=new int [archive];
	int* numeri=new int[archive];// contains the indices of the solution in the most crowded region

	int numObj=s.getSizeObjAlg();

	chromosome* tmp=new chromosome [archive];
	for (i=0;i<archive;i++)
    	tmp[i].inizializewithoutrule();
  	for (i = 0; i < archive; i++)
    	tag[i]=0;
   	if (arclength == 0)
  	{	 add_to_archive(s,arc,arclength);
      	delete[] tmp;
	  	delete[] tag;
		delete[] numeri;
	  	return;
  	}
  	i = 0;
  	result = 0;

  while((i < arclength)&&(result!=-1))
 {
      result = equal<double>(s.getObjAlg(), arc[i].getObjAlg(), numObj);
      if (result == 1)
		break;
      result = compare_min(s.getObjAlg(), arc[i].getObjAlg(), numObj);
      //  printf("%d\n", result);

    if ((result == 1)&&(join == 0))
	{
		arc[i].copyChrom(s);
	  	join = 1;
	}
    else if (result == 1)
	{
	  tag[i]=1;
	  set = 1;
	}
    i++;
 }

  old_arclength = arclength;
  if (set==1)
  {
    for (i = 0; i < arclength; i++)
	{
	  tmp[i].copyChrom(arc[i]);
	}
    arclength = 0;

    for (i = 0; i < old_arclength; i++)
	{
	  if (tag[i]!=1)
	    {
	      arc[arclength].copyChrom(tmp[i]);
	      arclength++;
	    }
	}
   }

  if ((join==0)&&(result==0))  // ie solution is non-dominated by the list
  	{	if (arclength == archive)
		{	most = grid_pop[s.getGrid_loc()]+1;
			for (i = 0; i < arclength; i++)
	    	{	if (grid_pop[arc[i].getGrid_loc()] >= most)
				{
			  		most = grid_pop[arc[i].getGrid_loc()];
			  		repl = i;
		  			yes = 1;
		  //   printf("i = %d\n", i);
				}
	    	}
	  		if (yes)
	  		{	/*if (s.getObj(1)<arc[repl].getObj(1))
		 		{*/	arc[repl].copyChrom(s);
          			grid_pop[s.getGrid_loc()]++;
				//}
	   		}
		}
    	else
			add_to_archive(s,arc,arclength);
	}

	if (arclength<2)
	{	arc[1].inizialize(fis);
		arc[1].evaluateChrom(fis,matrPoints,numPatt);
		arclength++;
	}




	delete[] numeri;
  	delete[] tmp;
 	delete[] tag;
}

/*****************************************************************
******************************************************************
Add the chromosome s to the archive updating the length
******************************************************************
******************************************************************/
void paes::add_to_archive(chromosome& s, chromosome* arc, int& arclength)
{
  arc[arclength].copyChrom(s);
  arclength++;
}


/*****************************************************************
******************************************************************
Finds the grid location of a solution given a vector of its objective
values
	- eval objectives vector
******************************************************************
******************************************************************/
int paes::find_loc(double * eval,int numObj)
{
  	int loc = 0;
	int n = 1;

	int *inc = new int[numObj];
	double *width = new double[numObj];
	double *gl_offsetApp = new double[numObj];

	for (int i = 0; i < numObj; i++)
    	gl_offsetApp[i]=gl_offset[i];

  // if the solution is out of range on any objective, return the maximum possible grid location number
	for (int i = 0; i < numObj; i++)
	{
    	if ((eval[i] < gl_offset[i])||(eval[i] > gl_offset[i] + gl_range[i]))
		{	delete[] inc;
			delete [] width;
			delete[] gl_offsetApp;
			return((int)pow((double)2,(numObj*depth)));
    	}
  	}

	for (int i = 0; i < numObj; i++)
	{
    	inc[i] = n;
		n *=2;
		width[i] = gl_range[i];
	}

	for (int d = 1; d <= depth; d++)
	{
    	for (int i = 0; i < numObj; i++)
		{
			if(eval[i] < width[i]/2+gl_offset[i])
	    		loc += inc[i];
			else
				gl_offsetApp[i] += width[i]/2;
		}
    	for (int i = 0; i < numObj; i++)
		{	inc[i] *= (numObj *2);
			width[i] /= 2;
		}
	}
	delete[] inc;
	delete [] width;
	delete[] gl_offsetApp;
	return(loc);
}



