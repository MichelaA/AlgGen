/*
    Copyright (c) <year>, <copyright holder>
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the <organization> nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY <copyright holder> ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL <copyright holder> BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "evolutionDT.h"


EvolutionDT::EvolutionDT(int dimch,int numPTR)
{

	double nP=(double)(numPTR*Perc)/100;

	dimChrom=dimch;

	int pointXBlock=nP/dimChrom;
	numBlock=numPTR/pointXBlock;



	//int dimBlock=numPTR/nBl;
	indTr=new block[numBlock];
	int* index=RandintDistinct(0,numPTR-1,numPTR);
	int indice=0;

	for (int i=0;i<numBlock-1;i++)
	{	indTr[i].dimBlock=pointXBlock;
		indTr[i].setpoint=new int[indTr[i].dimBlock];
		for (int j=0;j<indTr[i].dimBlock;j++)
		{	indTr[i].setpoint[j]=index[indice];
			indice++;
		}
	}
	indTr[numBlock-1].dimBlock=numPTR/numBlock+numPTR%numBlock;
		indTr[numBlock-1].setpoint=new int[indTr[numBlock-1].dimBlock];
	for (int j=0;j<indTr[numBlock-1].dimBlock;j++)
	{	indTr[numBlock-1].setpoint[j]=index[indice];
		indice++;
	}

	delete[] index;
}



EvolutionDT::~EvolutionDT()
{
  for (int i=0;i<numBlock;i++)
		delete[] indTr[i].setpoint;
	delete[] indTr;


}


double** EvolutionDT::buildMatrix(dataset& dt, int& numPoints)
{
	int* indBlock=dt.getBlocks();
	numPoints=0;
	for (int i=0;i<dt.getDim();i++)
		numPoints+=indTr[indBlock[i]].dimBlock;

	double** mat=new double*[numPoints];
	int k=0;
	for (int i=0;i<dimChrom;i++)
		for (int j=0;j<indTr[indBlock[i]].dimBlock;j++)
			mat[k++]=inOutTr[indTr[indBlock[i]].setpoint[j]];
	return mat;
}


double** EvolutionDT::buildMatrix(int* blocco, int dimblocco)
{
	double** mat=new double*[dimblocco];
	for (int i=0;i<dimblocco;i++)
			mat[i]=inOutTr[blocco[i]];
	return mat;
}




//inizialize the population of chromosomes
void EvolutionDT::inizializePop(dataset* dt, int dimP,int dimChrom ,int numBlock)
{
	for (int i=0;i<dimP;i++)
	  dt[i].inizialize(dimChrom,numBlock);
}
/******************************************************************
*******************************************************************
Generates two offspring crossing all the chromosome part
	-	p0 first parent to be crossed
    -	p1 second parent to be crossed
	-	f0 first offspring generated by crossover
	-	f1 second offspring generated by crossover
 	-	MutProbAtt returns the value of the RB mutation probability (it depends if the RB crossover is applied or not (mutProbAtt=1) )
******************************************************************
*******************************************************************/
void EvolutionDT::crossover(dataset& p0, dataset& p1, dataset& f0, dataset& f1)
{
	if (Rand()< crossProb)
	  cross(p0,p1,f0,f1,numBlock);

}


/******************************************************************
*******************************************************************
Permforms the mutation operation
	-	 m return the muteted chromosome
    -	MutProbAtt RB mutation probability
	-	MutProbGRAtt Granularity mutation probability
******************************************************************
*******************************************************************/
void EvolutionDT::mutation(dataset& m,double MutProb)
{
	if (Rand()<MutProb)
	  m.Mutacion(numBlock);


}

/******************************************************************
*******************************************************************
Save the results of the population on files
     -  chrom vector of chromosomes to be saved
    -	popsize dimension of the population
	- 	conc if equal to 1 save the results of the concrete RB
		otherwise save the results of the virtual RB
******************************************************************
*******************************************************************/
/*void evolution::saveresults(chromosome* chrom, int popsize, bool conc)
{
	fstream fpObj,fpM, fpGR, fpT;
	double app;
	int noZero;
	int** rule;


	if (!conc)
	{	fpObj.open("NoDomV.txt",ios::out | ios::app);	//Containes for each solution in the pareto a row with the complexity, the training MSE,
		fpM.open("myMatV.txt",ios::out | ios::app);  	//Containes for each solution in the pareto the virtual RB
		fpT.open("TunV.txt",ios::out | ios::app);
		fpGR.open("GranV.txt",ios::out | ios::app);
	}
	else							//the integrity, the test MSE, the concrete and the virtual number of rules and the actual
	{	fpObj.open("NoDomC.txt",ios::out | ios::app);						//number of variable used
		fpM.open("myMatC.txt",ios::out | ios::app);  	//Containes for each solution in the pareto the virtual RB  	//Containes for each solution in the pareto the virtual RB
		fpT.open("TunC.txt",ios::out | ios::app);
		fpGR.open("GranC.txt",ios::out | ios::app); 		//Containes for each solution in the pareto the granularity of each variable
	}

	if (!fpObj || !fpM || !fpT || !fpGR)
	{	printf("Error opening file\n");
       	exit(1);
    }

	for (int i=0;i<6;i++)
	{	fpObj<<"0 ";
		fpT<<"0 ";
		fpM<<"0 ";
		fpGR<<"0 ";
	}


	double comp=0,coovTs=0;
	int* appopart;
	float* pieces=new float[numpieces];
	int objectives;

	for (int i=0;i<popsize;i++)
	{	fpObj<<endl;
		fpT<<endl;
		fpM<<endl<<i+1<<')'<<endl;
		fpGR<<endl;

		objectives=chrom[i].getSizeObj();
		app=chrom[i].ECM_ts(fis,coovTs,conc);
		noZero=chrom[i].giveActive();
		rule=chrom[i].decode();
		int numR;
		if(!conc) //scrivo il virtuale
			numR=chrom[i].getNumRule();
		else
		{	chrom[i].convertimat (rule,fis);
			numR=chrom[i].getRealNumRule();
		}
		for (int k=0;k<numR;k++)
		{	for(int j=0;j<numVar;j++)
				fpM<<rule[k][j]<<'\t';
			fpM<<endl;
		}

		for (int k=0;k<objectives;k++)
			if (k==2)
				fpObj<<1-(2*chrom[i].getObj(k))/(MAXPART-2)<<' ';
			else
				fpObj<<chrom[i].getObj(k)<<' ';
		fpObj<<app<<' ';
		fpObj<<numR<<' ';
		fpObj<<noZero<<' ';
	//	fpObj<<1-(2*chrom[i].getInterp())/(MAXPART-2)<<' ';
		appopart=chrom[i].getRealPart();

		if (ProbVarPartion!=0)
			for (int j=0;j<numVar;j++)
			{	fpGR<<appopart[j]<<' ';
				fpObj<<appopart[j]<<' ';
			}
		for (int j=0;j<numVar;j++)
		{	if (TuningSF!=0)
				fpT<<chrom[i].getLambda_i(j)<<' '<<chrom[i].getKappa_i(j)<<' ';
			if (TuningPW!=0)
			{	pieces=chrom[i].getPwLT(j);
				for (int k=0;k<numpieces;k++)
				{	fpT<<pieces[k]<<' ';
					fpObj<<pieces[k]<<' ';
				}
			}
		}
		//fpGR<<endl;
		//fpT<<endl;

		if (Tuning2T)
		{	float * ap=chrom[i].getTun();
			for (int j=0;j<chrom[i].getNumPar();j++)
				fpT<<ap[i]<<' ';
      	}
		deleteMat(rule,numR);
	}

	fpObj<<endl;
	fpM<<endl;
	fpGR<<endl;
	fpT<<endl;
}




/******************************************************************
*******************************************************************
Print the solution with best accuracy
	-	arc0 vchromosome with the best accuracy
    -	it number of iteration
	-	archlenght archive dimension
*******************************************************************
*******************************************************************/
/*void evolution::stampaBest(chromosome& arc0,int it,int arclength,bool conc)
{
	float* pieces=new float[numpieces];
	int numR,numP=arc0.getNumPar();
	float* tuple=new float[numP];
	double d;

	cout<<it<<' ';
	for (int i=0;i<arc0.getSizeObj();i++)
		if (i==2)
			cout<<1-(2*arc0.getObj(i))/(MAXPART-2)<<' ';
		else
			cout<<arc0.getObj(i)<<' ';
	if (conc)
		numR=arc0.getRealNumRule();
	else
		numR=arc0.getNumRule();
	cout<<arc0.ECM_ts(fis,d,conc)<<' '<<numR<<' ';

	//cout<<1-(2*arc0.getInterp())/(MAXPART-2)<<' ';
	for (int j=0;j<numVar;j++)
		cout<<arc0.getRealPart(j)<<' ';
	for (int j=0;j<numVar;j++)
	{	if (TuningSF!=0)
			cout<<arc0.getLambda_i(j)<<' '<<arc0.getKappa_i(j)<<' ';
		if (TuningPW!=0)
		{	pieces=arc0.getPwLT(j);
			for (int k=0;k<numpieces;k++)
				cout<<pieces[k]<<' ';
		}
		if (Tuning2T!=0)
		{	tuple=arc0.getTun();
			for (int i=0;i<numP;i++)
				cout<<tuple[i]<<' ';
		}
	}
	cout<<"dim_arc:"<<arclength<<endl;
}
*/

/******************************************************************
*******************************************************************
Save the population after numEval evaluationed
    -	popsize dimension of the population
	-	numEval number of evaluation
*******************************************************************
*******************************************************************/
/*
void evolution::savestep(chromosome* chrom, int popsize,int numEval,bool conc)
{
	int noZero,numR;
	char numst[100];
	char nomefile[100]="Nd";
	char nomefileGr[100]="Gr";
	char nomefileTu[100]="Tu";
	char nomefileMat[100]="Mat";


	fstream fpObj,fpM, fpGR, fpT;
	double app;
	int** rule;


	sprintf(numst,"%i",numEval);
	strcat(nomefile,numst);
	strcat(nomefileGr,numst);
	strcat(nomefileTu,numst);
	strcat(nomefileMat,numst);


	fpObj.open(nomefile,ios::out | ios::app);
	fpM.open(nomefileMat,ios::out | ios::app);
	fpT.open(nomefileTu,ios::out | ios::app);
	fpGR.open(nomefileGr,ios::out | ios::app);

	if (!fpObj || !fpM || !fpT || !fpGR)
	{	printf("Error opening file\n");
       	exit(1);
    }
	for (int i=0;i<6;i++)
	{	fpObj<<"0 ";
		fpT<<"0 ";
		fpM<<"0 ";
		fpGR<<"0 ";
	}

	double comp=0,coovTs=0;
	int* appopart;
	float* pieces=new float[numpieces];
	int objectives;

	for (int i=0;i<popsize;i++)
	{	fpObj<<endl;
		fpT<<endl;
		fpM<<endl<<i+1<<')'<<endl;
		fpGR<<endl;

		objectives=chrom[i].getSizeObj();
		app=chrom[i].ECM_ts(fis,coovTs,conc);
		noZero=chrom[i].giveActive();
		rule=chrom[i].decode();
		int numR;
		if(!conc) //scrivo il virtuale
			numR=chrom[i].getNumRule();
		else
		{	chrom[i].convertimat (rule,fis);
			numR=chrom[i].getRealNumRule();
		}
		for (int k=0;k<numR;k++)
		{	for(int j=0;j<numVar;j++)
				fpM<<rule[k][j]<<'\t';
			fpM<<endl;
		}

		for (int k=0;k<objectives;k++)
			if (k==2)
				fpObj<<1-(2*chrom[i].getObj(k))/(MAXPART-2)<<' ';
			else
				fpObj<<chrom[i].getObj(k)<<' ';
		fpObj<<app<<' ';
		fpObj<<numR<<' ';
		fpObj<<noZero<<' ';
		//fpObj<<1-(2*chrom[i].getInterp())/(MAXPART-2)<<' ';

		appopart=chrom[i].getRealPart();

		if (ProbVarPartion!=0)
			for (int j=0;j<numVar;j++)
			{	fpGR<<appopart[j]<<' ';
				fpObj<<appopart[j]<<' ';
			}
		for (int j=0;j<numVar;j++)
		{	if (TuningSF!=0)
				fpT<<chrom[i].getLambda_i(j)<<' '<<chrom[i].getKappa_i(j)<<' ';
			if (TuningPW!=0)
			{	pieces=chrom[i].getPwLT(j);
				for (int k=0;k<numpieces;k++)
				{	fpT<<pieces[k]<<' ';
					fpObj<<pieces[k]<<' ';
				}
			}
		}
		//fpGR<<endl;
		//fpT<<endl;

		if (Tuning2T)
		{	float * ap=chrom[i].getTun();
			for (int j=0;j<chrom[i].getNumPar();j++)
				fpT<<ap[i]<<' ';
      	}
		deleteMat(rule,numR);
	}

	fpObj<<endl;
	fpM<<endl;
	fpGR<<endl;
	fpT<<endl;




/*


	fp.open(nomefile,ios::out | ios::app);
   	if(!fp)
	{	printf("Error opening file\n");
		exit(1);
     }
	fp<<endl;
    for (int i=0;i<5;i++)
		fp<<"0 ";
		fp<<endl;

    for (int i=0;i<popsize;i++)
	{   noZero=chrom[i].giveActive();
		int numP=chrom[i].getNumPar();
		for (int k=0;k<chrom[i].getSizeObj();k++)
			if (k==2)
				fp<<1-(2*chrom[i].getObj(k))/(MAXPART-2)<<' ';
			else
				fp<<chrom[i].getObj(k)<<' ';
		if (conc)
			numR=chrom[i].getRealNumRule();
		else
			numR=chrom[i].getNumRule();

		fp<<chrom[i].ECM_ts(fis,coovTs,conc)<<' '<<numR<<' '<<noZero<<' ';
		appopart=chrom[i].getRealPart();
		for (int j=0;j<numVar;j++)
			fp<<appopart[j]<<' ';


		float* pieces;
		for (int j=0;j<numVar;j++)
		{	if (TuningSF!=0)
				fp<<chrom[i].getLambda_i(j)<<' '<<chrom[i].getKappa_i(j)<<' ';
			if (TuningPW!=0)
			{	pieces=chrom[i].getPwLT(j);
				for (int k=0;k<numpieces;k++)
					fp<<pieces[k]<<' ';
			}
			float* tuple;
			if (Tuning2T!=0)
			{	tuple=chrom[i].getTun();
				for (int i=0;i<numP;i++)
					fp<<tuple[i]<<' ';
			}
		}
		fp<<endl;
	}
	*/
//}


/******************************************************************
*******************************************************************
Sort the chromosome population by ascendent accuracy
    -	chrom population of chromosomes
	-	dimPop population dimension
*******************************************************************
*******************************************************************/
/*void evolution::ordinaPopChrom(chromosome* chrom,int dimPop)
{    qsort(chrom, dimPop, sizeof(chromosome),&compareChrom);
}*/










