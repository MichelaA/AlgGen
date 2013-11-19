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
#include <fstream>
#include "sogaDT.h"



sogaDT::sogaDT(int numPTR,int dimch):EvolutionDT(dimch,numPTR)
{

	numElit=3;

	double nP=(double)(numPTR*Perc)/100;

	dimChrom=dimch;

	int pointXBlock=nP/dimChrom;
	int numBl=numPTR/pointXBlock;

	//=numPTR/numBl;
	//dimChrom=nP/pointXBlock;
	dimPopol=SIZE_POP_PR;

	parent_pop = new dataset[SIZE_POP_PR];
	child_pop = new dataset[SIZE_POP_PR];
	mixed_pop =new dataset[2*SIZE_POP_PR];

	inizializePop(parent_pop,SIZE_POP_PR,dimChrom,numBl); //nbl==massimo valore di un gene
	inizializePop(child_pop,SIZE_POP_PR,dimChrom,numBl);
	inizializePop(mixed_pop,2*SIZE_POP_PR,dimChrom,numBl);
}

sogaDT::~sogaDT()
{	delete[] parent_pop;
	delete[] child_pop;
	delete[] mixed_pop;
}

void sogaDT::evaluate_pop(chromosome* arcV, int archlen,frbs& fis,int type)
{	dataset* app;
     bool* flag=new bool[numBlock];
	double** errors=new double*[archlen]; //nell'ultima riga memorizzo la somma;
	double** errors1=new double*[archlen];
//	double** appo=new double*[archlen];
	int maxnum=0, indBlock;
	double err;
	double maxErr;
    switch (type) {
            case 0: app= parent_pop; maxnum=SIZE_POP_PR; break;
            case 1: app= child_pop; maxnum=SIZE_POP_PR; break;
            case 2: app= mixed_pop; maxnum=2*SIZE_POP_PR; break;
    }
	int apponum;
	double** mat=0;
	//int comp;
	double* dist=new double[numVar];
	for (int i=0;i<archlen;i++)
	{	errors[i]=new double[maxnum];
		errors1[i]=new double[maxnum];
		//appo[i]=new double[maxnum];
		for (int j=0;j<numBlock;j++)
			flag[j]=false;
		for (int k=0;k<maxnum;k++)
		{	errors[i][k]=0;
			errors1[i][k]=0;
			//appo[i][k]=0;
		}
		int realDimChrom;
		for (int k=0;k<maxnum;k++) //for each block in the chromosome the error is calculated
		{	apponum=0;
			realDimChrom=dimChrom;
			for (int j=0;j<dimChrom;j++)  //ciclo sui blocchi nel cromosoma
			{
			  indBlock=app[k].getBlock(j);
			  if (indBlock!=-1)
			  { 	apponum+=indTr[indBlock].dimBlock;
					if (!flag[indBlock])
					{	if (mat!=0)
							delete[] mat;

						mat=buildMatrix(indTr[indBlock].setpoint,indTr[indBlock].dimBlock);

					//	mat=buildMatrix(app[k],numPoints);
						indTr[indBlock].errBlock=*arcV[i].ECM(fis,mat,indTr[indBlock].dimBlock,false) ;

					//app[k].evaluateDataset(mat,indTr[indBlock].dimBlock,arcV[i],fis);
				//	cout<<indBlock<<' '<<indTr[indBlock].dimBlock<<endl;

						flag[indBlock]=true;
					}


				//maxErr=indTr[indBlock].errBlock < arcV[i].getObj(1)?arcV[i].getObj(1):indTr[indBlock].errBlock;
			//	if (maxErr!=0)
				errors[i][k]+=indTr[indBlock].errBlock;
				errors1[i][k]+=indTr[indBlock].errBlock;
				//appo[i][k]+=indTr[indBlock].errBlock;
			  }
			  else
				realDimChrom--;
			}
			errors1[i][k]/=realDimChrom;
			errors[i][k]/=realDimChrom;
		//	appo[i][k]/=dimChrom;



			double aa=arcV[i].getObjTot(5);



			/*	if (arcV[i].getNoCoverage())
			   aa/=apponum;*/
			maxErr=errors[i][k]< aa?aa:errors[i][k];
			//maxErr=errors[i][k]< arcV[i].getObj(1)?arcV[i].getObj(1):errors[i][k];
			errors1[i][k]=fabs(errors[i][k]-aa);
			if (maxErr!=0)
			   errors[i][k]=fabs(errors[i][k]-aa)/maxErr;
			else
			  errors[i][k]=fabs(errors[i][k]-aa);
			//errors[i][k]/=dimChrom;
		}


	}


	delete[] dist;
	/*for (int i=0;i<archlen;i++)
	{	for (int j=0;j<maxnum;j++)
			cout<<errors[i][j]<<' ';
		cout<<endl;

	}*/
	double err1;
	for (int i=0;i<maxnum;i++)
	{	err=0;
		err1=0;
		for (int j=0;j<archlen;j++)
		{	err+=errors[j][i];
			err1+=errors1[j][i];

		}
		err/=archlen;
		app[i].setObj(err,1);
		app[i].setObj(err1,0);
	}


	/*for (int j=0;j<archlen;j++)
	{	double aa=arcV[j].getObj(1);
		cout<<"errori "<<aa<<' '<<appo[j][0]<<' '<<endl;
		cout<<"diff "<<errors1[j][0]<<' '<<errors[j][0]<<' '<<endl;
	}*/

	for (int i=0;i<archlen;i++)
	{	delete[] errors[i];
	 	delete[] errors1[i];
	}
	delete[] errors;
	delete[] errors1;


}


int sogaDT::tournament(dataset& ind1, dataset& ind2)
{
    int flag;
    flag = compare_min (ind1.getObj(), ind2.getObj(),ind1.getnObj());

	if (flag==1)  // ind1 dominates ind2
	{	if (Rand()<=0.9)
			return 0;
		return 1;
	}
    if (flag==-1)  // ind2 dominates ind1
    {	if (Rand()<=0.9)
    		return 1;
    	return 0;
    }


	if (Rand() <= 0.5)		// random selection is performed
		return 0 ;

    return 1;

}

/****************************************************************
*****************************************************************
Generates the offspring population exploiting mutation and crossover
*****************************************************************
************************************\
***************************/
/*
void sogaDT::generate_child()
{
    int *a1, *a2;
    int temp,parent1,parent2;
    int rand;
    a1=new int[SIZE_POP_PR];
    a2=new int[SIZE_POP_PR];
	double MutProbAtt=0.2;
    for (int i=0; i<SIZE_POP_PR; i++)
		a1[i] = a2[i] = i;

    for (int i=0; i<SIZE_POP_PR; i++)// Generating random indices of parents to be selected for the binary tournament
    {
        rand = Randint (i, SIZE_POP_PR-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = Randint (i, SIZE_POP_PR-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    for (int i=0; i<SIZE_POP_PR; i+=4)
    {
		parent1 = tournament (parent_pop[a1[i]], parent_pop[a1[i+1]]);// Selection the first parent by binary tournament selection
        parent2 = tournament (parent_pop[a1[i+2]], parent_pop[a1[i+3]]);// Selection the second parent by binary tournament selection

		child_pop[i].copyDataset(parent_pop[a1[i+parent1]]);
		child_pop[i+1].copyDataset(parent_pop[a1[i+2+parent2]]);

		crossover(parent_pop[a1[i+parent1]], parent_pop[a1[i+2+parent2]], child_pop[i], child_pop[i+1]);

		for (int k=0; k<2; k++)
			mutation(child_pop[i+k],MutProbAtt);

		//second half
        parent1 = tournament (parent_pop[a2[i]], parent_pop[a2[i+1]]);// Selection the first parent by binary tournament selection
        parent2 = tournament (parent_pop[a2[i+2]], parent_pop[a2[i+3]]);// Selection the second parent by binary tournament selection

		child_pop[i+2].copyDataset(parent_pop[a2[i+parent1]]);
		child_pop[i+3].copyDataset(parent_pop[a2[i+2+parent2]]);

		crossover(parent_pop[a2[i+parent1]], parent_pop[a2[i+2+parent2]], child_pop[i+2], child_pop[i+3]);

       	for (int k=2; k<4; k++)
			mutation(child_pop[i+k],MutProbAtt);
	}

    delete[] a1;
    delete[] a2;

}

*/
int trova(double* a1)
{	double numero=Rand();
	int j=0;
	for (;j<SIZE_POP_PR && numero>a1[j] ;j++);
		//cout<<j<<' '<<a1[j]<<endl;
	if (j==SIZE_POP_PR)
		j--;
	return j;
}

void sogaDT::generate_child()
{
    int parent1,parent2;

    double * a1=new double[SIZE_POP_PR];

	double MutProbAtt=0.2;
	double sum=0;
    for (int i=0; i<SIZE_POP_PR; i++)
    	sum+=parent_pop[i].getObj(0);
    double sumProb=0;
    for (int i=0; i<SIZE_POP_PR; i++)
    {	if (i==0)
    		a1[i]=parent_pop[i].getObj(0)/sum;
    	else
    		a1[i]=a1[i-1]+parent_pop[i].getObj(0)/sum;
    }

    /*for (int i=0;i<SIZE_POP_PR;i++)
    	cout<<a1[i]<<' ';
    cout<<endl;*/
    double numero;
    for (int i=0; i<SIZE_POP_PR; i+=4)
    {

		parent1=trova(a1);
		do
			parent2=trova(a1);
		while (parent1==parent2);
	//	cout<<parent1<<' '<<parent2<<endl;
    	//parent1 = tournament (parent_pop[a1[i]], parent_pop[a1[i+1]]);// Selection the first parent by binary tournament selection
        //parent2 = tournament (parent_pop[a1[i+2]], parent_pop[a1[i+3]]);// Selection the second parent by binary tournament selection

		child_pop[i].copyDataset(parent_pop[parent1]);
		child_pop[i+1].copyDataset(parent_pop[parent2]);

		crossover(parent_pop[parent1], parent_pop[parent2], child_pop[i], child_pop[i+1]);

		for (int k=0; k<2; k++)
			mutation(child_pop[i+k],MutProbAtt);

		//second half
		//parent1 = tournament (parent_pop[a2[i]], parent_pop[a2[i+1]]);// Selection the first parent by binary tournament selection
        //parent2 = tournament (parent_pop[a2[i+2]], parent_pop[a2[i+3]]);// Selection the second parent by binary tournament selection
        parent1=trova(a1);
        parent1=trova(a1);
     //   cout<<parent1<<' '<<parent2<<endl;
        child_pop[i+2].copyDataset(parent_pop[parent1]);
		child_pop[i+3].copyDataset(parent_pop[parent2]);

		crossover(parent_pop[parent1], parent_pop[parent2], child_pop[i+2], child_pop[i+3]);

       	for (int k=2; k<4; k++)
			mutation(child_pop[i+k],MutProbAtt);
	}

    delete[] a1;


}


/****************************************************************
*****************************************************************
Evolution of the traditional soga
*****************************************************************
****************************************************************/
void sogaDT::evol(chromosome* arcV, int archlen, frbs& fis, int gener)
{
	//dataset* appo;


	evaluate_pop(arcV,archlen,fis,0); //Evaluating the initial parent population
    qsort(parent_pop,SIZE_POP_PR, sizeof(dataset),&compareChromBlock);



  /*  for (int i=0;i<SIZE_POP_PR;i++)
    {	appo=getDt(i);
		cout<<appo->conta()<<endl;

    	//appo=getDt(i);
		//appo->writeDT((fstream&)cout,1);
    }*/



    cout<<endl<<" Generations starting"<<endl;

	//int numVal;
    for (int i=1; i<=gener; i++)
    {   generate_child();//generating the offspring population
        evaluate_pop(arcV, archlen,fis,1); // evaluating the offspring population
        qsort(child_pop,SIZE_POP_PR, sizeof(dataset),&compareChromBlock);
        extractPopulation();
	   /*appo=getDt(0);
	   appo->writeDT((fstream&)cout,1);*/
    }

    /*for (int i=0;i<SIZE_POP_PR;i++)
       {	appo=getDt(i);
   		appo->writeDT((fstream&)cout,1);
       }
	cout<<endl<<" Generations finished, now reporting solutions"<<endl;
*/
   /* for (int i=0;i<SIZE_POP_PR;i++)
       {	appo=getDt(i);
   		cout<<appo->conta()<<endl;*/

       	//appo=getDt(i);
   		//appo->writeDT((fstream&)cout,1);
      // }

    /*appo=getDt(0);
    cout<<appo->conta()<<endl;*/
}


void sogaDT::writePop(int num)
{
	char nomefile[100]="EvolDataset.txt";
	fstream fpObj;

	fpObj.open(nomefile,ios::out | ios::app);
	if (!fpObj)
	{	cout<<("Error opening file\n");
       	exit(1);
    }
	dataset* appo;

	/*for (int i=0;i<dimPopol;i++)
	{ */   appo=getDt(0);
	     appo->writeDT(fpObj,num);
//	}
  fpObj<<endl;

}


/****************************************************************
*****************************************************************
Updates the current population substituting the worst numElit of
offspring  with the best numElit of the parent population
*****************************************************************
****************************************************************/
void sogaDT::extractPopulation (){

  /* for (int i=0;i<numElit;i++)
		parent_pop[SIZE_POP_PR-1-i].copyDataset(child_pop[i]);*/
      for (int i=numElit;i<SIZE_POP_PR;i++)
		parent_pop[i].copyDataset(child_pop[i-numElit]);
	  qsort(parent_pop,SIZE_POP_PR, sizeof(dataset),&compareChromBlock);
}



