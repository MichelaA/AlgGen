/*
 * soga.cpp
 *
 *  Created on: Sep 28, 2012
 *      Author: miki
 */
//#include <map>
#include "soga.h"

	double** soga::oldTr=0;
	double** soga::oldTs=0;
	int soga::oldMaxterms=0;
	int soga::oldNumVar=0;
	int* soga::oldNumParts=0;
	double* soga::oldMinVal=0;
	double * soga::oldMaxVal=0;
	//int soga::indMat=0;
	//map<string,soga::matGran>
	map<string,soga::matGran> soga::matSal;




void soga::inizializePop(crom* p,int n)
{	for (int i=0;i<n;i++)
	{	p[i].genera(numGen);
		p[i].setGen(numParts);
	}
}


soga::soga() {

	numPop=SIZE_POP_PR;
	numGen=numVar-1;
	parent_pop = new crom[numPop];
	child_pop = new crom[numPop];

	//mixed_pop = new crom[2*numPop];

	inizializePop(parent_pop,numPop);
	inizializePop(child_pop,numPop);
	//inizializePop(mixed_pop,2*numPop);

	if (oldTr==0)
	{	oldTr=new double*[numPatterTr];
		oldTs=new double*[numPatterTs];
		for (int i=0;i<numPatterTr;i++)
		{	oldTr[i]=new double[numVar];
			if (i<numPatterTs)
				oldTs[i]=new double[numVar];
		}
		oldNumParts=new int[numVar];
		oldMinVal=new double[numVar];
		oldMaxVal=new double[numVar];

	}
}



/*soga::soga(const soga& s)
{	numPop=s.numPop;
	numGen=s.numGen;
	popol=new crom[numPop];
	for (int i=0;i<numPop;i++)
	{	popol[i].gen=new int[numGen];
		for (int j=0;j<numGen;j++)
		{	popol[i].gen[j]=s.popol[i].gen[j];
			popol[i].fitness=s.popol[i].fitness;
		}
	}
}*/

soga::~soga() {
	/*for (int i=0;i<numPop;i++)
	{	//cout<<parent_pop[i].gen<<endl;
		//delete[] parent_pop[i].gen;
	//	delete[] child_pop[i].gen;
	}
	/*for (int i=0;i<2*numPop;i++)
		delete[] mixed_pop[i].gen;*/

	delete[] parent_pop;
	delete[] child_pop;
	//delete[] mixed_pop;

	for (int i=0;i<numPatterTr;i++)
	{	delete[] oldTr[i];
		if (i<numPatterTs)
			delete[] oldTs[i];

	}
	delete[] oldTr;
	delete[] oldTs;
	delete[] oldNumParts;
	delete[] oldMinVal;
	delete[] oldMaxVal;
}

void soga::ordinaPopol(crom* popol, int n) ///< Orders a vector in ascendent way
{
	bool ordinato=false;
	crom appo(numGen);

	for (int i=0;i<n-1 && !ordinato ;i++)
	{	ordinato = true;
		for (int j = n-1; j >= i+1; j--)
			if(	popol[j].getFitness() < popol[j-1].getFitness())
			{	appo=popol[j];
				popol[j]=popol[j-1];
				popol[j-1]=appo;
				ordinato = false;
			}
	}


}

void soga::generate_child()
{
    int *a1, *a2;
    int temp,parent1,parent2;
    int rand;
    a1=new int[numPop];
    a2=new int[numPop];
	double MutProbAtt;
    for (int i=0; i<numPop; i++)
		a1[i] = a2[i] = i;

    for (int i=0; i<numPop; i++)// Generating random indices of parents to be selected for the binary tournament
    {
        rand = Randint (i, numPop-1);
        temp = a1[rand];
        a1[rand] = a1[i];
        a1[i] = temp;
        rand = Randint (i, numPop-1);
        temp = a2[rand];
        a2[rand] = a2[i];
        a2[i] = temp;
    }
    for (int i=0; i<numPop; i+=4)
    {
		parent1 = tournament (parent_pop[a1[i]], parent_pop[a1[i+1]]);// Selection the first parent by binary tournament selection
        parent2 = tournament (parent_pop[a1[i+2]], parent_pop[a1[i+3]]);// Selection the second parent by binary tournament selection

        child_pop[i] = parent_pop[a1[i+parent1]];
		child_pop[i+1]=parent_pop[a1[i+2+parent2]];

		onePointCross(parent_pop[a1[i+parent1]], parent_pop[a1[i+2+parent2]], child_pop[i], child_pop[i+1]);

		for (int k=0; k<2; k++)
			mutation(child_pop[i+k]);

		//second half
        parent1 = tournament (parent_pop[a2[i]], parent_pop[a2[i+1]]);// Selection the first parent by binary tournament selection
        parent2 = tournament (parent_pop[a2[i+2]], parent_pop[a2[i+3]]);// Selection the second parent by binary tournament selection
		child_pop[i+2]=parent_pop[a2[i+parent1]];
		child_pop[i+3]=parent_pop[a2[i+2+parent2]];
		onePointCross(parent_pop[a2[i+parent1]], parent_pop[a2[i+2+parent2]], child_pop[i+2], child_pop[i+3]);

       	for (int k=2; k<4; k++)
       		mutation(child_pop[i+k]);
	}



    /* for (int i=0;i<numPop;i++)
    		cout<<parent_pop[i].gen<<endl;
    cout<<endl;*/

 /*   for (int i=0;i<numPop;i++)
        		cout<<child_pop[i].gen<<endl;*/

    delete[] a1;
    delete[] a2;

}

int soga::tournament(crom& ind1, crom& ind2)
{
    //int flag;
    //flag = compare_min (ind1.getObjAlg(), ind2.getObjAlg(),ind1.getSizeObjAlg());

	if (ind1.getFitness()<ind2.getFitness())  // ind1 dominates ind2
		return 0;
    else
        return 1;
    /*if (flag==-1)  // ind2 dominates ind1
    {*/

    //}

}

void soga::extractPopulation ()
{
	ordinaPopol(child_pop,numPop);
   for (int i=0;i<10;i++)
		parent_pop[numPop-1-i]=child_pop[i];

   /*for (int i=0;i<numPop;i++)
   		cout<<parent_pop[i].gen<<endl;*/

   ordinaPopol(parent_pop,numPop);




}

void soga::evolution()
{
	evaluatePop(parent_pop);
	ordinaPopol(parent_pop,numPop);
	cout<<"Initial population evalutated"<<endl;
	int* app;
	for (int i=2; i<=50; i++)
	 {	//numVal=i*popsize;
	    //cout<<i<<endl;
		generate_child();//generating the offspring population
	    extractPopulation();
	    evaluatePop(parent_pop);  // evaluating the offspring population
	    ordinaPopol(parent_pop,numPop);
   		if (!(i%10))
   		{
   			for (int j=0;j<numPop;j++)
			{	cout<<"Gran: ";
				app=parent_pop[j].getGen();
				for (int k=0;k<numVar-1;k++)
					cout<<app[k]<<' ';
				cout<<endl<<"Fitness: "<<parent_pop[j].getFitness()<<' '<<"NumRegole: "<<parent_pop[j].getNregole()<<' '<<"Errore: "<<(parent_pop[j].getFitness()-0.2*parent_pop[j].getNregole())/0.8<<endl;

			}
   			cout<<endl<<endl;
   		}

		//if (!(numVal%25000))
		//savestep(parent_pop,popsize,numVal,"nd");
	 }
	cout<<endl<<" Generations finished, now reporting solutions"<<endl;
	for (int j=0;j<numPop;j++)
	{	app=parent_pop[j].getGen();
		cout<<"Gran: ";
		for (int k=0;k<numVar-1;k++)
			cout<<app[k]<<' ';
		cout<<endl<<" Fitness: "<<parent_pop[j].getFitness()<<" NumRegole: "<<parent_pop[j].getNregole()<<endl; //" Errore: "<<(parent_pop[j].getFitness()-0.2*parent_pop[j].getNregole())/0.8<<endl;
	}
}

void soga::copiaValori(double** Tr,double** Ts,int NVar,int Maxt,int* NParts,double* mVal,double* mxVal)
{
	oldNumVar=NVar;
	oldMaxterms=Maxt;
	for (int i=0;i<numVar;i++)
	{	oldMinVal[i]=mVal[i];
		oldMaxVal[i]=mxVal[i];
		oldNumParts[i]=NParts[i];
	}

	for (int i=0;i<numPatterTr;i++)
		for (int j=0;j<NVar;j++)
		{	oldTr[i][j]=Tr[i][j];
			if (i<numPatterTs)
				oldTs[i][j]=Ts[i][j];
		}

}

void soga::setc45Best()
{	int* indici=0;
	int numAtt=0;
	int* app=parent_pop[0].getGen();
	map<string,matGran>::iterator it;
	for (int k=0;k<numGen;k++)
		numParts[k]=app[k];
	matC45=calcolaMatC45(dimmatC45,numAtt,indici);
	cambiaVariabili(numAtt,indici,numParts);
	delete[] indici;
	for (it=matSal.begin();it!=matSal.end();it++)
		for (int i=0;i<it->second.nRighe;i++)
			delete[] it->second.mat[i];
		delete[] it->second.mat;


}



void soga::ripristinaValori()
{
	numVar=oldNumVar;
//numVar-1=oldMaxterms;
	for (int i=0;i<numVar;i++)
	{	minVal[i]=oldMinVal[i];
		maxVal[i]=oldMaxVal[i];
		numParts[i]=oldNumParts[i];
	}

	for (int i=0;i<numPatterTr;i++)
		for (int j=0;j<numVar;j++)
		{	inOutTr[i][j]=oldTr[i][j];
			if (i<numPatterTs)
				inOutTs[i][j]=oldTs[i][j];
		}

}

string soga::faiStringa(int* vett)
{	ostringstream oss(ostringstream::out);
	for (int temp = 0; temp < numGen; temp++)
		oss << vett[temp];
	return oss.str();
}

void soga::evaluatePop(crom* popol)
{
	int** matrice=0;
	int dimMat;
	int* indici;
	int numAtt;
	string codice;
	matGran appo;
	//int comp;
	int* app;
	map<string,matGran>::iterator it;
	copiaValori(inOutTr,inOutTs,numVar,numVar-1,numParts,minVal,maxVal);
	for (int i=0;i<numPop;i++)
	{
		app=popol[i].getGen();
		codice=faiStringa(app);
		it=matSal.find(codice);
		if (it==matSal.end())  //non c'e' l'entrata corrispondente devo calcolare la matrice con C4.5 e l'errore
		{	app=parent_pop[i].getGen();
			for (int k=0;k<numGen;k++)
				numParts[k]=app[k];
			matrice=calcolaMatC45(dimMat,numAtt,indici);
			cambiaVariabili(numAtt,indici,numParts);
			parent_pop[i].setFitness(valutamatrice(matrice, inOutTr,dimMat)*100*0.8+0.2*dimMat);
			parent_pop[i].setNregole(dimMat);
			appo.fitness=parent_pop[i].getFitness();
			appo.mat=matrice;
			appo.nRighe=dimMat;
			appo.ncolonne=numAtt;
			matSal.insert(pair<string,matGran>(codice,appo));
			ripristinaValori();
			delete[] indici;
		}
		else
		{	popol[i].setFitness(it->second.fitness);
			popol[i].setNregole(it->second.nRighe);
		}

	}
}


double soga::valutamatrice(int** matrice,double** training,int dimMat)
{
	frbs fis(numVar);
	double* pesi=new double[dimMat];

	double* usc=new double[numPatterTr];
	TuningPW=0;
	fis.createFuzzySet(numParts,0);
	fis.faiMatAttivazione(training,numPatterTr,numParts,matrice,dimMat,pesi,0);


	for (int i = 0; i < numPatterTr; i++)
		usc[i]=matrice[fis.calcolaClasse(pesi,i,(int**)matrice,dimMat)][numVar-1];

	int* numPuntixClasse=new int[numParts[numVar-1]];
	int** matConfusione=new int*[numParts[numVar-1]];

	for (int i=0;i<numParts[numVar-1];i++)
	{	numPuntixClasse[i]=0;
		matConfusione[i]=new int[numParts[numVar-1]];
		for (int j=0;j<numParts[numVar-1];j++)
			matConfusione[i][j]=0;
	}
	int noClass=0;

	for (int i=0;i<numPatterTr;i++)
	{	numPuntixClasse[(int)training[i][numVar-1]-1]++;
		if (usc[i]!=-1)
			matConfusione[(int)(training[i][numVar-1])-1][(int)usc[i]-1]++;
		else
		{		while(true)
				{	noClass=Randint(1,numParts[numVar-1]);
					if (noClass!=training[i][numVar-1])
						break;
				}
				matConfusione[(int)(training[i][numVar-1])-1][noClass-1]++;
		}
	}

	int numPunti=0;
	for (int i=0;i<numParts[numVar-1];i++)
		for (int j=0;j<numParts[numVar-1];j++)
			numPunti+=matConfusione[i][j];
	double uscita=0;
	for (int i=0;i<numParts[numVar-1];i++)
		uscita+=matConfusione[i][i];

	uscita/=(double)numPatterTr;
	TuningPW=1;

	for (int i=0;i<numParts[numVar-1];i++)
		delete[] matConfusione[i];
	delete[] matConfusione;
	delete[] numPuntixClasse;

	fis.deleteFuzzySet(numParts);
	return 1-uscita;
}


void soga::onePointCross(const crom& par1, const crom& par2, crom& suc1,crom& suc2)
{	int cut;
	cut = Randint(1,numGen);

	/*for (i=cut; i < numGen; i++)
	{	suc1.gen[i] = par2.gen[i];
		suc2.gen[i] = par1.gen[i];
	}
	int conta1=0;
	int conta2=0;
	for (int i=0;i<numGen;i++)
	{	if (suc1.gen[i]!=1)
			conta1++;
		if (suc2.gen[i]!=1)
			conta2++;
	}
	if (conta1<MINNFEAT || conta1>MAXNFEAT)
		suc1 = par1;
	if (conta2<MINNFEAT || conta2>MAXNFEAT)
			suc2 = par2;*/
}

void soga::mutation(crom& par) {

	int g, k;
	int* appo=par.getGen();
	g = Randint(0,numGen-1);
	int oldGran=appo[g];
	if (minTerms==numParts[g])
		return;
	if (appo[g] == minTerms && appo[g]<numParts[g])
	{	appo[g]++;
		par.setGen(appo);
		return;
	}
	if (appo[g] == numParts[g] && appo[g] >= minTerms)
		appo[g]--;
	else
	{	if (Rand() < 0.85) //increasing granularity
		{ if(appo[g]!=numParts[g])
			{	k = Randint(appo[g] + 1, numParts[g]);
				appo[g] = k;
			}
		}
		else //decreasing granularity
		{
			k = Randint(minTerms, appo[g] - 1);
			appo[g] = k;
		}
	}
	int conta=0;
	for (int i=0;i<numGen;i++)
	{	if (appo[g]!=1)
			conta++;
	}

	if (conta>=minTerms)
		par.setGen(appo);
}
