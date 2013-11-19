/***************************************************************************
 *   Copyright (C) 2008 by miki   *
 *   michela.antonelli@iet.unipi.it   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the numParts of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.             nocop1             *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "chromosome.h"
#include <fstream>
#include <vector>
#include <algorithm>
# define INF 1.0e14

float semilla, ProbVarPartion;

int maxRule, minRule, numPatterTr, numPatterTs, seed, maxMut, minTerms,	numVar, popsize, IndexObj, GrDig, Perc, DimChromR, limVat, MAXPARTO;
int TuningPW,TYPE;
unsigned int Seed;

bool vinc,CLASSIFICATION;

double** inOutTr;//training set
double** inOutTs;//test set
char file_test[300];
char MAXPARTI[100];

double* maxVal;
double* minVal;

double crossProb;//crossover probability
double mutProb;//mutation probability
double pmAdd;//addRule mutation probability
double pmRem;//remRule mutation probability
double pmMod; //modRule mutation probability

double probCrossTun;
double probCrossPart;
double probMutPart;
double mutProbTun;

int** matC45;
int dimmatC45;
int* numParts;

int Totaltrials,dimVettObj;

int MINPART = 2;
int SIZE_POP_PR; //size of the predictor population
int indiceO=5;  //5accuratezza  //3 auc //6 TPR Per l'ordinamento delle soluzioni in archivio


unsigned int** chromosome::matWM = 0;
unsigned int chromosome::dimMatWM = 0;
int* chromosome::indiciCL = 0;

//Funzione per ordinare il vettore degli indici
void chromosome::ordinaIndex() {
	int appo;
	bool* appoAn;
	bool trovato = false;
	for (int i = 0; i < maxRule - 1 && !trovato; i++)
	{	trovato = true;
		for (int j = maxRule - 1; j >= i + 1; j--)
			if (vettR[j].index < vettR[j - 1].index) {
				appo = vettR[j].index;
				vettR[j].index = vettR[j - 1].index;
				vettR[j - 1].index = appo;
				appoAn = vettR[j].vectAnt;
				vettR[j].vectAnt = vettR[j - 1].vectAnt;
				vettR[j - 1].vectAnt = appoAn;
				trovato = false;
			}
	}
	int j = maxRule - 1;
	for (; j >= 0 && vettR[j].index == dimMatWM; j--)	;
	numRul = j+1;
}

//Ordina la matrice in base agli indici della variabile di uscita (CLASSE)
void chromosome::ordinaxClasse(int** mat,int numR,int* index)
{	bool trovato=false;
	int* appo;
	int* appoIndex=new int[numParts[numVar-1]];
	for (int i=0;i<numParts[numVar-1];i++)
		appoIndex[i]=0;
	for (int i = 0; i < numR - 1 && !trovato; i++) //ordino per il campo classe
	{	trovato = true;
		for (int j = numR - 1; j >= i + 1; j--)
			if (mat[j][numVar-1] < mat[j - 1][numVar-1])
			{		appo = mat[j];
					mat[j]= mat[j-1];;
					mat[j-1]= appo;
					trovato = false;
			}
	}
	for (int i=0;i<numR;i++)
		appoIndex[mat[i][numVar-1]-1]++;
	index[0]=appoIndex[0];
	for (int i=1;i<numParts[numVar-1];i++)
		index[i]=index[i-1]+appoIndex[i];
	delete[] appoIndex;
}

///ordino la matrice di regole in base ai pesi
void chromosome::ordinaxPesi(int** mat,int numR,double* pesi)
{
	int* appoReg;
	bool trovato=false;
	double appopesi;

	for (int i = 0; i < numR - 1 && !trovato; i++) //ordino la matrice di regole in base ai pesi
		{	trovato = true;
			for (int j = numR - 1; j >= i + 1; j--)
				if (pesi[j] > pesi[j - 1])
				{	appoReg = mat[j];
					mat[j] = mat[j - 1];
					mat[j - 1] = appoReg;
					appopesi = pesi[j];
					pesi[j] = pesi[j - 1];
					pesi[j - 1] = appopesi;
					trovato = false;
				}
		}
}


//Elimina le regole uguali, in conflitto e quelle con peso inferiore a soglia
void chromosome::eliminaConliffitto(double* appopesi,double soglia, int** mat, int& dimMat)
{
	int* rule;
	double appoP;
	for (int k=0;k<dimMat;k++)  //elimino le regole con pesi minori di soglia
	{	if (appopesi[k]<=soglia)
		{	rule=mat[k];
			mat[k]=mat[dimMat-1];
			mat[dimMat-1]=rule;
			appoP=appopesi[k];
			appopesi[k]=appopesi[dimMat-1];
			appopesi[dimMat-1]=appoP;
			dimMat--;
			k--;
		}
	}

	int s;
	bool trovato; //,eliminatok;
	for (int k=0;k<dimMat-1;k++)  //elimino le regole uguali e quelle in conflitto
	{	s=k+1;
		while (s<dimMat)
		{	trovato=false;
			for (int j=0;j<numVar-1;j++)
				if (mat[k][j]!=mat[s][j])
				{	trovato=true;
					break;
				}
			if (!trovato) //antecedenti uguali
			{	if (mat[s][numVar-1]==mat[k][numVar-1]) //la regola k-esima e' in conflitto con la regola i-esima
				{	rule=mat[s];   //regole uguali elimino la regola s;
					mat[s]=mat[dimMat-1];
					mat[dimMat-1]=rule;
					appopesi[s]=appopesi[dimMat-1];
					dimMat--;
				}
				else
					s++;
			}
			else
				s++;
		}
	}
}

//Calcola l'errore della matrice matWM, lo scrive su file ed esce
void chromosome::calcolateWMError(frbs& fis)
{
	fstream f;
	f.open("WMErr", ios::out | ios::app);
	numRul=dimMatWM;
	matR=new unsigned int*[numRul];
	for (int i=0;i<numRul;i++)
	{	matR[i]=new unsigned[numVar];
		for (int j=0;j<numVar;j++)
			matR[i][j]=matWM[i][j];
	}

	double* usc=ECM(fis,inOutTr,numPatterTr,false);
	f<<numRul<<' '<<usc[0]<<' ';
	delete[] usc;
	usc=ECM(fis,inOutTs,numPatterTs,true);
	f<<usc[0]<<endl;
	f.close();
	exit(1);
}



void chromosome::generateMatCL(frbs& fis)
{
	//da scommentare per generazione wm
	/*double* appopesi=fis.calcolaPesi((int**)matWM,dimMatWM,realPart,inOutTr,numPatterTr,pwLT,tun, scaling);
	eliminaConliffitto(appopesi,0, (int**)matWM, (int&)dimMatWM);
	delete[] appopesi;

	*/
	if (dimMatWM<maxRule)
		maxRule=dimMatWM;


	indiciCL=new int[numParts[numVar-1]];
	for (int i=0;i<numParts[numVar-1];i++)
		indiciCL[i]=0;
	ordinaxClasse((int**)matWM,dimMatWM,indiciCL);
	#ifdef stampa
		stampaMatrice((int**)matWM,dimMatWM,numVar);
	#endif
}

//Genero il vettore delle regole selezionando dalla matrice matWM
void chromosome::generateIndiciCL(frbs& fis)
{
	int* reg;
	int numRegxClasse;
	int maxDim;
	int indmin;

	for (int i=0;i<maxRule;i++)
		vettR[i].index=dimMatWM;

	if (dimMatWM <= maxRule)  //Considero tutta la matrice matWM
	{	for (int i = 0; i < maxRule; i++)
		{	vettR[i].index = i;
			for (int j = 0; j < numVar-1; j++)
				if (matWM[i][j]!=0)
					vettR[i].vectAnt[j] = true;   //All'inizio tutti gli antecedenti della regola sono a 1
				else
					vettR[i].vectAnt[j] = false;
		}
		numRul=maxRule;
	}
	else
	{	//seleziono maxRule dalla matrice matWM
		numRegxClasse=maxRule/numParts[numVar-1];
		numRul=0;
		int numReg;
		for (int i=0;i<numParts[numVar-1];i++)
		{	if (i==0)
			{	numReg=indiciCL[i];
				indmin=0;
			}
			else
			{	numReg=indiciCL[i]-indiciCL[i-1];
				indmin=indiciCL[i-1];
			}
			if (numRegxClasse<=numReg)
				maxDim=numRegxClasse;
			else
				maxDim=numReg;
			reg=RandintDistinct(indmin, indiciCL[i] - 1, maxDim);
			for (int k = 0; k <maxDim ; k++)
			{	vettR[numRul].index = reg[k];
				for (int j = 0; j < numVar-1; j++)
					if (matWM[reg[k]][j]!=0)
						vettR[numRul].vectAnt[j] = true;   //All'inizio tutti gli antecedenti della regola sono a 1
					else
						vettR[numRul].vectAnt[j] = false;
				numRul++;
			}
			delete[] reg;
		}
	}

	ordinaIndex();
	faiMatdaVec();
	#ifdef stampa
		cout<<endl<<"Matrice Iniziale"<<endl;
		stampaMatrice((int**)matR,numRul,numVar);
	#endif

}


//genera la matrice matR dopo aver generato la matrice di WM (Per la regressione)
void chromosome::generateMatR(frbs& fis)
{	matR = new unsigned*[maxRule];
	for (int i=0;i<maxRule;i++)
		matR[i]=new unsigned [numVar];
	int maxDim;
	numRul=0;
	double peso;
	if (maxRule>=dimMatWM)   //la matrice di WM ci va tutta!!
	{	for (int j=0;j<dimMatWM;j++)
		{	for (int k=0;k<numVar;k++)
				matR[numRul][k]=matWM[j][k];
			//addDC(matR[numRul]);
			peso=fis.calcolaPeso((int*)matR[numRul],inOutTr,numPatterTr,realPart,pwLT);
			if (peso<=0)
				for (int k=0;k<numVar;k++)
					matR[numRul][k]=matWM[j][k];
			else
				pesi[numRul]=peso;
			numRul++;

		}
	}
	else
	{	int numRegxClasse=maxRule/numParts[numVar-1];
		for (int i=0;i<numParts[numVar-1];i++)
		{	if (indiciCL[i]>0)
			{	if (numRegxClasse<indiciCL[i])
					maxDim=numRegxClasse;
				else
					maxDim=indiciCL[i];

				int indice;
				int indmin;
				if (i==0)
					indmin=0;
				else
					indmin=indiciCL[i-1];

				for (int j=0;j<maxDim;j++)
				{	indice=Randint(0,indiciCL[i]-1);
					for (int k=0;k<numVar;k++)
						matR[numRul][k]=matWM[indice+indmin][k];
					//addDC(matR[numRul]);
					peso=fis.calcolaPeso((int*)matR[numRul],inOutTr,numPatterTr,realPart,pwLT);
					if (peso<=0)
						for (int k=0;k<numVar;k++)
							matR[numRul][k]=matWM[indice+indmin][k];
					else
						pesi[numRul]=peso;
					numRul++;
				}
			}
		}
	}
}


void chromosome::salvaDatiC45(frbs& fis)
{	fstream fout;
	fout.open("ErrTrIni.txt",ios::out | ios::app);

	double uscita=calcolaErroreC45(fis, inOutTr, numPatterTr,0, matWM,dimMatWM);
	fout<<uscita<<endl;
	fout.close();
	fout.clear();

	fout.open("ErrTsIni.txt",ios::out | ios::app);
	uscita=calcolaErroreC45(fis, inOutTs, numPatterTs,1, matWM,dimMatWM);
	fout<<uscita<<endl;
	fout.close();
	fout.clear();

	fout.open("C45Matrice.txt",ios::out | ios::app);
	fout<<"MatriceWM"<<endl;
	for(int i=0;i<dimMatWM;i++)
	{	for (int j=0;j<numVar;j++)
			fout<<matWM[i][j]<<' ';
		fout<<endl;
	}

	fout.close();
	fout.clear();

	fout.open("RegoleIni.txt",ios::out | ios::app);
	fout<<dimMatWM<<endl;

	fout.close();
	fout.clear();

	int comp=0;
	for (int i = 0; i < dimMatWM; i++)
		comp += countNoZero(matWM[i], numVar - 1);
	fout.open("CompIniziale.txt",ios::out | ios::app);
	fout<<comp<<endl;

	fout.close();
	fout.clear();

	exit(1);

}


//Creo la matrice di WM
void chromosome::inizializeWM(frbs& fis) {
	int numR;
	int** mat;
	if (CLASSIFICATION==0)
	{	fis.createFuzzySet(realPart, pwLT);
		mat = fis.generateWM(numR, realPart, inOutTr, numPatterTr); //, pwLT, tun,scaling);
		matWM = new unsigned*[numR];
		for (int i = 0; i < numR; i++)
			matWM[i] = (unsigned*) mat[i];
		dimMatWM = numR;
		fis.deleteFuzzySet(realPart);
	}
	else
	{	if (true)  // per la generazione della matrice con c45
		{	matWM=(unsigned int**)matC45;
			dimMatWM=dimmatC45;
		}
		else
		{	mat = fis.generateWMCL(numR,realPart,inOutTr, numPatterTr, pwLT); //Per la generazione con WM
			matWM = new unsigned*[numR];
			for (int i = 0; i < numR; i++)
				matWM[i] = (unsigned*) mat[i];
			dimMatWM = numR;
			delete[] mat;
		}
	}

	//salvaDatiC45(fis);
	if (maxRule > dimMatWM)
		maxRule = dimMatWM;
}

void chromosome::generateIndici(int maxIndex,int dimRules)
{
	int* reg = RandintDistinct(0, maxIndex - 1, dimRules);
	for (int i = 0; i < dimRules; i++)
	{	vettR[i].index = reg[i];
		for (int j = 0; j < numVar-1; j++)
			vettR[i].vectAnt[j] = true;   //All'inizio tutti gli antecedenti della regola sono a 1
	}
	//int appo;
	//bool trovato = false;

	ordinaIndex();
	faiMatdaVec();

	delete[] reg;
}

void chromosome::deleteDup() {
	bool* appoAn;
	unsigned int* appoReg;
	int i = 0;

	while (i<numRul) //Elimina le regole con  tutti gli antecedenti a zero
	{	if(countNoZero(matR[i], numVar - 1)==0)
		{	appoReg = matR[i];
			appoAn = vettR[i].vectAnt;
			for (int j = i; j < numRul - 1; j++)
			{
				matR[j] = matR[j + 1];
				vettR[j].index = vettR[j + 1].index;
				vettR[j].vectAnt = vettR[j + 1].vectAnt;
			}
			numRul--;
			matR[numRul] = appoReg;
			vettR[numRul].vectAnt = appoAn;
			vettR[numRul].index=dimMatWM;
		}
		else
			i++;
	}
	i=0;
	while (i < numRul) //elimina le regole duplicate
	{	int a = trovaUguali(i, matR[i], matR, numRul);
		if (a != -1) {
			appoReg = matR[i];
			appoAn = vettR[i].vectAnt;
			for (int j = i; j < numRul - 1; j++) {
				matR[j] = matR[j + 1];
				vettR[j].index = vettR[j + 1].index;
				vettR[j].vectAnt = vettR[j + 1].vectAnt;
			}
			numRul--;
			matR[numRul] = appoReg;
			vettR[numRul].vectAnt = appoAn;
			vettR[numRul].index=dimMatWM;

		} else
			i++;

	}
}

void chromosome::faiVecdaMat()
{	for (int i = 0; i < maxRule; i++)
		for (int j = 0; j < numVar-1; j++)
			if (matR[i][j]!=0)
				vettR[i].vectAnt[j] = true;
			else
				vettR[i].vectAnt[j] = false;
}


void chromosome::faiMatdaVec() {
	for (int i = 0; i < numRul; i++) {
		for (int j = 0; j < numVar-1; j++)
		{	if (!vettR[i].vectAnt[j])
				matR[i][j] = 0;
			else
				matR[i][j] = matWM[vettR[i].index][j];
		}
		matR[i][numVar-1] = matWM[vettR[i].index][numVar-1];
	}
}

/*************************************************************
 *************************************************************
 Allocates the necessary space for the chromosome variables
 except for tun which is allocated in inizializewithoutrule
 *************************************************************
 **************************************************************/
void chromosome::allocateSpace() {

	matR = new unsigned*[maxRule];
	pesi=0;

	if(CLASSIFICATION)
	{	pesi=new double[maxRule];
		matConfusione =new int*[numParts[numVar-1]];
		for (int i=0;i<numParts[numVar-1];i++)
			matConfusione[i]=new int[numParts[numVar-1]];
	}
	for (int i = 0; i < maxRule; i++)
		matR[i] = new unsigned[numVar];


	realPart = new int[numVar];
	int ingressi=numVar;

	if (TuningPW)
	{	if (CLASSIFICATION)
			ingressi=numVar-1;
		pwLT = new float*[ingressi];
		pwLT_lim[0] = new float*[ingressi];
		pwLT_lim[1] = new float*[ingressi];
		pwLT_lim[2] = new float*[ingressi];
		for (int i = 0; i < ingressi; i++)
		{	pwLT[i] = new float[numParts[i]-2];
			pwLT_lim[0][i] = new float[numParts[i]-2];
			pwLT_lim[1][i] = new float[numParts[i]-2];
			pwLT_lim[2][i] = new float[numParts[i]-2];
		}
	}
	sizeObjTot=dimVettObj;
	objTot = new double[sizeObjTot];
	objIndex = new int[sizeObjTot];
	objAlg=new double[sizeObjTot];
	sizeObjAlg=0;

	//seleziono gli indici degli obiettivi che regolano l'evoluzione
	int cifre=IndexObj;
	while (cifre!=0)
	{	objIndex[sizeObjAlg++]=cifre%10;
		cifre/=10;
	}

	if (TYPE == 1) //RULE SELECTION
	{	vettR = new gen[maxRule];
		for (int i = 0; i < maxRule; i++)
		{	vettR[i].index=dimMatWM;
			vettR[i].vectAnt = new bool[numVar-1];
			for (int j=0;j<numVar-1;j++)
				vettR[i].vectAnt[j]=true;
		}
	}

}

void chromosome::deleteChrom() {

	int ingressi=numVar;
	if (CLASSIFICATION)
	{	ingressi=numVar-1;
		delete[] pesi;
		for (int i=0;i<numParts[numVar-1];i++)
			delete[] matConfusione[i];
		delete[] matConfusione;
	}

	for (int i = 0; i < maxRule; i++)
		delete[] matR[i];
	delete[] matR;

	delete[] realPart;

	if (TuningPW) {
		for (int i = 0; i < ingressi; i++) {
			delete[] pwLT_lim[0][i];
			delete[] pwLT_lim[1][i];
			delete[] pwLT_lim[2][i];
			delete[] pwLT[i];
		}
		delete[] pwLT_lim[0];
		delete[] pwLT_lim[1];
		delete[] pwLT_lim[2];
		delete[] pwLT;
	}
	delete[] objTot;
	delete[] objAlg;
	delete[] objIndex;

	if (TYPE == 1) {
		for (int i = 0; i < maxRule; i++)
			delete[] vettR[i].vectAnt;
		delete[] vettR;

	}
}

//
 bool chromosome::controllaPresenzaClassi(int** mat,int Nreg)
 {
	 bool* okclasse = new bool[numParts[numVar-1]];
	 for (int i=0;i<numParts[numVar-1];i++)
		 okclasse[i]=false;
	 for (int i=0;i<Nreg;i++)
		 okclasse[mat[i][numVar-1]-1]=true;
	 int k;
	 for (k=0;k<numParts[numVar-1];k++)
	 if (!okclasse[k])
	 	break;
	 delete[] okclasse;
	 if (k!=numParts[numVar-1])
	 {
	 	return false;
	 }
	 return true;

 }


void chromosome::mutSelection(double probMut) {
	int ind, num;
	chromosome app(*this);
	int realNumR;

	bool trovato=true;
/*	if (numRul>minRule && Rand()<0.1)
	{	ind = Randint(0, numRul - 1);
		app.vettR[ind].index=dimMatWM;
		app.ordinaIndex();  //ordina gli indici e calcola il numRul
		app.faiMatdaVec();
		app.deleteDup();
		realNumR=app.numRul;
		if (realNumR >= minRule)
		{	if (!CLASSIFICATION)
				copyChrom(app);
			else
				if (controllaPresenzaClassi((int**)app.matR,app.numRul))
					copyChrom(app);
		}
	}*/

	if (dimMatWM > maxRule && Rand() < 0.1)  //modifico o aggiungo una regola
	{
		ind = Randint(0, maxRule - 1);	//elemento del vettore da modificare
		while(trovato)
		{	num = Randint(0, dimMatWM - 1);	//seleziono una regola che non e' gia' presente
			trovato=false;
			for (int i=0;i<numRul;i++)
			if (app.vettR[i].index==num)
			{	trovato=true;
				break;
			}
		}
		app.vettR[ind].index = num;
		/*if (app.vettR[ind].index==dimMatWM)	//sto aggiungendo
		{	*/
		for (int i=0;i<numVar-1;i++)
			if(matWM[num]!=0)
				app.vettR[ind].vectAnt[i]=true;
			else
				app.vettR[ind].vectAnt[i]=false;

		//}

		app.ordinaIndex();  //ordina gli indici e calcola il numRul

		app.faiMatdaVec();
		app.deleteDup();

		realNumR=app.numRul;
		if (realNumR >= minRule)
		{	if (!CLASSIFICATION)
				copyChrom(app);
			else
				if (controllaPresenzaClassi((int**)app.matR,app.numRul))
					copyChrom(app);
		}
	}

	int conta = 0;
	if (SIZE_POP_PR && Rand() < 0.7)
	{	ind = Randint(0, app.numRul - 1);
		for (int i = 0; i < numVar-1; i++)
			if (matWM[ind][i]!=0 && Rand() < (double) 2 / (numVar-1))
				app.vettR[ind].vectAnt[i] = !app.vettR[ind].vectAnt[i];

		for (int i = 0; i < numVar-1; i++)
			if (app.vettR[ind].vectAnt[i] != 0)
				conta++;
	}
	if (conta >= minTerms)
	{	app.faiMatdaVec();
		app.deleteDup();
		realNumR=app.numRul;
		if (realNumR >= minRule)
		{	if (!CLASSIFICATION)
				copyChrom(app);
			else
				if (controllaPresenzaClassi((int**)app.matR,app.numRul))
					copyChrom(app);
		}
	}
}
/*************************************************************
 *************************************************************
 Initializes a chromosome with random values for the DB parameters
 adding numInitRules random rules
 *************************************************************
 **************************************************************/
void chromosome::inizialize(frbs& fis) {

	inizializewithoutrule();
		//calcolaErroreC45(fis, inOutTr, numPatterTr,0, (unsigned**)matC45,dimmatC45);

	if (CLASSIFICATION)
	{	if (matWM==0)
		{	inizializeWM(fis);
			generateMatCL(fis);
		}
		if (TYPE==0)
			generateMatR(fis);
		else
			generateIndiciCL(fis);
	}
	else
	{	if (TYPE == 0)
		{	for (int i = 0; i < maxRule; i++)
				addRule(fis);
		}
		else
		{	if (matWM==0)
				inizializeWM(fis);
			generateIndici(dimMatWM,maxRule);
		}
	}
}

/*************************************************************
 *************************************************************
 Initializes a chromosome with random values without adding rule
 *************************************************************
 **************************************************************/
void chromosome::inizializewithoutrule() {
	rank = 1;
	crowd_dist = 0;
	//unsigned maxpart;
	int ingressi=numVar;
	numFeat=0;
	if (CLASSIFICATION)
		ingressi--;

	allocateSpace();

	float space;
	bool trovato;
	for (int i = 0; i < ingressi; i++)
	{	if (TuningPW)
		{	space = (double) 1 / (numParts[i]-2 + 1);// half of support of a fuzzy set in a uniform partition
			for (int j = 0; j < numParts[i]-2 ; j++) {
				if (vinc == true) //constrained version
				{
					pwLT_lim[0][i][j] = space * (j + 1) - space / 2;
					pwLT_lim[1][i][j] = space * (j + 1) + space / 2
							- (double) 1 / GrDig;
					//pwLT[i][j] = Randint(pwLT_lim[0][i][j] * GrDig,	pwLT_lim[1][i][j] * GrDig);
					pwLT[i][j]= space*(j+1)*GrDig;
					pwLT[i][j] /= (double) GrDig;
				} else {
					pwLT_lim[0][i][j] = 0.01;
					pwLT_lim[1][i][j] = 0.99;
					do {
						pwLT[i][j] = Randint(1, GrDig - 1);
						pwLT[i][j] /= (double) GrDig;
						trovato = false;
						for (int k = 0; k < j; k++)
							if (pwLT[i][j] == 0 || pwLT[i][j] == 1 || (j != 0
									&& abs(pwLT[i][j] - pwLT[i][k])
											< (double) (1 / GrDig))) //to ensure that two genes cannot be equal and
							{
								trovato = true;
								break;
							}
					} while (trovato);

				}

				pwLT_lim[2][i][j] = space * (j + 1);// centroids of the uniform partition in [0,1]
			}
			if (!vinc)
				ordina(pwLT[i], numParts[i]-2 );
		}
		for (int i=0;i<numVar;i++)
			realPart[i]=numParts[i];

	}

	numRul = realNumRule = 0;
	for (int i = 0; i < maxRule; i++)
		for (int j = 0; j < numVar;j++)
			matR[i][j] = 0;

	for (int i = 0; i < sizeObjTot; i++)
	{	objTot[i] = 56456465;
		objAlg[i]=56456465;
	}
}

chromosome::chromosome(const chromosome& source) {
	inizializewithoutrule();
	copyChrom(source);
}

/*************************************************************
 *************************************************************
 Copy the chromosome source in *this
 *************************************************************
 **************************************************************/
void chromosome::copyChrom(const chromosome& surc) {

	int ingressi=numVar;
	if (TYPE == 1) {
		for (int i = 0; i < maxRule; i++) {
			vettR[i].index = surc.vettR[i].index;
			for (int j = 0; j < numVar-1; j++)
				vettR[i].vectAnt[j] = surc.vettR[i].vectAnt[j];
		}
	}
	if (CLASSIFICATION)
	{	ingressi--;
		for (int i=0;i<maxRule;i++)
			pesi[i]=surc.pesi[i];
		for (int i=0;i<numParts[numVar-1];i++)
			for (int j=0;j<numParts[numVar-1];j++)
				matConfusione[i][j]=surc.matConfusione[i][j];
	}

	comp=surc.comp;
	interp = surc.interp;
	noCovered = surc.noCovered;
	for (int i = 0; i < maxRule; i++)
		for (int j = 0; j < numVar; j++)
			matR[i][j] = surc.matR[i][j];

	for (int i = 0; i < ingressi; i++) {
		realPart[i] = surc.realPart[i];
		if (TuningPW)
			for (int j = 0; j < numParts[i]-2; j++) {
				pwLT[i][j] = surc.pwLT[i][j];
				pwLT_lim[0][i][j] = surc.pwLT_lim[0][i][j];
				pwLT_lim[1][i][j] = surc.pwLT_lim[1][i][j];
				pwLT_lim[2][i][j] = surc.pwLT_lim[2][i][j];
			}
	}

	numRul = surc.numRul;
	realNumRule = surc.realNumRule;
	numFeat=surc.numFeat;
	sizeObjTot = surc.sizeObjTot;
	sizeObjAlg=surc.sizeObjAlg;

	for (int i = 0; i < sizeObjTot; i++)
	{	objTot[i] = surc.objTot[i];
		objIndex[i]=surc.objIndex[i];
	}

	for (int i=0;i<sizeObjAlg;i++)
		objAlg[i]=surc.objAlg[i];

	rank = surc.rank;
	crowd_dist = surc.crowd_dist;
	grid_loc = surc.grid_loc;
}

chromosome& chromosome::operator=(const chromosome& c1) {
	if (&c1 != this) {
		deleteChrom();
		inizializewithoutrule();
		copyChrom(c1);
	}
	return *this;
}

void chromosome::remRule()
{
	unsigned** matappo=new unsigned*[numRul];
	int realRule;
	int numero;
	unsigned* appo;
	for (int i=0;i<numRul;i++)
	{	matappo[i]=new unsigned[numVar];
		for (int j=0;j<numVar;j++)
			matappo[i][j]=matR[i][j];
	}

	realRule = deleteDup(matappo, numRul); //Eliminating the duplicates form the chromosome appo
	if (realRule>minRule)
	{	numero=Randint(0,numRul-1);
		appo=matR[numero];
		matR[numero]=matR[numRul-1];
		matR[numRul-1]=appo;
		numRul--;
	}

	for (int i=0;i<numRul;i++)
		delete[] matappo[i];
	delete[] matappo;
}



/*************************************************************
 *************************************************************
 Add a random rule
 *************************************************************
 **************************************************************/
void chromosome::addRule(frbs& fis) {

	int numNoZero;
	int flag = 1;
	unsigned** matappo;
	int num,num1;

	realNumRule = numRul;

	if (CLASSIFICATION)
	{	num=Randint(0,dimMatWM-1);
		if (numRul<maxRule)
		{	for (int i=0;i<numVar;i++)
				matR[numRul][i]= matWM[num][i];
			//addDC(matR[numRul]);
			numRul++;

		}
		else
		{	num1=Randint(0,maxRule-1);
			for (int i=0;i<numVar;i++)
				matR[num1][i]= matWM[num][i];
			//addDC(matR[num1]);
		}
		numRul = deleteDup(matR, numRul);
		return;
	}

	if (numRul != 0 && ProbVarPartion) //calcolo il concrete RB
	{	matappo = new unsigned*[numRul];
		for (int i = 0; i < numRul; i++) {
			matappo[i] = new unsigned[numVar];
			for (int j = 0; j < numVar; j++)
				matappo[i][j] = matR[i][j];
		}
		convertimat(matappo, fis);
		realNumRule = deleteDup(matappo, numRul); //Eliminating the duplicates form the chromosome appo

	} else
		matappo = matR;

	if ((ProbVarPartion && realNumRule == maxRule) || (ProbVarPartion == 0
			&& numRul == maxRule)) {
		if (matappo != matR) {
			for (int i = 0; i < numRul; i++)
				delete[] matappo[i];
			delete[] matappo;
		}
		return;
	}

	unsigned* candidate;

	while (flag == 1) {
		candidate = RandintRule(realPart, numVar); //extracting a  candidate rule
		numNoZero = countNoZero(candidate, numVar - 1);

		if (((numNoZero < minTerms) || (numNoZero > numVar-1)) || (numRul != 0	&& trovaVett(candidate, matappo, realNumRule))) {
			delete[] candidate;
			continue;
		};
		flag = 0;
	}
	if (matappo != matR) {
		for (int i = 0; i < numRul; i++)
			delete[] matappo[i];
		delete[] matappo;
	}

	if (numRul < maxRule) {
		matR[numRul] = candidate;
		numRul++;
	} else {
		num = Randint(0, numRul - 1);
		delete[] matR[num];
		matR[num] = candidate;
	}


}

/*************************************************************
 *************************************************************
 Modifies a term in the virtual RB
 *************************************************************
 **************************************************************/
void chromosome::modRule() {
	int m, f, numNoZero;
	chromosome app(*this);

	m = Randint(0, numRul - 1); //selecting the rule index

	if (!CLASSIFICATION)
		f = Randint(0, numVar - 1); //selecting the linguistic variable
	else
		f = Randint(0, numVar - 2);

	numNoZero = countNoZero(app.matR[m], numVar - 1);

	if (f == (numVar - 1)) //the consequent was selected
		app.matR[m][f] = Randint(1,realPart[numVar-1]);
	else //we are in an antecedent condition
	{
		if (numNoZero > minTerms && numNoZero < numVar-1)//no problems with the constrains
		{	if (Rand()<=0.7)
					app.matR[m][f]=0;
			else
				app.matR[m][f] = Randint(1, realPart[f]);
		}
		else
		{
			if (numNoZero == minTerms)//limit condition on the minimum number of conditions
			{	if (app.matR[m][f] == 0)
					app.matR[m][f] = Randint(0, realPart[f]);
				else
					app.matR[m][f] = Randint(1, realPart[f]);//the minimum number of rules will be saved
			}
			else //numNoZero==maxTerms: limit condition on the maximum number of conditions
			{	if (Rand()<=0.5)
					app.matR[m][f]=0;
				else
					app.matR[m][f] = Randint(1, realPart[f]);
			}
		}
	}// end extern else

	app.numRul = app.deleteDup(app.matR, app.numRul);

	if (app.numRul >= minRule) {
		copyChrom(app);

	}
}


/*************************************************************
 *************************************************************
 Mutation for the granularity part
 *************************************************************
 **************************************************************/
void chromosome::MutacionGR() {
	int g, k; //, maxpart;

	if (CLASSIFICATION)
		g = Randint(0, numVar - 2); //Selecting the linguistic variable
	else
		g = Randint(0, numVar - 1);

	int oldGran=realPart[g];
	if (realPart[g] == MINPART) {
		realPart[g]++;
		return;
	}
	if (realPart[g] == numParts[g]) {
		realPart[g]--;
	}
	else
	{
		if (Rand() < 0.85) //increasing granularity
		{ if(realPart[g]!=numParts[g])
			{	k = Randint(realPart[g] + 1, numParts[g]);
				realPart[g] = k;
			}
		}
		else //decreasing granularity
		{
			k = Randint(MINPART, realPart[g] - 1);
			realPart[g] = k;
		}
	}
	if (g==numVar-1 && realPart[g]==1)
	{	realPart[g]=oldGran;
		return;
	}
	//int conta=0;

}

/*************************************************************
 *************************************************************
 Mutation for the piecewise
 *************************************************************
 **************************************************************/
void chromosome::MutacionPW() {
	int g, h;
	int min1, max1; //range of variability of pwLT[i][j]
	int appo;
	double appo1;
	int ingressi=numVar;

	if (CLASSIFICATION)
		ingressi--;

	g = Randint(0, ingressi - 1); //selecting the variable to mutate
	if (numParts[g]==2)
	return;

	h = Randint(0, numParts[g]-3); //selecting the piece to mutate

	if (h == 0 || vinc) //the second fuzzy set
		min1 = pwLT_lim[0][g][h] * GrDig;
	else
		min1 = pwLT[g][h - 1] * GrDig;
	if (h == numParts[g]-2 - 1 || vinc) //the next-to-last fuzzy set
		max1 = pwLT_lim[1][g][h] * GrDig;
	else
		max1 = pwLT[g][h + 1] * GrDig;

	appo = Randint(min1, max1);
	appo1 = pwLT[g][h];
	pwLT[g][h] = (double) appo / GrDig;

	//to eliminate duplicate pwLT[][]
	for (int i = 0; i < numParts[g]-2; i++)
		if (i != h && abs(pwLT[g][h] - pwLT[g][i]) < (double) (1 / GrDig)) {
			pwLT[g][h] = appo1;
			break;
		}
}

/*************************************************************
 *************************************************************
 One-point crossover for the virtual RB
 *************************************************************
 **************************************************************/
void crossRB(chromosome& par1, chromosome& par2, chromosome& suc1,
		chromosome& suc2) {
	int i, minVal, cut;

	suc1.numRul = par2.numRul;
	suc2.numRul = par1.numRul;

	for (i = 0; i < maxRule; i++)
		for (int j = 0; j < numVar; j++) {
			suc1.matR[i][j] = par1.matR[i][j];
			suc2.matR[i][j] = par2.matR[i][j];
		}

	if (par1.numRul != 1) {
		if (par2.numRul != 1)
			minVal = minN(par1.numRul - 1, par2.numRul - 1);
		else
			//par2->numRul==1
			minVal = minN(par1.numRul - 1, par2.numRul);
	} else //par1->numRul==1
	{
		if (par2.numRul != 1)
			minVal = minN(par1.numRul, par2.numRul - 1);
		else
			//par2->numRul==1
			minVal = minN(par1.numRul, par2.numRul);
		// These nested if should ensure that the sons are different from parents if the number of genes of both parents are >1
	}
	cut = Randint(minRule, minVal);

	for (i = cut; i < maxRule; i++)
		for (int j = 0; j < numVar; j++) {
			suc1.matR[i][j] = par2.matR[i][j];
			suc2.matR[i][j] = par1.matR[i][j];
		}

	suc1.numRul = suc1.deleteDup(suc1.matR, suc1.numRul);
	suc2.numRul = suc2.deleteDup(suc2.matR, suc2.numRul);

	if (suc1.numRul < minRule)
		suc1.copyChrom(par1);


	if (suc2.numRul<minRule)
		suc2.copyChrom(par2);


	if ((suc1.numRul == 0) || (suc2.numRul == 0))
		cout << "One-point crossover error: son with 0 rule" << endl;

}

/*
void chromosome::checkChromosome(chromosome& parent,frbs& fis) {

	int realNumRule;
	unsigned** matAppo = new unsigned*[numRul];
	int oldnumRul = numRul;
	for (int i = 0; i < numRul; i++) {
		matAppo[i] = new unsigned[numVar];
		for (int j = 0; j < numVar; j++)
			matAppo[i][j] = matR[i][j];
	}
	convertimat(matAppo, fis);
	realNumRule = deleteDup(matAppo, numRul); //Eliminating the duplicates form the chromosome
	for (int i = 0; i < oldnumRul; i++)
		delete[] matAppo[i];
	delete[] matAppo;

	if (realNumRule < minRule) //ripristino la parte RB del padre
	{	matAppo=new unsigned*[parent.numRul];
		for (int i = 0; i < parent.numRul; i++)
		{	matAppo[i]=new unsigned[numVar];
			for (int j = 0; j < numVar; j++)
			{	matR[i][j] = parent.matR[i][j];
				matAppo[i][j] = parent.matR[i][j];
			}
		}
		numRul = parent.numRul;
		convertimat(matAppo,fis);
		realNumRule = deleteDup(matAppo, parent.numRul);
		if (realNumRule < minRule) //ripristino la parte FT del padre (caso peggiore figlio=padre)
			for (int i = 0; i < numVar; i++)
				realPart[i] = parent.realPart[i];
		for (int i = 0; i < parent.numRul; i++)
			delete[] matAppo[i];
		delete[] matAppo;
	}
}
*/


void crossWM(chromosome& par1, chromosome& par2, chromosome& suc1,chromosome& suc2)
{	int i, minVal, cut;

	suc1.numRul = par2.numRul;
	suc2.numRul = par1.numRul;

	for (i = 0; i < maxRule; i++)
	{	suc1.vettR[i].index=par1.vettR[i].index	;
		suc2.vettR[i].index=par2.vettR[i].index;
		for (int j = 0; j < numVar-1; j++)
		{	suc1.vettR[i].vectAnt[j]=par1.vettR[i].vectAnt[j];
		 	suc2.vettR[i].vectAnt[j]=par2.vettR[i].vectAnt[j];
		}
	}
	if (par1.numRul != 1 && par2.numRul != 1)
		minVal = minN(par1.numRul - 1, par2.numRul - 1);
	else
		return;

	cut = Randint(minRule, minVal);

	for (i=cut; i < maxRule; i++)
	{
		suc1.vettR[i].index = par2.vettR[i].index;
		suc2.vettR[i].index = par1.vettR[i].index;

		for (int j = 0; j < numVar-1; j++) {
			suc1.vettR[i].vectAnt[j] = par2.vettR[i].vectAnt[j];
			suc2.vettR[i].vectAnt[j] = par1.vettR[i].vectAnt[j];
		}
	}

	suc1.ordinaIndex();
	suc1.faiMatdaVec();
	suc1.deleteDup();

	suc2.ordinaIndex();
	suc2.faiMatdaVec();
	suc2.deleteDup();

	if (suc1.numRul<minRule || (CLASSIFICATION && !suc1.controllaPresenzaClassi((int**)suc1.matR,suc1.numRul)))
		suc1 = par1;

	if (suc2.numRul<minRule || (CLASSIFICATION && !suc2.controllaPresenzaClassi((int**)suc2.matR,suc2.numRul)))
		suc2 = par2;
}



/*************************************************************
 *************************************************************
 One-point crossover for the partition part
 *************************************************************
 **************************************************************/
void crossPart(chromosome& par1, chromosome& par2, chromosome& suc1,
		chromosome& suc2)//crosses the partition cromosome
{
	int num = Randint(1, numVar - 1);
	int i;
	//crossing the partition cromosome
	for (i = 0; i < num; i++) {
		suc1.realPart[i] = par1.realPart[i];
		suc2.realPart[i] = par2.realPart[i];
	}
	for (; i < numVar; i++) {
		suc1.realPart[i] = par2.realPart[i];
		suc2.realPart[i] = par1.realPart[i];
	}

	int conta1=0,conta2=0;
	for (int i=0;i<numVar-1;i++)
	{	if (suc1.realPart[i]>1)
			conta1++;
		if (suc2.realPart[i]>1)
				conta2++;
	}
}

void CruceBLXPW(chromosome& par1, chromosome& par2, chromosome& suc1,
		chromosome& suc2) {
	double temp, px, py, I, alfa;
	int ingressi=numVar;

	if (CLASSIFICATION)
		ingressi--;
	alfa = 0.5;
	//double space=(double)1/(numpieces+1);
	for (int j = 0; j < ingressi; j++)
	{
		for (int i = 0; i < numParts[j]-2; i++) {
			px = par1.pwLT[j][i];
			py = par2.pwLT[j][i];
			if (px > py) {
				temp = px;
				px = py;
				py = temp;
			}
			I = py - px;
			px = px - I * alfa;
			py = py + I * alfa;

			if (px < par1.pwLT_lim[0][j][i])
				px = par1.pwLT_lim[0][j][i];
			if (py >= par1.pwLT_lim[1][j][i])
				py = par1.pwLT_lim[1][j][i];

			suc1.pwLT[j][i] = px + Rand() * (py - px);
			suc2.pwLT[j][i] = px + Rand() * (py - px);
		}
		bool trovato = false;
		ordina(suc1.pwLT[j], numParts[j]-2);
		for (int i = 1; i < numParts[j]-2 && !trovato; i++)
			if (suc1.pwLT[j][i] - suc1.pwLT[j][i - 1] < (double) 1 / GrDig)
				trovato = true;
		if (trovato)
			for (int i = 0; i < numParts[j]-2; i++)
				suc1.pwLT[j][i] = par1.pwLT[j][i];

		ordina(suc2.pwLT[j], numParts[j]-2);
		trovato = false;
		for (int i = 1; i < numParts[j]-2 && !trovato; i++)
			if (suc2.pwLT[j][i] - suc2.pwLT[j][i - 1] < (double) 1 / GrDig)
				trovato = true;
		if (trovato)
			for (int i = 0; i < numParts[j]-2; i++)
				suc2.pwLT[j][i] = par2.pwLT[j][i];
	}

}



void chromosome::eliminaPesiNegativi()
{	unsigned* appomat;
	double appopesi;

	for (int i=0;i<numRul;i++)  //elimina le regole con pesi negativi
	{	if (pesi[i]<=0)
		{	appomat=matR[i];
			matR[i]=matR[numRul-1];
			matR[numRul-1]=appomat;
			appopesi=pesi[i];
			pesi[i]=pesi[numRul-1];
			pesi[numRul-1]=appopesi;
			for (int j=0;j<numVar-1;j++)
				if (matR[numRul-1][j]!=0)
					comp--;
			numRul--;
			realNumRule--;
			i--;
		}
	}
}

int chromosome::scegliClasse()
{	int* contaClasse=new int[numParts[numVar-1]];
	int max,indiceMax, indicePeso=0;
	double pesoMin=2;

	for (int i=0;i<numParts[numVar-1];i++)
		contaClasse[i]=0;
	for (int i=0;i<numRul;i++)
		contaClasse[matR[i][numVar-1]-1]++;
	max=contaClasse[0];
	indiceMax=0;
	for (int i=1;i<numParts[numVar-1];i++)
	if (max < contaClasse[i])
	{	max=contaClasse[i];
		indiceMax=i;
	}

	for (int i=0;i<numRul;i++)
		if (matR[i][numVar-1]==indiceMax || pesoMin>pesi[i])
		{	pesoMin=pesi[i];
			indicePeso=i;
		}
	delete[] contaClasse;
	return indicePeso;
}


double chromosome::calcolaErroreC45(frbs& fis, double** inOutCamp, int numCamp,bool test, unsigned int** matAppo,int nRegole)
{
	double uscita;

	fis.createFuzzySet(realPart, pwLT);
	fis.faiMatAttivazione(inOutCamp,numCamp,realPart,(int**)matAppo,nRegole,pesi,test);
	double* usc = new double[numCamp];

	for (int i = 0; i < numCamp; i++)
	{	int classe = fis.calcolaClasse(pesi,i,(int**)matAppo,nRegole);
		usc[i]=matAppo[classe][numVar-1];
		if (fis.getTutto0() == true)
		{	usc[i]=-1;

		}
	}


	int* numPuntixClasse=new int[numParts[numVar-1]];
	int ** matConfusione=new int*[numParts[numVar-1]];
	for (int i=0;i<numParts[numVar-1];i++)
	{	numPuntixClasse[i]=0;
		matConfusione[i]=new int[numParts[numVar-1]];
		for (int j=0;j<numParts[numVar-1];j++)
			matConfusione[i][j]=0;
	}
	int noClass=0;
	for (int i=0;i<numCamp;i++)
	{	numPuntixClasse[(int)inOutCamp[i][numVar-1]-1]++;
		if (usc[i]!=-1)
			matConfusione[(int)(inOutCamp[i][numVar-1])-1][(int)usc[i]-1]++;
		else
		{	while(true)
			{	noClass=Randint(1,numParts[numVar-1]);
				if (noClass!=inOutCamp[i][numVar-1])
					break;
			}
			matConfusione[(int)(inOutCamp[i][numVar-1])-1][noClass-1]++;
		}
	}

	int numPunti=0;
	for (int i=0;i<numParts[numVar-1];i++)
	{	for (int j=0;j<numParts[numVar-1];j++)
			numPunti+=matConfusione[i][j];
	}
	uscita=0;

	printMat<int>(matC45,dimmatC45,numVar,cout);
	cout<<endl;

	printMat<int>(matConfusione,numParts[numVar-1],numParts[numVar-1],cout);

	for (int i=0;i<numParts[numVar-1];i++)
		uscita+=matConfusione[i][i];
	uscita/=(double)numCamp;

	fis. deleteFuzzySet(realPart);

	return (uscita);
}

/*************************************************************
 *************************************************************
 Calculates the error on the training set and the complexity of
 the concrete RB
 - comp returns the complexity of the concrete RB
 - fis the fis is needed to convert the virtual RB into
 the concrete RB and to calculate the output
 - corr return the correlation if activated
 *************************************************************
 **************************************************************/
double* chromosome::ECM(frbs& fis, double** inOutCamp, int numCamp,bool test)
{
	double suma;
	unsigned** matAppo;
	bool nocop;
	double* usc;
	double output;
	double* uscCL=new double[10];
	int* numPuntixClasse;
	int** matConfusione;

	if (CLASSIFICATION);
		usc=new double[numCamp];

	comp = 0;
	realNumRule = numRul;

	if (ProbVarPartion != 0 ) //Partition learning
	{
		matAppo = new unsigned*[numRul];
		for (int i = 0; i < numRul; i++) {
			matAppo[i] = new unsigned[numVar];
			for (int j = 0; j < numVar; j++)
				matAppo[i][j] = matR[i][j];
		}
		convertimat(matAppo, fis);
		realNumRule = deleteDup(matAppo, numRul);
	}
	else
		matAppo = matR;

	numFeat=giveActive(matAppo, realNumRule);

	for (int i = 0; i < realNumRule; i++)
		comp += countNoZero(matAppo[i], numVar - 1); //calculating the complexity as sum of the antecedents != 0

	if (realNumRule < minRule)
		return 0;
	/*cout<<"Matrice da valuare"<<endl;
	for (int i=0;i<numRul;i++)
	{		for (int j=0;j<numVar;j++)
			cout<<matR[i][j]<<' ';
		cout<<endl;*/

	fis.createFuzzySet(realPart, pwLT);
	if (CLASSIFICATION)  //calcola i pesi per le regole
	{	fis.faiMatAttivazione(inOutCamp,numCamp,realPart,(int**)matAppo,realNumRule,pesi,test);
		int k;

		for (k=0;k<realNumRule;k++)
			if (pesi[k]> 0)
				break;

		if (k==realNumRule)
		{	uscCL[0]=1;
			uscCL[1] = 0; //TPR
			uscCL[2] = 1; //FPR
			uscCL[3] = 0;//AUC
			return uscCL;
		}
		if (TYPE==0)
		{	eliminaPesiNegativi();
			//if (numParts[numVar-1]!=2)
				//controllaClassi(fis);
		}
	}

	suma = 0;
	int nocop1=0;
	for (int i = 0; i < numCamp; i++)
	{	nocop = false;
		if (CLASSIFICATION)
		{	//usc[i]=fis.calcolaClasse(pesi,i,(int**)matAppo,realNumRule);
			usc[i]=matAppo[fis.calcolaClasse(pesi,i,(int**)matAppo,realNumRule)][numVar-1];
			if (fis.getTutto0() == true)
			{	usc[i]=-1;
				nocop1++;
			}

		}
		else
		{	output = FLC(inOutCamp[i], matAppo, realNumRule, fis, nocop); //calculating the output for each input pattern using the concrete RB
			suma += 0.5 * pow(inOutCamp[i][numVar - 1] - output, 2.0); //calculating the half of the MSE
		}
	}

	//fis.deleteFuzzySet(realPart);

	if (matAppo != matR) {
		for (int i = 0; i < numRul; i++)
			delete[] matAppo[i];
		delete[] matAppo;
	}

	if (CLASSIFICATION)
	{	numPuntixClasse=new int[numParts[numVar-1]];
		matConfusione=new int*[numParts[numVar-1]];

		for (int i=0;i<numParts[numVar-1];i++)
		{	numPuntixClasse[i]=0;
			matConfusione[i]=new int[numParts[numVar-1]];
			for (int j=0;j<numParts[numVar-1];j++)
				matConfusione[i][j]=0;
		}
		int noClass=0;

		for (int i=0;i<numCamp;i++)
		{	numPuntixClasse[(int)inOutCamp[i][numVar-1]-1]++;
			if (usc[i]!=-1)
				matConfusione[(int)(inOutCamp[i][numVar-1])-1][(int)usc[i]-1]++;
			else
			{	/*if (!test)
				{*/	while(true)
					{	noClass=Randint(1,numParts[numVar-1]);
						if (noClass!=inOutCamp[i][numVar-1])
							break;
					}
					matConfusione[(int)(inOutCamp[i][numVar-1])-1][noClass-1]++;
				//}

			}
		}

		int numPunti=0;
		for (int i=0;i<numParts[numVar-1];i++)
			for (int j=0;j<numParts[numVar-1];j++)
				numPunti+=matConfusione[i][j];

		uscCL[0]=0;
		for (int i=0;i<numParts[numVar-1];i++)
			uscCL[0]+=matConfusione[i][i];

		uscCL[0]/=(double)numCamp;
		uscCL[0]=1-uscCL[0];
		double TPR,FPR,prob,AUC=0;

		#ifdef stampa
			if (test)
			{	cout<<"Errore Test: "<<1-uscCL[0]<<endl<<"Matrice Confusione"<<endl;
				for (int i=0;i<MAXPARTO;i++)
				{	for (int j=0;j<MAXPARTO;j++)
						cout<<matConfusione[i][j]<<' ';
					cout<<endl;
				}

			}
		#endif

		if (numParts[numVar-1]==2)  //0 positivi 1 negativi
		{	uscCL[2] = (double)matConfusione[1][0]/numPuntixClasse[1]; //FPR
			uscCL[1] = (double)matConfusione[0][0]/numPuntixClasse[0]; //TPR
			uscCL[3] = (1+uscCL[1]-uscCL[2])/2;				//AUC
		}
		else
		{	prob=(double)1/numParts[numVar-1];
			TPR=0;
			FPR=0;
			AUC=0;
			double FP=0;
			double numF=0;

			for (int i=0;i<numParts[numVar-1];i++)
			{	TPR+=(double)matConfusione[i][i]/numPuntixClasse[i];
				FP=0;
				numF=0;
				for (int k=0;k<numParts[numVar-1];k++)
					if (i!=k)
					{	FP+=matConfusione[k][i];
						numF+=numPuntixClasse[k];
					}
				FPR+=(double)FP/numF;
			}
			uscCL[2]=FPR*prob;
			uscCL[1]=TPR*prob;
			AUC+=(1+uscCL[1]-uscCL[2])/2;
		}

		
		delete[] usc;
		delete[] numPuntixClasse;
		for (int i=0;i<numParts[numVar-1];i++)
			delete[] matConfusione[i];
		delete[] matConfusione;
	}
	else
		uscCL[0]= suma/(double)numCamp;

	fis.deleteFuzzySet(realPart);
	return uscCL;
	}

	/****************************************************************
	 *****************************************************************
	 Calculates the output of the FRBS given an input pattern and a RB
	 - var input pattern
	 - mat concrete RB
	 - dim number of rules of the concrete RB
	 - fis the fis is needed to calculate the output
	 - nocop returns true if the input pattern is not covered by the RB
	 *************************************************************
	 **************************************************************/
	double chromosome::FLC(double* var, unsigned** mat, int dim, frbs& fis,
			bool& nocop) {
		double (*tipoaggr)(double*, int, double);
		tipoaggr = minimo; // function for the AND operator

		return fis.calcOutput(var, *tipoaggr, realPart, (int**) mat, dim, nocop);
	}


	void chromosome::evaluateChrom(frbs& fis, double** inOutCamp, int numCamp) {

		double* uscite;
		//double* distTot = new double[numVar];

		uscite = ECM(fis, inOutCamp, numCamp,false);

		if (TuningPW)
			objTot[1] = 1-evaluateDistPW(); //objTot[1], fis);//questa va scommentata per il tuning normale
		else
			objTot[1]=0;

		objTot[0]=comp; //complessita
		objTot[2] =realNumRule; //numero di regole

		objTot[3]=objTot[6]=objTot[7]=0;
		objTot[4] = numFeat; //numero feature
		objTot[5] = uscite[0]; //accuratezza training

		if (CLASSIFICATION)
		{	objTot[3]=1-uscite[3]; //AUC
			objTot[6] = 1-uscite[1];  //TPR
			objTot[7] = uscite[2];  //FPR
		}

		for (int i=0;i<sizeObjAlg;i++)
			objAlg[i]=objTot[objIndex[i]-1];

		//delete[] distTot;
		delete[] uscite;

	}




	/*************************************************************
	 *************************************************************
	 Converts the matrix matReg of the virtual RB into matrix of
	 the concrete DB
	 - matReg the virtual RB.converted into the concrete RB
	 - fis a fis containing the structures of the DB.
	 *************************************************************
	 **************************************************************/
	void chromosome::convertimat(unsigned** matReg, frbs& fis) { //int ruota;
		for (int i = 0; i < numRul; i++)
			matReg[i] = convertiRule(matReg[i], fis);
	}


	/**************************************************************
	 ***************************************************************
	 Converts a rule of the virtual RB into a rule of the concrete DB
	 - rule the virtual rule
	 - fis a fis containing the structures of the DB.
	 *************************************************************
	 **************************************************************/
	unsigned* chromosome::convertiRule(unsigned* rule, frbs& fis) {
		unsigned * apporule = new unsigned[numVar];
		double* cent_partitions;
		double cent_max;
		double min_dist;
		int min_ind;

		for (int j = 0; j < numVar; j++)
		{	apporule[j] = 0;
			if (rule[j] != 0 && realPart[j]!=1)
			{
				cent_max = fis.getCentroids(realPart[j] - 2, rule[j] - 1);
				cent_partitions = fis.getVettCentroids(realPart[j]);
				min_ind = 0;
				min_dist = abs(cent_max - cent_partitions[0]);
				for (int s = 1; s < realPart[j]; s++)
					if (abs(cent_max - cent_partitions[s]) <= min_dist)
					{	min_dist = abs(cent_max - cent_partitions[s]);
						min_ind = s;
					}
				apporule[j] = min_ind + 1;
				delete[] cent_partitions;
			}
		}
		return apporule;
	}



	/**************************************************************
	 ***************************************************************
	 Deletes duplicate rules from the vector cc, reorders cc and
	 counts the number of rules
	 *************************************************************
	 **************************************************************/
	int chromosome::deleteDup(unsigned** mat, int numR) {
		int i = 0;
		int j;
		unsigned* appo;

		while (i < numR) //elimina le righe con tutti gli antecedenti a 0
		{
			j = 0;
			for (j = 0; j < numVar-1 && mat[i][j] == 0; j++)
				;
			if (j == numVar-1) {
				appo = mat[i];
				mat[i] = mat[numR - 1];
				mat[numR - 1] = appo;
				numR--;

			} else
				i++;
		}
		i = 0;
		while (i < numR)
			if (trovaUguali(i, mat[i], mat, numR) != -1) //controllo che la regola nn sia in RB
			{
				appo = mat[i];
				mat[i] = mat[numR - 1];
				mat[numR - 1] = appo;
				numR--;
			} else
				i++;

		return numR;
	}




	/*************************************************************
	 *************************************************************
	 Returns a vector that for each input linguistic label counts
	 the number of rules associated to
	 *************************************************************
	 **************************************************************/

	int chromosome::giveActive() {
		int q;
		unsigned * active = new unsigned[numVar - 1];
		for (int f = 0; f < numVar - 1; f++) {
			active[f] = 0;
			for (int m = 0; m < numRul; m++)
				if (matR[m][f] != 0)
					active[f]++;
		}
		q = countNoZero(active, numVar - 1);
		delete[] active;

		return q;
	}

	/*************************************************************
	 *************************************************************
	 Returns a vector that for each input linguistic label counts
	 the number of rules associated to
	 *************************************************************
	 **************************************************************/

	int chromosome::giveActive(unsigned** mat, int numR) {
		int q;
		unsigned * active = new unsigned[numVar - 1];

		for (int f = 0; f < numVar - 1; f++) {
			active[f] = 0;
			for (int m = 0; m < numR; m++)
				if (mat[m][f] != 0)
					active[f]++;
		}
		q = countNoZero(active, numVar - 1);
		delete[] active;

		return q;
	}



	/*************************************************************
	 *************************************************************
	 Returns the interpretability measure for the piecewise
	 transformation considering also the complexity of the RB
	 *************************************************************
	 **************************************************************/
	double chromosome::evaluateDistPW() {
		double somma = 0;
		double appo;
		int ingressi=numVar;
			if (CLASSIFICATION)
				ingressi--;

		for (int i = 0; i < ingressi; i++) {
			appo = 0;
			for (int j = 0; j < numParts[i]-2 ; j++)
				appo += fabs(pwLT[i][j] - pwLT_lim[2][i][j]);
			somma += appo;
		}

		return 2 * somma / ((double) ingressi * (numParts[0] - 2));

	}


	chromosome::~chromosome() {
		deleteChrom();

	}

	void chromosome::stampaPW() {
		int ingressi=numVar;
			if (CLASSIFICATION)
				ingressi--;
		for (int i = 0; i < ingressi; i++) {
			for (int j = 0; j < numParts[i]-2 ; j++)
				cout << pwLT[i][j] << ' ';
			cout << endl;
		}

	}

	//compare the chromosome taking into account the accuracy
	int compareChrom(const void * _a, const void * _b) {
		chromosome* elem1;
		chromosome* elem2;

		elem1 = (chromosome*) _a;
		elem2 = (chromosome*) _b;
		if ((*elem1).getObjTot(indiceO) < (*elem2).getObjTot(indiceO))
			return -1;
		else if ((*elem1).getObjTot(indiceO) > (*elem2).getObjTot(indiceO))
			return 1;
		else
			return 0;
	}

	void chromosome::printChrom(frbs& fis) {
		/*int * appo=new int[maxRule];
		 for (int i=0;i<maxRule;i++)
		 appo[i]=vettR[i].index;

		 order(appo, maxRule);*/
		for (int i = 0; i < numRul; i++)
		{	cout << vettR[i].index <<" ant:";
			for (int j=0;j<numVar-1;j++)
				cout<<vettR[i].vectAnt[j]<<' ';
			cout << endl;
		}

	/*	cout<<"Matrice"<<endl;
		for (int i = 0; i < numRul; i++) {
			for (int j = 0; j < maxTerms + 1; j++)
				cout << matR[i][j] << ' ';
			cout << endl;

		}*/

		for (int i=0;i<sizeObjAlg;i++)
			cout<<objAlg[i]<<' ';
		cout<<endl<<"pesi"<<endl;

		for (int i=0;i<numRul;i++)
			cout<<pesi[i]<<' ';
		cout<<endl;

	}








