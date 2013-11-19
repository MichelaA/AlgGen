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
#include "frbs.h"
#include "supportchrom.h"

#include <iostream>
#include <cmath>

# define MAXINT 2147483647
using namespace std;

//extern int numpieces;
extern int maxPW;
extern int* numParts;



frbs::frbs(int nVariable)
{
	fuzzyset=0;
	/*if (MAXPARTI>MAXPARTO)
		maxpart=MAXPARTI;
	else
		maxpart = MAXPARTO;*/

	dimMatAtt=maxRule;
	//cout<<maxpart<<endl;
	//partVar=new double**[maxpart-1];
	nVar =nVariable;

	matAttxPesi=new double*[numPatterTr+2];
	for (int i=0;i<numPatterTr+2;i++)
		matAttxPesi[i]=new double[dimMatAtt];


	dontCare=new double*[nVar];
	for (int i=0;i<nVar;i++)
	{	dontCare[i]=new double[2];  //0 per il minimo, 1 oper il max
		dontCare[i][0]=2*minVal[i]-maxVal[i];
		dontCare[i][1]=2*maxVal[i]-minVal[i];
	}
	/*for (int i=0;i<maxpart-1;i++)
		partVar[i]=	genTriUnifPartition(i+2,0,1);*/
	/*fuzzyset =new double**[nVar];
	for (int i=0;i<nVar;i++)
		fuzzyset[i]=0;*/
}


frbs::frbs(const frbs& fis)
{
	/*unsigned maxpart;
		if (MAXPARTI>MAXPARTO)
			maxpart=MAXPARTI;
		else
			maxpart = MAXPARTO;*/
	dimMatAtt=fis.dimMatAtt;
	//partVar=new double**[maxpart-1];
	nVar =fis.nVar;
	matAttxPesi=new double*[numPatterTr];
	for (int i=0;i<numPatterTr;i++)
	{	matAttxPesi[i]=new double[dimMatAtt];
		for (int j=0;j<maxRule;j++)
			matAttxPesi[i][j]=fis.matAttxPesi[i][j];
	}
	dontCare=new double*[nVar];
	for (int i=0;i<nVar;i++)
	{	dontCare[i]=new double[2];  //0 per il minimo, 1 oper il max
		dontCare[i][0]=2*minVal[i]-maxVal[i];
		dontCare[i][1]=2*maxVal[i]-minVal[i];
	}
	/*for (int i=0;i<maxpart-1;i++)
		partVar[i]=	genTriUnifPartition(i+2,0,1);*/
	/*fuzzyset =new double**[nVar];
	for (int i=0;i<nVar;i++)
		fuzzyset[i]=0;*/
}



frbs::~frbs()
{	for (int i=0;i<numPatterTr;i++)
		delete[] matAttxPesi[i];
	delete[] matAttxPesi;

	for (int i=0;i<nVar;i++)
		delete[] dontCare[i];
	delete[] dontCare;
	//delete[] fuzzyset;

}


void frbs::createFuzzySet(int* realPart,float** pwLT)
{
	fuzzyset =new double**[nVar];
	for (int i=0;i<nVar;i++)
		fuzzyset[i]= genTriUnifPartition(realPart[i],minVal[i],maxVal[i]);

	if (TuningPW && pwLT!=0)
		calcolafuzzyset(realPart,pwLT);
}

void frbs::deleteFuzzySet(int* realPart)
{
  for (int i=0;i<nVar;i++)
  {		for (int j=0;j<realPart[i];j++)
		  delete[] fuzzyset[i][j];
		delete[] fuzzyset[i];
  }
  delete[] fuzzyset;

}

int frbs::calcolaClasse(double* pesi, int indicePunto,int** matAppo,int nReg)
{
	if (nReg<=0)
		cout<<"ZERO"<<endl;
	double* att=new double[nReg];
	double maxAtt;
	int indiceRegola;
	tutto0=true;

	for (int i=0;i<nReg;i++)
	{	if (pesi[i]<=0)
			att[i]=-1;
		else
			att[i]=matAttxPesi[indicePunto][i]*pesi[i];
		if (att[i]>0)
			tutto0=false;
	}

	maxAtt=att[0];
	indiceRegola=0;  //cerco la regola con la massima attivazione
	for (int i=1;i<nReg;i++)
		if (att[i]>maxAtt)
		{	maxAtt=att[i];
			indiceRegola=i;
		}


	delete[] att;

	return indiceRegola;
}

void frbs::faiMatAttivazione(double** inOutCamp,int numCamp,  int* realPart, int** matReg,int nReg,double*& pesi,bool test)
{
    double* CFm = new double[nReg];
	double* NCFm = new double[nReg];


	if (nReg>dimMatAtt)
	{	for (int i=0;i<numPatterTr;i++)
			delete[] matAttxPesi[i];
		delete[] matAttxPesi;
		matAttxPesi=new double*[numPatterTr+2];
		for (int i=0;i<numPatterTr+2;i++)
			matAttxPesi[i]=new double[nReg];
		dimMatAtt=nReg;
		if (!test)
		{	if (pesi!=0)
				delete[] pesi;
			pesi = new double[nReg];
		}

	}

	for (int i=0;i<nReg;i++)
    {	CFm[i]=0;
		NCFm[i]=0;
    }

	for (int k=0;k<numCamp;k++)
	{	for(int i=0;i<nReg;i++) 							//Calculates in att[i] the activation of rule i
		{	matAttxPesi[k][i]=fuzzyInputCL(inOutCamp[k], realPart,matReg[i]);
			if (inOutCamp[k][nVar-1]==matReg[i][nVar-1])
				CFm[i]+=matAttxPesi[k][i];  //somma attivazione dei campioni di classe=classe della regola
			else
				NCFm[i]+=matAttxPesi[k][i];  //somma attivazionidei campioni di classe diversa dalla classe della regola
		}
	}

	if (!test)
	{	for (int i=0;i<nReg;i++)
			if ((CFm[i]+NCFm[i])==0)
				pesi[i]=0;
			else
				//pesi[i]=(CFm[i]-NCFm[i])/(CFm[i]+NCFm[i]);
				pesi[i]=(CFm[i])/(CFm[i]+NCFm[i]);
	}



	delete[] CFm;
	delete[] NCFm;
}

/**********************************************************************
***********************************************************************
Calculate the output of the FRBS
	- var input pattern
	- tipoaagr implication operator (set always as min)
	- realPart vector containing the granularity of each variable
	- matReg  matrix of rules
	- nReg	number of rules
	- nocop return parameter, specifies if the pattern is covered by the RB (nocop=false)
	- scaling,tun, pwLT parameters for the tuning phase
	- return  the output of the FRBS
**********************************************************************
**********************************************************************/

double frbs::calcOutput(double* var, double (*tipoaggr)(double*,int,double), int* realPart, int** matReg,int nReg,bool& nocop)
{
    double** in = fuzzyInput(var,realPart); 	// Input fuzzification
    double* appo=new double[nVar];
    double* att= new double[nReg];
    double sommapesi=0;
    double uscita=0;

	for(int i=0;i<nReg;i++) 							//Calculates in att[i] the activation of rule i
    {	for (int j=0;j<nVar-1;j++)
		{    appo[j]=1111;
		     if(matReg[i][j]!=0)
               appo[j]= in[j][matReg[i][j]-1];
		}
		att[i]=tipoaggr(appo,nVar-1,0);
     }

    int k;
	int kom;
    for (k=0;k<nReg && att[k]==0;k++);
    if (k==nReg)	//if there is no actived rule
	{	nocop=true;
		nearestRule(att,matReg,nReg,realPart,var,tipoaggr);
	}

	//uscita=normalDef(att,matReg,nReg,kom,fuzzyset[nVar-1]);
	uscita=WECOA(att,matReg,nReg,realPart[nVar-1],fuzzyset[nVar-1]);

   	delete[] att;
    delete[] appo;

    for (int i=0;i<nVar;i++)
	   delete[] in[i];

    delete[] in;

    return (uscita);
}


/*********************************************************************************
/********************************************************************************
	Fuzzifies the input pattern
 		- var input pattern
    	- fuzzyset contains the scaled and Tuned (if enabled) DB parameters
		- realPart vector containing the granularity of each variable
    	- return  for each variable and for each fuzzy set contains the membership degree to the fuzzy set
**********************************************************************************
*********************************************************************************/
double** frbs::fuzzyInput(double* varInput, int* vettPart)
{
	double** appart=new double*[nVar];
	for (int i=0;i<nVar;i++)
  	{   appart[i]=new double[vettPart[i]];
	  	for (int j=0;j<vettPart[i];j++)
      	{   appart[i][j]=0;
      		if (varInput[i]>=fuzzyset[i][j][1] && varInput[i]<=fuzzyset[i][j][2])
				appart[i][j]=1;
		  	else
				if (varInput[i]>=fuzzyset[i][j][0] && varInput[i]<=fuzzyset[i][j][1])
					appart[i][j]=(varInput[i]-fuzzyset[i][j][0])/(fuzzyset[i][j][1]-fuzzyset[i][j][0]);
				else
					if (varInput[i]>=fuzzyset[i][j][2] && varInput[i]<=fuzzyset[i][j][3])
						appart[i][j]=(varInput[i]-fuzzyset[i][j][3])/(fuzzyset[i][j][2]-fuzzyset[i][j][3]);
      	}
	}
	return appart;
}

//Return the matching degree of the rule on varInput
double frbs::fuzzyInputCL(double* varInput, int* vettPart,int* rule)
{
	double matDeg =1; //matching degree of the rule on the pattern varInput

	for (int i=0;i<nVar-1;i++)
		if (rule[i]!=0)
	  		matDeg*=evalMF(varInput[i],rule[i],i,vettPart[i]);
	  return matDeg;
}

double frbs::evalMF(double xf,int Af,int f,int F) //restituisce il grado di appartenenza del campione xf al fuzziset Af della var f (F num Fuzzy set)
{

	if (Af==1)
	{	if (xf>=fuzzyset[f][0][1] && xf<=fuzzyset[f][0][3])
			return (xf-fuzzyset[f][0][3])/(fuzzyset[f][0][1]-fuzzyset[f][0][3]);
		return 0;
	}

	if (Af==F)
	{	if (xf>=fuzzyset[f][F-1][0] && xf<=fuzzyset[f][F-1][1])
			return ((xf-fuzzyset[f][F-1][0])/(fuzzyset[f][F-1][1]-fuzzyset[f][F-1][0]));
		return 0;

	}

	if (xf>=fuzzyset[f][Af-1][0] && xf<=fuzzyset[f][Af-1][1])
		return ((xf-fuzzyset[f][Af-1][0])/(fuzzyset[f][Af-1][1]-fuzzyset[f][Af-1][0]));
	if (xf>fuzzyset[f][Af-1][1] && xf<=fuzzyset[f][Af-1][3])
		return ((xf-fuzzyset[f][Af-1][3])/(fuzzyset[f][Af-1][1]-fuzzyset[f][Af-1][3]));
	return 0;
}


/********************************************************************
*********************************************************************
Returns the centroid of the fuzzyset number index of granularity gran
	- gran granularity of the variable
    - index fuzzyset index (0-->gran-1)
It works only for triangular membership
/********************************************************************
*********************************************************************/
double frbs::getCentroids(int gran,int index)
{
	return partVar[gran][index][1];
}


/********************************************************************
*********************************************************************
Returns the vector of centroid of the MF of a linguistic variable
range [0--1] of granularity gran
	- gran granularity of the variable
/********************************************************************
*********************************************************************/
double* frbs::getVettCentroids(int gran)
{
	double* cent=new double[gran];
	for (int i=0;i<gran;i++)
		cent[i]=getCentroids(gran-2,i);
	return cent;
}




/**********************************************************************
***********************************************************************
Calculates the activation vector, choosing the two nearest rules of a
non-covered input pattern or move the input pattern in order to be covered
at least by two rules
	- att return the activation vector
    - matReg  matrix of rules
	- nReg	number of rules
	- realPart vector containing the granularity of each variable
 	- var input pattern
	- tipoaagr implication operator (set always as min)
	- fuzzyset contains the scaled and Tuned (if enabled) DB parameters
**********************************************************************
**********************************************************************/
void frbs::nearestRule(double* att, int** matReg, int nReg, int* realPart, double* var,double (*tipoaggr)(double*,int,double))
{

	// rango old function
	double* rangemin=new double[nVar-1];	// vector containing the minimum values of each linguistic variable
	double* rangemax=new double[nVar-1];	// vector containing the maximum values of each linguistic variable

	double* newvar = new double[nVar];
	newvar[nVar-1]=var[nVar-1];

	int first, second;
	for (int j=0;j<nVar-1;j++)
		if (matReg[0][j]==0)
		{	rangemin[j]=dontCare[j][0];
			rangemax[j]=dontCare[j][0];
		}
		else
		{	rangemin[j]=fuzzyset[j][matReg[0][j]-1][1];
			rangemax[j]=fuzzyset[j][matReg[0][j]-1][1];
		}

	for (int i=1;i<nReg;i++)
		for (int j=0;j<nVar-1;j++)
			if (matReg[i][j]!=0)
			{	if (fuzzyset[j][matReg[i][j]-1][1]<rangemin[j])
					rangemin[j]=fuzzyset[j][matReg[i][j]-1][1];
				if (fuzzyset[j][matReg[i][j]-1][1]>rangemax[j])
					rangemax[j]=fuzzyset[j][matReg[i][j]-1][1];
			}
			else
			{	if (dontCare[j][0]<rangemin[j])
					rangemin[j]=dontCare[j][0];
				if (dontCare[j][0]>rangemax[j])
					rangemax[j]=dontCare[j][0];
			}
	for (int i=0;i<nVar-1;i++)
		if ((rangemax[i]-rangemin[i])<=0)
			rangemax[i]=rangemin[i]+1;

	//end rango old function


	// calculating the two nearest rules to the input pattern
	// first nearest rule
	first=0;
	att[0]=euclidea(var,matReg[0],rangemin,rangemax,realPart);

	for (int i=1;i<nReg;i++)
	{	att[i]=euclidea(var,matReg[i],rangemin,rangemax,realPart);

		if (att[i]<att[first])
			first=i;

	}

	for (int i=0;i<nVar-1;i++)
	{
			newvar[i]=var[i];
			if (matReg[first][i]!=0)
			{	if (var[i]<=fuzzyset[i][matReg[first][i]-1][0])
					newvar[i]=fuzzyset[i][matReg[first][i]-1][0]+
					(fuzzyset[i][matReg[first][i]-1][3]
					-fuzzyset[i][matReg[first][i]-1][0])/2;
				else
					if (var[i]>=fuzzyset[i][matReg[first][i]-1][3])
						newvar[i]=fuzzyset[i][matReg[first][i]-1][3]-
						(fuzzyset[i][matReg[first][i]-1][3]
						-fuzzyset[i][matReg[first][i]-1][0])/2;
			}


	}

	//second nearest rule
	if (first==0)
		second=1;
	else
		second=0;

	att[second]=euclidea(newvar,matReg[second],rangemin,rangemax,realPart);
	for (int i=second+1;i<nReg;i++)
	{	att[i]=euclidea(newvar,matReg[i],rangemin,rangemax,realPart);
		if (att[i]<att[second]&& i!=first)
			second=i;

	}

	//if the consequents overlap we selected the nearest rule (the first)
	if ((fuzzyset[nVar-1][matReg[first][nVar-1]-1][1] >=
		fuzzyset[nVar-1][matReg[second][nVar-1]-1][1] &&
		fuzzyset[nVar-1][matReg[first][nVar-1]-1][0] <
		fuzzyset[nVar-1][matReg[second][nVar-1]-1][3])||
		(fuzzyset[nVar-1][matReg[first][nVar-1]-1][1] <=
		fuzzyset[nVar-1][matReg[second][nVar-1]-1][1] &&
		fuzzyset[nVar-1][matReg[first][nVar-1]-1][3] >
		fuzzyset[nVar-1][matReg[second][nVar-1]-1][0]))
	{	for (int i=0;i<nReg;i++)
			att[i]=0;
		att[first]=1;
		delete[] rangemin;
		delete[] rangemax;
		delete[] newvar;
		return;
	}


	double appo0, appo1, appo0_2;
	for (int i=0;i<nVar-1;i++)
	{
		if (matReg[first][i]!=0)
		{	appo0=fuzzyset[i][matReg[first][i]-1][0];
			appo1=fuzzyset[i][matReg[first][i]-1][3];
		}
		else
		{	appo0=dontCare[i][0];
			appo1=dontCare[i][1];
		}
		if (var[i]<=appo0)
			newvar[i]=appo0+(appo1-	appo0)/2;
		else
			if (var[i]>=appo1)
				newvar[i]=appo1-(appo1-appo0)/2;
			else
			{	if (matReg[second][i]!=0)
					appo0_2=fuzzyset[i][matReg[second][i]-1][1];
				else
					appo0_2=dontCare[i][0];
				if (var[i]<=appo0_2)
				{	newvar[i]=appo1-(appo1-appo0)*0.1;
					if (newvar[i]<var[i])
						newvar[i]=var[i];
				}
				else
				{	newvar[i]=appo0+(appo1-appo0)*0.1;
					if (newvar[i]>var[i])
						newvar[i]=var[i];
				}
			}

	}

	double** in = fuzzyInput(newvar,realPart);
    double* appo=new double[nVar];

	 //calculating the activation vector

	 for(int i=0;i<nReg;i++)
     { for (int j=0;j<nVar-1;j++)
		{    appo[j]=1111;
		     if(matReg[i][j]!=0)
               appo[j]= in[j][matReg[i][j]-1];
		}

		att[i]=tipoaggr(appo,nVar-1,0);

	 }
 	 for (int i=0;i<nVar;i++)
        delete[] in[i];

	delete[] in;
	delete[] appo;
	delete[] rangemin;
	delete[] rangemax;
	delete[] newvar;

}

/**********************************************************************
***********************************************************************
Calculates the output of the MFRBS, using the centre of gravity
	- att contains the activation vector
    - matReg  matrix of rules
	- nReg	number of rules
	- numPart granularity of the output variable
 	- fs fuzzyfied input pattern
**********************************************************************
**********************************************************************/
double frbs::WECOA(double* att,int** matReg,int nReg,int numPart,double** fs)
{
	double** cons=calccons(att, nReg, numPart, matReg,fs);
	double num=0,den=0;

	for (int i=0; i<nReg; i++)
    	if (att[i] != 0)
	  	{	num +=att[i]*(AreaTrapecioX(cons[i],att[i])/AreaTrapecio(cons[i],att[i]));
			den += att[i];
   	  	}
	if (den==0.0)
	{	num=0.0;
  		for (int i=0;i<nReg;i++)
   			num+=fs[matReg[i][nVar-1]-1][1];
		den = nReg;
 	}

 	for (int i=0;i<nReg;i++)
		delete[] cons[i];
	delete[] cons;

	if (den==0)
   		return (0);
	return (num/den);
}


/**********************************************************************
***********************************************************************
Calculate the 4 fuzzyset parameters for each rule
	- att contains the activation vector
	- nReg	number of rules
	- numPart granularity of the output variable
    - matReg  matrix of rules
 	- fs fuzzyfied input pattern
**********************************************************************
**********************************************************************/
double** frbs::calccons(double* att, int nReg, int numPart, int** matReg,double** fs)
{
	double** cons=new double* [nReg];

	for (int i=0;i<nReg;i++)
	{	cons[i]=new double[4];
		if (att[i] != 0)
		{	if (att[i] == 1.0)
			{	cons[i][0] = fs[matReg[i][nVar-1]-1][0];
				cons[i][1] = fs[matReg[i][nVar-1]-1][1];
				cons[i][2] = fs[matReg[i][nVar-1]-1][2];
				cons[i][3] = fs[matReg[i][nVar-1]-1][3];
			}
			else
			{
				cons[i][0] = fs[matReg[i][nVar-1]-1][0];
				cons[i][1] = cons[i][0] + (fs[matReg[i][nVar-1]-1][1]-cons[i][0])*att[i];
				cons[i][2] = fs[matReg[i][nVar-1]-1][3]+(fs[matReg[i][nVar-1]-1][2]-
							fs[matReg[i][nVar-1]-1][3])*att[i];
				cons[i][3] =fs[matReg[i][nVar-1]-1][3];
			}
		}
	}
	return cons;
}


/*********************************************************************************
/*********************************************************************************
 Generates a uniform fuzzy partition between maxVal and minVal
	-	indice the fuzzy set number of the partition (2-->MAXPART)
    - 	valMin the minimum value of the variable
	-	valMax the maximum value of the variable
    return  a  [indicex4] matrix each row contains the parameter the correspondent fuzzy set
**********************************************************************************
*********************************************************************************/
double** frbs::genTriUnifPartition(int indice, double valMin,double valMax)
{
	double interv=(valMax-valMin)/(indice-1);
	double** appo=new double*[indice];
	for (int i=0;i<indice;i++)
		appo[i]=new double [4];
	for(int i=0;i<indice-1;i++)
    	if (i==0)
	 	{	appo[0][0]=valMin-interv;
			appo[0][1]=appo[0][2]=valMin;
			appo[0][3]=valMin+interv;
	 	}
     	else
     	{	appo[i][0]=appo[i-1][1];
		 	appo[i][1]=appo[i][2]=appo[i-1][3];
         	appo[i][3]=appo[i][1]+interv;
    	}
	appo[indice-1][0]=appo[indice-2][1];
	appo[indice-1][1]=appo[indice-1][2]=valMax;
	appo[indice-1][3]=valMax+interv;
	return appo;
}



/*********************************************************************************
/*********************************************************************************
Performs the scaling and tuning of the database parameters
	- fuzzyset return the scaled and Tuned (if enabled) DB parameters
	- realPart vector containing the granularity of each variable
	- pwLT  parameters for the piecewise linear transormation
	- scaling  parameters for the nonlinear transformation of the fuzzy sets (if activeted)
	- tun  parameters for the 2tuples tuning
**********************************************************************************
*********************************************************************************/
void frbs::calcolafuzzyset(int* realPart,float** pwLT)
{
	if (TuningPW)  //piecewise linear transformation
	{	double space;
		int indice;

		int ingressi=nVar;
		if (CLASSIFICATION)
			ingressi--;
		double valore;
		for (int i=0;i<ingressi;i++)
		{	space=(double)1/(numParts[i]-1);
			if (realPart[i]!=2)
			{	retta* rette = calcolarette(pwLT[i],numParts[i]);
				for (int j=0;j<realPart[i];j++)
				{	for (int k=0;k<4;k++)
					{	valore=(fuzzyset[i][j][k]-minVal[i])/(maxVal[i]-minVal[i]);
						if ((j==0 && k==3)|| (j==1 && k!=0) || (j==realPart[i]-1 && k==0) || (j==realPart[i]-2 && k!=4) || (j!=0 && j!=1 && j!=realPart[i]-2 && j!=realPart[i]-1))
						{	indice=0;
							for (int s=0;s<numParts[i]-1;s++)
								if (valore>s*space && valore<=(s+1)*space)
								{	indice=s;
									break;
								}
							 // indice del fuzzyset normalizzato
							if (rette[indice].p==MAXINT)
								fuzzyset[i][j][k]=rette[indice].q;
							else
								fuzzyset[i][j][k]=(valore-rette[indice].q)/rette[indice].p;//linear transofrmation using the selected segment
							fuzzyset[i][j][k]=minVal[i]+(maxVal[i]-minVal[i])*fuzzyset[i][j][k];	//scaling the fuzzy sets
						}
					}
				}
				delete[] rette;
			}
		}
	}
}

/**********************************************************************
**********************************************************************
Calculates the coefficients p and q for each linear segment
	- pwLTi parameters for the piecewise linear transformation of the i-th variable
	return the p and q coefficients for each linear segment
**********************************************************************
**********************************************************************/

retta* frbs::calcolarette(float* pwLTi,int npart)
{
	double space=(float)(1)/(npart-1);

	int numpieces=npart-2;
	retta* rette= new retta[numpieces+1];


	rette[0].p=space/pwLTi[0];
	rette[0].q=0;

	rette[numpieces].p=(1-numpieces*space)/(1-pwLTi[numpieces-1]);
	rette[numpieces].q=1-rette[numpieces].p;

	for (int i=1;i<numpieces;i++)
	{	if (pwLTi[i]==pwLTi[i-1])
		{	rette[i].p=MAXINT;
			rette[i].q = space*i;
		}
		else
		{	rette[i].p = space/(pwLTi[i]-pwLTi[i-1]);
			rette[i].q = space*i -rette[i].p*pwLTi[i-1];
		}
	}
	return rette;
}

double frbs::normalDef(double* att,int** matReg,int nReg,int part,double** fs)
{
	double uscita=0, sommapesi=0;

	for (int i=0;i<nReg;i++)
	{   sommapesi+=att[i];
		uscita+=att[i]*	fs[matReg[i][nVar-1]-1][1];
	}

	if (sommapesi==0.0)
	{	uscita=0.0;
  		for (int i=0;i<nReg;i++)
   			uscita+=fs[matReg[i][nVar-1]-1][1];
  		sommapesi=nReg;
 	}
	if (sommapesi==0)
		uscita=0;
	else
		uscita/=sommapesi;

	return uscita;
}


double AreaTrapecioX(double* fuzzyset, double att)
{
	double izq, centro, der;
	double x0=fuzzyset[0];
	double x1=fuzzyset[1];
	double x2=fuzzyset[2];
	double x3=fuzzyset[3];
   if (fuzzyset[1] != fuzzyset[0])
      izq = (2*x1*x1*x1 - 3*x0*x1*x1 + x0*x0*x0) / (6*(x1-x0));
   else
      izq = 0;

   centro = (x2*x2-x1*x1)/2;

   if (fuzzyset[3] != fuzzyset[2])
      der = (2*x2*x2*x2 - 3*x3*x2*x2 + x3*x3*x3) / (6*(x3-x2));
   else
      der = 0;

   return (att * (izq+centro+der));
}

double AreaTrapecio(double* fuzzyset, double att)
{
	double izq, centro, der;
	double x0=fuzzyset[0];
	double x1=fuzzyset[1];
	double x2=fuzzyset[2];
	double x3=fuzzyset[3];

	if (fuzzyset[1]!=fuzzyset[0])
      izq = (x1*x1 - 2*x0*x1 + x0*x0) / (2*(x1-x0));
   else
      izq = 0;

   centro = x2-x1;

   if (fuzzyset[3] != fuzzyset[2])
      der = (x3*x3- 2*x2*x3 + x2*x2) / (2*(x3-x2));
   else
      der = 0;

   return (att * (izq+centro+der));

}


double frbs::euclidea(double* var,int reg[],double rangemin[], double rangemax[],int* realPart)
{
	double dist=0;
	for (int i=0;i<nVar-1;i++)
	{
		if (reg[i]!=0)
			dist+=pow((var[i]-fuzzyset[i][reg[i]-1][1])/(rangemax[i]-rangemin[i]),2);
		else
			dist+=pow((var[i]-dontCare[i][0])/(rangemax[i]-rangemin[i]),2);

	}
	return sqrt(dist);

}


double minimo(double* vett,int dim, double gamma)
{
       double m;
       for (int i=0;i<dim;i++)
           if (i==0)
              m=vett[0];
           else
               if(m>vett[i])
                   m=vett[i];
       return m;
}


int** frbs::generateWM(int& NRule, int* vettPart, double** pattern, int nPattern) //,float** pwLT) //,float* tun, scalingFunction& scaling)
{
	int** mat=new int*[nPattern];
	double*** fuzzySet=new double**[numVar];
	int* campioni=new int[nPattern];

	int indice;
	double** appo;
	int* rule;
	int numP, indMax;
	double valMax;
	bool trovato;
	double* attiv=new double[nPattern];

	for (int i=0;i<nPattern;i++)
	{	campioni[i]=0;
		mat[i]=0;
		attiv[i]=0;

	}


	double* dist=0;
//createFuzzySet(vettPart,pwLT);

	NRule=0;
	numP=0;
	double att=1;

	while(numP<nPattern)
	{
		rule=new int[numVar];
		appo=fuzzyInput(pattern[numP],vettPart);
		att=1;
		for (int i=0;i<numVar;i++)   //genero la regola
		{	indMax=0;
			valMax=appo[i][0];
			for(int j=1;j<vettPart[i];j++)
				if (valMax<appo[i][j])
				{	valMax=appo[i][j];
					indMax=j;
				}
			rule[i]=indMax+1;
			att*=appo[i][indMax];
		}




		trovato=true;
		int k;
		double p1,p2;
		for(k=0;k<NRule && trovato; k++)
		{	trovato=false;
			for (int j=0;j<numVar-1;j++)
				if (rule[j]!=mat[k][j])
				{	trovato=true;
					break;
				}
			if (!trovato)
			{	if (rule[numVar-1]!=mat[k][numVar-1]) //la regola k-esima e' in conflitto con rule
				{	if (attiv[k]<att) //cambio la regole in mat con rule perche' ha grado piu alto
					{   mat[k]=rule;
						attiv[k]=att;
					}
					else
					   delete[] rule;
				}
				else
				  delete[] rule;
			}
		}
		if (trovato) //non ci sono regole uguali o in conflitto
		{	mat[NRule]=rule;
		    attiv[NRule]=att;
			NRule++;
		}

		for (int i=0;i<numVar;i++)
			delete[] appo[i];
		delete[] appo;
		numP++;
	}
    //deleteFuzzySet(vettPart);

	delete[] campioni;

	return mat;
}



int** frbs::generateWMCL(int& NumR, int* vettPart, double** pattern, int nPattern,float** pwLT)
{
	int** mat=new int*[nPattern];
	//double*** fuzzySet=new double**[numVar];
	double** appo;
	int* rule;
	int indMax;
	double valMax, CF,NCF;

	//double* dist=0;
	createFuzzySet(vettPart,pwLT);

	for (int k=0;k<nPattern;k++) //creo tutta la matrice
	{	//pesi[k]=2;
		mat[k]=new int[numVar];
		appo=fuzzyInput(pattern[k],vettPart);

		for (int i=0;i<numVar-1;i++)   //genero la regola
		{	indMax=0;
			valMax=appo[i][0];
			for(int j=1;j<vettPart[i];j++)
				if (valMax<appo[i][j])
				{	valMax=appo[i][j];
					indMax=j;
				}
			mat[k][i]=indMax+1;
		}
		mat[k][numVar-1]=pattern[k][numVar-1];

		for (int i=0;i<numVar;i++)
			delete[] appo[i];
		delete[] appo;
	}
	NumR=nPattern;


	/*for (int i=0;i<nPattern;i++)
	{	for (int j=0;j<numVar;j++)
			cout<<mat[i][j]<<' ';
		cout<<endl;
	}*/



	return mat;
}

double* frbs::calcolaPesi(int** mat,int NumR,int* vettPart, double** pattern, int nPattern,float** pwLT) //calcolo i pesi della matrice mat
{
	createFuzzySet(vettPart,pwLT);
	double p;
	double* pesi=new double[NumR];
	double CF,NCF;
	for (int k=0;k<NumR;k++)
	{	CF=0;
		NCF=0;
		for (int j=0;j<nPattern;j++)  //per ogni punto
		{	p=fuzzyInputCL(pattern[j], vettPart,mat[k]);
			if (pattern[j][numVar-1]==mat[k][numVar-1])
				CF+=p;  //somma attivazione dei campioni di classe=classe della regola
			else
				NCF+=p;  //somma attivazionidei campioni di classe diversa dalla classe della regola
		}
		//pesi[k]=(CF-NCF)/(CF+NCF);
		pesi[k]=(CF)/(CF+NCF);
	}
	return pesi;
}


double frbs::calcolaPeso(int* rule, double** pattern, int nPattern,int* vettPart,float** pwLT) //calcolo i pesi della matrice mat
{
	createFuzzySet(vettPart,pwLT);
	double p;

	double CF,NCF;
	CF=0;
	NCF=0;
	for (int j=0;j<nPattern;j++)  //per ogni punto
	{	p=fuzzyInputCL(pattern[j], vettPart,rule);
		if (pattern[j][numVar-1]==rule[numVar-1])
			CF+=p;  //somma attivazione dei campioni di classe=classe della regola
		else
			NCF+=p;  //somma attivazionidei campioni di classe diversa dalla classe della regola
	}
	if (CF+NCF!=0)
		//return (CF-NCF)/(CF+NCF);
		return (CF)/(CF+NCF);
	return 0;
}


