/*
 * supportchrom.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef SUPPORTCHROM_H_
#define SUPPORTCHROM_H_

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "frbs.h"
#include <cstring>


using namespace std;

#define MASK 2147483647
#define PRIME 65539
#define SCALE 0.4656612875e-9
/***************************************************************
****************************************************************
To read the algorithm parameters from the configuration file
***************************************************************
****************************************************************/

#define FORMATO_ENT "File with the training set = %s \
File with test set = %s \
Number of Evaluation = %d \
Population length = %d \
Maximum number of rules = %d \
Minimum number of rules = %d \
Minimum number of antecedents = %d \
Prob RBcross = %lf \
Prob RBmutation = %lf \
pmAdd = %lf \
pmRem = %lf \
pmMod = %lf \
maxMut = %d \
probCrossTun = %lf \
probCrossPart = %lf \
probMutTun = %lf \
probMutPart = %lf \
Seed = %f \
ProbVarPartion = %f \
TuningPW = %d \
IndexObj = %d \
dimVettObj = %d \
vinc = %d \
Granularity of digits = %d \
Reduction Perc = %d \
Evaluation Interval = %d \
DimChromR = %d \
MAXPARTI = %s \
MAXPARTO = %d \
TYPE = %d \
SIZE_POP_PR = %d \
CLASSIFICATION = %d"

#define FORMATO_SAL "Input parameters read :\n\n\
File with the training set = %s \n\
File with test seta = %s \n\
Number of Evaluation = %d \n\
Population length = %d \n\
Maximum number of rules  = %d \n\
Minimum number of rules  = %d \n\
Minimum number of antecedents  = %d \n\
Prob RBcross = %f \n\
Prob RBmutation = %f \n\
pmAdd = %f \n\
pmRem = %f \n\
pmMod = %f \n\
maxMut = %d \n\
probCrossTun = %f \n\
probCrossPart = %f \n\
probMutTun = %f \n\
probMutPart = %f \n\
Seed = %f \n\
ProbVarPartion = %f \n\
TuningPW = %d \n\
IndexObj = %d \n\
dimVettObj= %d \n\
vinc = %d \n\
Granularity of digits = %d \n\
Reduction Perc = %d \n\
Evaluation Interval = %d \n\
DimChromR = %d \n\
MAXPARTI = %s \n\
MAXPARTO = %d \n\
TYPE = %d \n\
SIZE_POP_PR = %d \n\
CLASSIFICATION = %d "

#define VAR_ENT  file_tra, file_test, \
                 &Totaltrials, &popsize, &maxRule, &minRule, &minTerms, &crossProb, &mutProb, &pmAdd, &pmRem, \
                 &pmMod, &maxMut, &probCrossTun,&probCrossPart, &mutProbTun, &probMutPart, &semilla, \
				 &ProbVarPartion, &TuningPW,&IndexObj,&dimVettObj, &vinc, &GrDig, &Perc, &limVat, &DimChromR, MAXPARTI, &MAXPARTO, &TYPE, \
				 &SIZE_POP_PR, &CLASSIFICATION

#define VAR_SAL  file_tra,file_test, \
                 Totaltrials, popsize, maxRule, minRule,  minTerms,  crossProb, mutProb, pmAdd, pmRem,\
                 pmMod, maxMut, probCrossTun,probCrossPart, mutProbTun,probMutPart,semilla, \
                 ProbVarPartion, TuningPW, IndexObj,dimVettObj,vinc, GrDig,Perc, limVat, DimChromR, MAXPARTI, MAXPARTO, TYPE, \
                 SIZE_POP_PR, CLASSIFICATION

extern char file_test[300];
extern char file_tra[300];
extern char MAXPARTI[100];
extern int** matC45;
extern int maxRule, minRule, numPatterTr, numPatterTs,minTerms, numVar, maxMut,IndexObj, GrDig, Perc,limVat,DimChromR;
extern float semilla, ProbVarPartion;
extern double**inOutTr;
extern double** inOutTs;
extern double* maxVal, * minVal;
extern double crossProb, mutProb ,pmAdd, pmRem, pmMod, probCrossTun, probMutPart, probCrossPart, mutProbTun;
extern int Totaltrials, popsize,MAXPARTO, TuningPW;
extern unsigned int Seed;
extern bool vinc, CLASSIFICATION;
extern int dimmatC45;

extern int* numParts;

extern int TYPE,dimVettObj;
extern int SIZE_POP_PR;

typedef struct lists
{
    int index;
    struct lists *parent;
    struct lists *child;
}
list;
void insert (list *, int);
list* del (list *);

void order(int*, int);
int compareDec(const void*, const void*);///< Used in quicksort to arrange element in ascending order
int compareCre(const void*, const void*);
//int compareDecSup(const void*, const void*);

bool controllaPresenzaClassi(int** mat,int Nreg);

int compare_min(double*, double*, int);///< find the minimum value in a double vector
int equal(double *, double *, int ); ///< Cecks to n-dimensional vectors of objectives to see if they are identical returns 1 if they are, 0 otherwise

void calcolaWMErr();
double calcolaIns(vector<int>,double**,vector<double>&,vector<double>&);

/**
 	*  \param min minimum value of each generated element
    *  \param max maximum value of each generatedelement
	*  \param num number of elements to be generated
    *  \return returns num distinct random integers in [min,max]
*/
int* RandintDistinct(int min,int max,int num); ///< Returns num distinct random integer between [min,max], allocating the space for the return vector
unsigned* RandintRule(int* max,int num);
int* RandintDistinctOrd(int min,int max,int num); ///< Returns num distinct random integer between [min,max], allocating the space for the return vector which is ordered

void calcolaWMrandom(vector<int>);
/**
 	*  \param min minimum allowed value for the generated element
    *  \param max maximum allowed value for the generated element
    *  \return returns a random integer in [min,max]
*/
	int Randint (int min, int max); ///< Returns a random integer between [min,max]

	int** calcolaMatC45(int&,int&,int*&);
	double** generaNewTr(double** inOut,int NumInOut,int* indici,int RealAtt);

	void seleziona();
	double H(int indice,frbs& fis,double** inOutCamp, int numCamp);
	double HH(int indice1,int indice2, frbs& fis,double** inOutCamp, int numCamp);
	void riduciVariabili(vector<int> indiciAtt); //, int*&);


	void costruiscArffFile(const char*);
	char* calcolasets(int);
	void cambiaVariabili(int RealAtt,int* indiciAtt, int*&);

	float Rand();  ///< Generates a random number in [0,1]

/**
 	*  \param a first number
    *  \param b  second number
    *  \return the minimu value between a and b
*/
	int minN(int a,int b);///< returns the minimun value between a and b

	int countNoZero(unsigned*,int);///< Counts the number of non-zeros in a vector

/**
 	*  \param METHOD name of the evolutionary algorithm
    *  \param ec_tra MSE error on the training set
    *  \param ec_tst MSE error on the test set
	*  \param num_rules number of rules
*/
	void WriteResults (char *METHOD, double ec_tra, double ec_tst, int num_rules,int comp); ///< Writes the training error, the test error and the number of rules of the most accurate solution


/**
 	*  \param mat1 name of the original dataset
    *  \param mat2 name of the reduced dataset
	*  \param num vector containing the indices of the pattern to select from mat1
	*  \param numRed dimension of vector num
*/
	void reduceDataSet(double** mat1, double** mat2, int* num, int numRed);//<Copies int mat2 the rows of mat1 corresponding to the indices contained in num
	int calcolaMI(double** MI, frbs& fis,vector <double>& valMI,vector<double>& indH,int* indici);
	void calcolaIndice(vector <int>& insF,vector <int>& insS,double** MI,vector<double>& indH,frbs& fis,vector <double>& valIndice,vector<double>&,vector<double>&);
	void calcolaErroreWM(vector <int>insS);
/**
 	*  \param nomefile1 name of configuration file
    *  \param nomefile2 name of the output file for checking the read of the configuration parameters
*/
	void inizializeVar(char* nomefile1,char* nomefile2,int); //,vector<int>,double*,double*); ///< Reads the configuration file nomefile1, inizializes the dataset matrices (training and test), and find the range of each linguistic variable

/**
 	*  \param fileInput name of configuration file
    *  \param filePar name of the output file for checking the read of the configuration parameters
*/
	void readFiles(const char* fileInput, const char* filePar,int); ///< Read the configuration file and inizializes the dataset matrices (training and test)
/**
 	*  \param fp name of the dataset file
    *  \param num returns the number of patterns
	*  \param nvar returns the number of variables
	*  \return the matrix dataset [num x nvar] allocating the memory space
*/
	double** readDataSet(FILE* fp,int& num, int& nvar); ///< Reads from file fp the dataset matrices (training and test)

	void deleteMat(int**,int);///< Deallocate matrix memory

	double calculateMean(double vett[], int dimVet); ///< Calculates the mean value of the elements in vector vett
	double standarDeviation(double* vet, int dimVet);  ///<Calculates the standard deviation of the elements in vector vett
	double covariance(double* vet1, double* vet2, int dimVet); ///<calculates the covariance between vet1 and vet2
	double quaredR (double* vet1, double* vet2, int dimVet); ///<Calculates the squared R coefficient between vet1 and vet2

/**
 	*  \param vext vector to be printed
    *  \param dim  dimension of the vector
    *  \param os  output stream
*/

bool trovaVett(int* vett,int dim,int num);

bool trovaVett(unsigned* vett, unsigned** mat, int num);
int trovaUguali(int index, unsigned* vett, unsigned** mat, int num);

template <class Tipo>
void order(Tipo* vett, int n)
{	bool trovato=false;
	Tipo appo;
	for (int i=0;i<n-1 && !trovato;i++)
	{	trovato=true;
		for (int j=n-1;j>=i+1;j--)
		  if (vett[j]<vett[j-1])
		  {	appo=vett[j];
			vett[j]=vett[j-1];
			vett[j-1]=appo;
			trovato=false;
		  }
	}
}

template <class Tipo>
int trovaMax(Tipo* vett,int dim)
{	int indMas=0;
	Tipo mas=vett[0];
	for (int i=1;i<dim;i++)
		if (vett[indMas]<vett[i])
			indMas=i;
	return indMas;
}


/******************************************************************
*******************************************************************
Cecks to n-dimensional vectors of objectives to see if they are identical
returns 1 if they are, 0 otherwise
*****************************************************************
*******************************************************************/
template <class Tipo>
int equal(Tipo *first, Tipo *second, int n)
{	int obj = 0;
	do
    {	if(*first!=*second)
			return(0);
      *first++;
      *second++;
      obj++;
    }
  	while(obj < n);
  	return(1);
}


template <class Tipo>
void printVect(Tipo* vect,int dim, ostream& os)///< Prints the vector vect of dimension dim on the specified stream os
{
     for (int i=0;i<dim;i++)
         os<<vect[i]<<('\t');
     os<<endl;
}

/**
 	*  \param mat matrix to be printed
    *  \param numRows  number of rows
	*  \param numCol  number of columns
    *  \param os  output stream
*/
template <class Tipo>
void printMat(Tipo**mat,int numRows, int numCol, ostream& os)///< Prints the matrix mat on the specified stream os
{	for (int i=0;i<numRows;i++)
	{	for(int j=0;j<numCol;j++)
			os<<mat[i][j]<<'\t';
        os<<endl;
     }
}


/**
 	*  \param vettore vector to be ordered
    *  \param n  dimension of the vector
*/
template <class T>
void ordina(T* vettore,int n) ///< Orders a vector in ascendent way
{
	bool ordinato=false;
	double appo;
	for (int i=0;i<n-1 && !ordinato ;i++)
	{	ordinato = true;
		for (int j = n-1; j >= i+1; j--)
			if(vettore[j] < vettore[j-1])
			{	appo=vettore[j];
				vettore[j]=vettore[j-1];
				vettore[j-1]=appo;
				ordinato = false;
			}
	}
}

/**
 	*  \param vect vector where looking for the element
	*	\param val element to be searched
    *  \param dim  dimension of the vector
	*  \return true if the value specified is present in the vector
*/
template<class T, class G>
int find(G* vect, T val ,int dim) ///<  Searches for element in a vector
{
     for (int i=0;i < dim; i++)
         if (vect[i]==val)return 1;
     return 0;
}


#endif /* SUPPORTCHROM_H_ */
