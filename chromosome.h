/*
 * chromosome.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef CHROMOSOME_H_
#define CHROMOSOME_H_

#define  _stampa ;


#include<cmath>
#include <iostream>
#include <cstdlib>
#include "supportchrom.h"
#include "frbs.h"

using namespace std;


/**
 * Class chromosome codifies the rule base (RB) and the database (DB) of an MFRBS
 */

struct gen{
		int index;
		bool* vectAnt;
	};

class chromosome
{

protected:

	unsigned** matR;

	int* realPart; 		///< codifies the granularity of each input-output variable as descibed in [1], each element is between 2 and the maximum number of partitions (MAXPART)
	float** pwLT; 		///< codifies the displacement of the core of each fuzzy set (except for the first and the last fuzzy set) for each variable. it is a (number_of_variable x (MAXPART-2)) matrix. Each element is between [0 1] and each row is ordered in a crescent way as described in [4]
	float** pwLT_lim[3];///< for each element of pwLT stores the minimum, median and maximum number allowed. It is exploited when we perform a constrained tuning
	//int* terms;  		///< codifies the maximum granularity for each variable
	int numRul;			///< current number of rules of the virtual RB
	int realNumRule; 	///< current number of rules of the concrete RB
	int numFeat;
	double comp;

	double* objTot;		///< contains the multi-objective fitness function
	int sizeObjTot;	   ///< number of specific measures calculated for each chromosome
	int* objIndex;   ///< contains the index of objectives to take into account for the evolution
	double* objAlg;  	///<containes the objective values that drive the genetic algorithm
	int sizeObjAlg; 	///< number of real objectives taken into account for the evolutionary process

	int** matConfusione;
	//pesi per la classificazione
	double* pesi;
	bool noCovered;

	//nsga-II
	int rank;			///< nsga2 ranking
	double crowd_dist;  ///< nsga2 crowding distance

	//paes
	int grid_loc;		///< index of the position of the chromosome in the paes grid

	double interp;

	gen* vettR;
	static unsigned int dimMatWM;
	static unsigned int** matWM;
	static int* indiciCL;

	void allocateSpace(); ///< allocates the space for all the variables of the chromosome except for tun which is allocated in inizializewithoutrule
	void deleteChrom();	///< deallocates all the memory of the chromosome
	void faiMatdaVec(); ////costruisco ma matrice matR a seconda dei valori di vattR
	void ordinaIndex(); //ordina gli indici di vettR
	void deleteDup();
	int deleteDup(unsigned** mat, int numR,int*);

	void ordinaxClasse(int** mat,int numR,int* index);
	void ordinaxPesi(int**,int,double*);
	void eliminaConliffitto(double* appopesi,double soglia, int** mat, int& dimMat);


	int scegliClasse();

	void generateIndici(int,int);
	void generateMatCL(frbs&);
	void generateIndiciCL(frbs&);
	void generateMatR(frbs&);

	bool controllaPresenzaClassi(int**,int);

	double calcolaErroreC45(frbs& fis, double** inOutCamp, int numCamp,bool test, unsigned int** matAppo,int);

public:
		chromosome(){}
        ~chromosome();

        void faiVecdaMat();
        void controllaPesi(chromosome& Par,frbs&);

       void controllaClassi(frbs& fis);
        void eliminaPesiNegativi();

        /**
    *  \param fis a fis containing the structures of the DB.
 	*/  void inizialize(frbs& fis); ///< Initializes a chromosome with random values adding numInitRules rules
		void inizializewithoutrule(); ///< Initializes a chromosome with random values for all the parameters without adding rule

		void inizializeWM(frbs&);

		chromosome(const chromosome&);
		chromosome& operator=(const chromosome&);
		void copyChrom(const chromosome&);///< Copy the chromosome source in *this

			//RB mutation operators
		/**
    	*  \param fis the fis is needed to extract centroids for converting the virtual RB in the concrete RB
 		*/
		void addRule(frbs& fis);///< Adds a rule to the virtual RB which is not present in both the virtual and the concrete RB
		void remRule();			///< Remove a random rule from the virtual RB
		void modRule();			///< Modifies a term in the virtual RB
		void elimFeatureWM();
		//DB mutation operators
		void MutacionGR();///< Mutation for the granularity part (modifies vector realPart)
		void MutacionPW(); ///<Mutation for the piecewise (modifies pwLT)
		void mutSelection(double); //mutation for rule selection


		//RB crossover
		friend void crossRB(chromosome&,chromosome&,chromosome&,chromosome&); ///< One-point crossover for the virtual RB (cc)
		friend void crossWM(chromosome&,chromosome&,chromosome&,chromosome&); ///< One-point crossover for the virtual RB (cc)

		friend void cross2(chromosome&,chromosome&,chromosome&,chromosome&);///< two-point crossover for the virtual RB (cc)
		//DB crossover
		friend void crossPart(chromosome&,chromosome&,chromosome&,chromosome&);		///< one point crossover for the granularity part (realPart)
		friend void CruceBLXPW(chromosome&,chromosome&,chromosome&,chromosome&);	///< BLX crossover for the piecewise part (pwLT)

		/**
		*  \param comp returns the complexity of the concrete RB
		*  \param fis the fis is needed to convert the virtual RB into the concrete RB and to calculate the output
    	*  \param corr return the correlation if activated
		*  \param conc if conc=true the mse is calculated using the concrete RB
		* \return the error on the training set
 		*/
		double* ECM(frbs& fis,double** inOutCamp,int numCamp,bool test) ;
		void generateRB_WM(frbs&, double**,int);
		void generateRuleWM(frbs&);
		void allocateSpaceWM(int,unsigned int**);

		void calcolateWMError(frbs& fis);
		void salvaDatiC45(frbs&);

	//	double* calcolaIndice(frbs& fis,double** inOutCamp, int numCamp);

		/**
		*  \param var input pattern
		*  \param mat concrete RB
		*  \param dim number of rules of the concrete RB
		*  \param fis the fis is needed to calculate the output
    	*  \param nocop returns true if the input pattern is not covered by the RB
		*  \param conc if conc=true the mse is calculated using the concrete RB
		* \return the output on the input pattern of the FRBS
 		*/
		double FLC(double* var,unsigned** mat, int dim,frbs& fis,bool& nocop); ///< Calculates the output of the FRBS given an input pattern and a RB

		/**
    	*  \param fis the fis is needed to extract centroids for converting the virtual RB in the concrete RB
		*  \param conc if conc=true the mse is calculated using the concrete RB
 		*/
		void evaluateChrom(frbs& , double** , int);///< Calculates the fitness functions for a chromosome for classification problems

		/**
	 	* \param matReg the virtual RB.
    	*  \param fis a fis containing the structures of the DB.
 		*/
		void convertimat(unsigned**,frbs&);///< Converts the matrix matReg of the virtual RB into matrix of the concrete DB

		/**
	 	* \param rule the virtual rule.
    	*  \param fis a fis containing the structures of the DB.
    	* \return the concrete rule
 		*/
		unsigned* convertiRule(unsigned* rule, frbs& fis);///< Converts the rule of a virtual RB into a rule of the concrete DB
		void countRule();  	///< Counts the number of rules of the virtual RB and updates numRul
		int deleteDup(unsigned** mat,int numR);
		int giveActive();///< Returns a vector that for each input linguistic label counts the number of rules associated to
		int giveActive(unsigned**, int);

		double evaluateDistPW();///Returns the interpretability measure for the piecewise transformation
		double evaluateDistPWC(double&,double*,frbs&);///Returns the interpretability measure for the piecewise transformation
		double evaluateDistPW(int);///<Returns the interpretability measure for variable var

		bool getNoCoverage(){return noCovered;}

		//Get functions
		float** getPwLT(){return pwLT;}
		float* getPwLT(int i){return pwLT[i];};
		unsigned** getMatR(){return matR;}
		int getSizeObjTot(){return sizeObjTot;}
		int getSizeObjAlg(){return sizeObjAlg;}
		double getObjAlg(int num){return objAlg[num];}
		double* getObjAlg(){return objAlg;}
		double getObjTot(int num){return objTot[num];}
		double* getObjTot(){return objTot;}
		int getNumRule(){return numRul;}
		int getRealNumRule(){return realNumRule;}
		int getNumVar(){return numVar;}
		double getCrowd_dist(){return crowd_dist;}
		int getGrid_loc(){return grid_loc;}
		int getRealPart(int index){return realPart[index];}
		int* getRealPart(){return realPart;}
		double getInterp(){return interp;}
		int getNumFeat(){return numFeat;}
		double* getPesi(){return pesi;}
		double getComp(){return comp;}
		gen* getVettR(){return vettR;};

		//Set functions
		void setRealPart(int* vett,int n){for(int i=0;i<n;i++)  realPart[i]=vett[i];}
		void setRank(int r){rank=r;}
		void setCrowd_dist(double d){crowd_dist=d;}
		void setGrid_loc(int gr){grid_loc=gr;}
		void setnumRul(int n){ numRul=n;}

		// Print functions
		void printChrom(frbs&);
		void stampaPW();
};

int compareChrom(const void* ,const void*); ///< Comparing functions used in quicksort to order chromosome by ascending accuracy

/*************************************************************
*************************************************************
One-point crossover
*************************************************************
**************************************************************/
template <class T>
void OnepointCross(T* par1,T* par2,T* suc1,T* suc2, int max)
{
	int num=Randint(1,max-1);
	int i;
	//crossing the partition cromosome
	for (i=0;i<num;i++)
	{	suc1[i]=par1[i];
		suc2[i]=par2[i];
	}
	for (;i<max;i++)
	{	suc1[i]=par2[i];
		suc2[i]=par1[i];
	}
}

/*************************************************************
*************************************************************
BLX crossover
*************************************************************
**************************************************************/
template <class T>
void CruceBLX(T* par1,T* par2,T* suc1,T* suc2, int numV,T min, T max)
{	double temp, px, py, I, alfa;
	alfa=0.5;
	for (int j=0; j<numV; j++)
	{	px=par1[j];  py=par2[j];
		if (px > py) {temp=px; px=py; py=temp;}
		I=py-px; px=px-I*alfa; py=py+I*alfa;
		if (px < min) px=min;
		if (py > max) py=max;
		suc1[j]=px + Rand()*(py-px);
		suc2[j]=px + Rand()*(py-px);
		if (suc1[j]<min)
			suc1[j]=min;
		if (suc1[j]>max)
			suc1[j]=max;
		if (suc2[j]<min)
			suc2[j]=min;
		if (suc2[j]>max)
			suc2[j]=max;
	}
}

#endif //CHROMOSOME_H

