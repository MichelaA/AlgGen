/*
 * frbs.h
 *
 *  Created on: May 20, 2010
 *      Author: miki
 */

#ifndef FRBS_H_
#define FRBS_H_


struct retta{
	float p;
	float q;
};

#include <fstream>

/**
* Class to implement an MRFBS, (fuzzyfication, inference and defuzzyfication) different granularity can be defined for each lingustic variable
*/
class frbs
{
	double*** partVar; 		///< three dimensional vector [maximumGranularity][numFuzzysets][4], for each granularity, for each fuzzyset contains the 4 points representing the fuzzy set. Each variable is defined in [0-1]
	double*** fuzzyset;
	int nVar; 			///< number of variables of input-output
	double** dontCare;	///< for each variable contains the minimum and maximum value of the dont Care term
	double** matAttxPesi;
	bool tutto0;
	int dimMatAtt;

	/**
 	*  \param att return the activation vector
    *  \param matReg  matrix of rules
	*  \param nReg	number of rules
	*  \param realPart vector containing the granularity of each variable
 	*  \param var input pattern
	*  \param tipoaggr implication operator (set always as min)
	*  \param fs contains the scaled and Tuned (if enabled) DB parameters
*/
	void nearestRule(double* att, int** matReg, int nReg, int* realPart, double* var,double (*tipoaggr)(double*,int,double)); ///< calculates the activation vector, choosing the two nearest rules of a non-covered input pattern or move the input pattern in order to be covered at least by two rules

/**
 	*  \param att contains the activation vector
    *  \param matReg  matrix of rules
	*  \param nReg	number of rules
	*  \param numPart granularity of the output variable
 	*  \param fs fuzzyfied input pattern
   	* \return  the output value
*/
	double WECOA(double* att,int** matReg,int nReg,int numPart,double** fs); ///< Calculates the output of the MFRBS using the centre of gravity method

/**
 	*  \param att contains the activation vector
    *  \param nReg	number of rules
	*  \param numPart granularity of the output variable
	*  \param matReg  matrix of rules
	*  \param fs fuzzyfied input pattern
   	* \return  the 4 fuzzy set parameters after the implication
*/
	double** calccons(double* att, int nReg, int numPart, int** matReg,double** fs);  ///< Calculate the 4 fuzzyset parameters for each rule

/**
 	* \param indice the fuzzy set number of the partition (2-->MAXPART)
    *  \param valMin the minimum value of the variable
	*  \param valMax the maximum value of the variable
    * \return  a  [indicex4] matrix each row contains the parameter the correspondent fuzzy set
*/
	double** genTriUnifPartition(int indice, double valMin,double valMax); ///< Generates a uniform fuzzy partition between maxVal and minVal

/**
 	*  \param fuzzyset returns the scaled and Tuned (if enabled) DB parameters
	*  \param realPart vector containing the granularity of each variable
	*  \param pwLT  parameters for the piecewise linear transormation
	*  \param scaling  parameters for the nonlinear transformation of the fuzzy sets (if activeted)
	*  \param tun  parameters for the 2tuples tuning
*/
   	//void calcolafuzzyset(int* realPart,float** pwLT,scalingFunction& scaling,float* tun,double*,bool,int*); ///Performs the scaling and tuning of the database parameters
   	void calcolafuzzyset(int* realPart,float** pwLT);
/**
 	*  \param pwLTi parameters for the piecewise linear transformation of the i-th variable
	*  \return the p and q coefficients for each linear segment
*/
	retta* calcolarette(float* pwLTi,int); ///< Calculates the coefficients p and q for each linear segment

	double normalDef(double*,int**,int,int,double**);
	double euclidea(double* var,int reg[],double rangemin[], double rangemax[],int*);

public:

	double evalMF(double xf,int Af,int f,int F);
/**
 	* \param nVariable number of variables
    */
	frbs(int nVariable);  ///< Constructor
	frbs(const frbs& fis);
	~frbs();
	//void createFuzzySet(int* ,float** , float* , scalingFunction& , double*,bool);
	void createFuzzySet(int* realPart,float** pwLT);
	void deleteFuzzySet(int*);
/**
 	* \param var input pattern
    *  \param tipoaggr implication operator (set always as min)
	*  \param realPart vector containing the granularity of each variable
	*  \param matReg  matrix of rules
	*  \param nReg	number of rules
	*  \param nocop return parameter, specifies if the pattern is covered by the RB (nocop=false)
	*  \param scaling  parameters for the nonlinear transformation of the fuzzy sets (if activeted)
	*  \param tun  parameters for the 2tuples tuning
	*  \param pwLT  parameters for the piecewise linear transormation
    * \return  the output of the FRBS
 	*/
	double calcOutput(double* var, double (*tipoaggr)(double*,int,double), int* realPart, int** matReg,int nReg,bool& nocop); ///< Calculate the output of the FRBS

	void faiMatAttivazione(double** inOutCamp,int numCamp, int* realPart, int** matReg,int nReg,double*&,bool);
	int calcolaClasse(double* pesi, int indicePunto,int**,int nReg);


/**
 	* \param varInput input pattern
    *  \param fuzzyset contains the scaled and Tuned (if enabled) DB parameters
	*  \param vettPart vector containing the granularity of each variable
    * \return  for each variable and for each fuzzy set contains the membership degree to the fuzzy set
 */
	double** fuzzyInput(double* varInput, int* vettPart); ///< Fuzzifies the input pattern
	double fuzzyInputCL(double* varInput, int* vettPart,int*);

/**
 	* \param gran granularity of the variable
    *  \param index fuzzyset index (0-->gran-1)
	* \return  the centroid of the fuzzyset specified by index of granularity gran
 */
	double getCentroids(int gran,int index); ///< returns the centroid of the fuzzyset number index of granularity gran

/**
 	* \param gran granularity of the variable
	* \return  Returns the vector of centroid of the MF of a linguistic variable range [0--1] of granularity gran
 */
	double* getVettCentroids(int gran); ///<Returns the vector of centroid of the MF of a linguistic variable range [0--1] of granularity gran

	int** generateWM(int& numRul, int* vettPart, double** pattern, int nPattern); //float** pwLT) //,float* tun, scalingFunction& scaling );
	bool getTutto0(){return tutto0;}

	double* calcolaPesi(int** mat,int NumR,int* vettPart, double** pattern, int nPattern,float** pwLT); //calcolo i pesi della matrice mat
	double calcolaPeso(int* rule, double** pattern, int nPattern,int* vettPart,float** pwLT); //calcolo i pesi della matrice mat
	int** generateWMCL(int& NumR, int* vettPart, double** pattern, int nPattern,float** pwLT);



};

	double minimo(double*,int, double);
	double AreaTrapecioX(double*, double);
	double AreaTrapecio(double*, double);


#endif /* FRBS_H_ */
