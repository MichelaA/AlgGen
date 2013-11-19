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
#include "supportchrom.h"
#include <fstream>
#include <string>
#include <sstream>



	int trovaUguali(int index, unsigned* vett, unsigned** mat, int num) {
		for (int i = index; i < num; i++)
			if (index != i) {
				int j = 0;
				for (; j < numVar && mat[i][j] == vett[j]; j++)
					;
				if (j == numVar)
					return i;
			}

		return -1;
	}

bool trovaVett(unsigned* vett, unsigned** mat, int num)
{
		for (int i = 0; i < num; i++) {
			int j = 0;
			for (; j < numVar && mat[i][j] == vett[j]; j++)
				;
			if (j == numVar)
				return true;
		}
		return false;
	}


bool trovaVett(int* vett,int dim,int num)
{
	for (int i=0;i<dim;i++)
		if (vett[i]==num)
			return true;
	return false;
}

char* calcolasets(int numeroPart)
{
	char* sets=new char[numeroPart*2+2];
	sets[0]='{';
	sets[numeroPart*2]='}';
	for (int i=1;i<numeroPart*2;i++)
		if (i%2==0)
			sets[i]=',';
		else
			sets[i]=(char)(((int)'0')+((i+1)/2));

	sets[numeroPart*2+1]='\0';

	//cout<<sets<<endl;
	return sets;
}

//costruisco il file in formato ARFF per WEKA

void costruiscArffFile(const char* nomefile)
{
	ofstream fout (nomefile);
	double** intervalli = new double* [numVar-1];
	for (int i=0;i<numVar-1;i++)
		intervalli[i]=new double[numParts[i]-1];

	//Costruisco gli intervalli sulle variabili di ingresso per discretizzare
	double diff;
	for (int i=0;i<numVar-1;i++)
	{	diff=(maxVal[i]-minVal[i])/(double)(2*(numParts[i]-1));
		intervalli[i][0]=diff+minVal[i];
		for (int j=1;j<numParts[i]-1;j++)
			intervalli[i][j]=2*diff+intervalli[i][j-1];
	}

	fout<<"@relation RR"<<endl;

	for (int i=0;i<numVar-1;i++)
	{	fout<<"@attribute \'"<<i+1<<"\' {";
		for (int j=0;j<numParts[i]-1;j++)
				fout<<j+1<<',';
		fout<<numParts[i]<<'}'<<endl;
	}

	//sets=calcolasets(MAXPARTO);
	fout<<"@attribute \'Classe\' {";
	for (int i=0;i<numParts[numVar-1]-1;i++)
		fout<<i+1<<',';
	fout<<numParts[numVar-1]<<'}'<<endl;;
	fout<<"@data"<<endl;

	//Discretizzo il file di training e lo salvo
	int indice;
	for (int i=0;i<numPatterTr;i++)
	{	for (int j=0;j<numVar-1;j++)
		{	for (indice=0;indice<numParts[j]-1 && inOutTr[i][j]>intervalli[j][indice];indice++);
			if (indice==numParts[j]-1)
				fout<<numParts[j]<<'\t';
			else
				fout<<indice+1<<'\t';
		}
		fout<<inOutTr[i][numVar-1]<<endl;
	}

	for (int i=0;i<numVar-1;i++)
		delete[] intervalli[i];
	delete[] intervalli;

}

bool controllaPresenzaClassi(int** mat,int Nreg)
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



int** calcolaMatC45(int& numRegole,int& RealNumAtt, int*& indici)
{
	string nometest(file_test);

	int** newMat;
	int num0;
	int** matrice;



	int indice=nometest.find_last_of('/');
	int indice2=nometest.find_last_of('.');
	int numAtt;
	string nomeWeka=nometest.substr(indice+1,indice2-indice-1);
	nomeWeka=nomeWeka.append(".arff");
	fstream fin;

	costruiscArffFile(nomeWeka.c_str());
	std::ostringstream convert;
	string iterazione;
	string nomefile("file");

	string comando("java -classpath /usr/local/bin/WEKA/weka1.jar weka.classifiers.trees.J48Mod -P file.txt -C 0.25 -M 1 -t ");
	//string comando("java -classpath /usr/local/bin/WEKA/weka1.jar weka.classifiers.trees.J48Mod -P file.txt -U -M 1 -t ");
	comando=comando.append(nomeWeka);
	comando=comando.append("> elimina.txt");
//	cout<<comando.c_str()<<endl;
	int di=system (comando.c_str());

	fin.open("file.txt",ios::in);
	if (!fin)
	{	cout<<"Errore nell'apertura del file"<<endl;
		exit(1);
	}
	fin>>numRegole>>numAtt;
	if (numAtt<=2 || numRegole<2) //non applico pruning
	{	fin.close();
		fin.clear();
		comando="java -classpath /usr/local/bin/WEKA/weka1.jar weka.classifiers.trees.J48Mod -P file.txt -U -M 1 -t ";
		comando=comando.append(nomeWeka);
		comando=comando.append("> elimina.txt");
		di=system (comando.c_str());
		fin.open("file.txt",ios::in);
		if (!fin)
		{	cout<<"Errore nell'apertura del file"<<endl;
			exit(1);
		}
		fin>>numRegole>>numAtt;
		indici=new int[numAtt];
		for (int i=0;i<numAtt;i++)
			indici[i]=0;
		newMat=new int*[numRegole];
		num0=0;
		for (int i=0;i<numRegole;i++)
		{	newMat[i]=new int[numAtt];
			for (int j=0;j<numAtt;j++)
			{	fin>>newMat[i][j];
				if (newMat[i][j]!=0 && indici[j]==0)
				{	indici[j]=1;
					num0++;
				}
			}
		}
	}
//}
	else
	{	indici=new int[numAtt];
		for (int i=0;i<numAtt;i++)
			indici[i]=0;
		newMat=new int*[numRegole];
		num0=0;
		for (int i=0;i<numRegole;i++)
		{	newMat[i]=new int[numAtt];
			for (int j=0;j<numAtt;j++)
			{	fin>>newMat[i][j];
				if (newMat[i][j]!=0 && indici[j]==0)
				{	indici[j]=1;
					num0++;
				}
			}
		}
		if (!controllaPresenzaClassi((int**)newMat,numRegole))
		{	fin.close();
			fin.clear();

			delete[] indici;
			for (int i=0;i<numRegole;i++)
				delete[] newMat[i];
			delete[] newMat;

			comando="java -classpath /usr/local/bin/WEKA/weka1.jar weka.classifiers.trees.J48Mod -P file.txt -U -M 1 -t ";
			comando=comando.append(nomeWeka);
			comando=comando.append("> elimina.txt");
			di=system (comando.c_str());
			fin.open("file.txt",ios::in);
			if (!fin)
			{	cout<<"Errore nell'apertura del file"<<endl;
				exit(1);
			}
			fin>>numRegole>>numAtt;
			indici=new int[numAtt];
			for (int i=0;i<numAtt;i++)
				indici[i]=0;
			newMat=new int*[numRegole];
			num0=0;
			for (int i=0;i<numRegole;i++)
			{	newMat[i]=new int[numAtt];
				for (int j=0;j<numAtt;j++)
				{	fin>>newMat[i][j];
					if (newMat[i][j]!=0 && indici[j]==0)
					{	indici[j]=1;
						num0++;
					}
				}
			}
		}
	}
	RealNumAtt=num0;
	matrice=new int*[numRegole];

	for (int i=0;i<numRegole;i++)
	{	matrice[i]=new int[RealNumAtt];
		int k=0;
		for (int j=0;j<numAtt;j++)
			if (indici[j]!=0)
				matrice[i][k++]=newMat[i][j];
	}

	fstream fpObj;
	fpObj.open("NumVar.txt",ios::out | ios::app);
	fpObj<<RealNumAtt-1<<endl;
	fpObj.close();

	for (int i=0;i<numRegole;i++)
		delete[] newMat[i];
				delete[] newMat;
	return matrice;
}
/*
double** generaNewTr(double** inOut,int NumInOut,int* indici,int RealAtt)
{
	double ** inOutTrNew=new double*[NumInOut];
	for (int i=0;i<NumInOut;i++)
	{	inOutTrNew[i]=new double[RealAtt];
		int k=0;
		for (int j=0;j<maxTerms;j++)
		{	if (indici[j]!=0)
			{	inOutTrNew[i][k]=inOut[i][j];
				k++;
			}
		}
		inOutTrNew[i][RealAtt-1]=inOut[i][maxTerms];
	}
	return inOutTrNew;
}
*/


void cambiaVariabili(int RealAtt,int* indici,int*& part)
{
	fstream fout;
	fout.open("C45.txt",ios::out | ios::app);
	fout<<"Variabili non Selezionate: ";
	for (int i=0;i<numVar-1;i++)
		if (indici[i]==0)
			fout<<i+1<<' ';
	fout<<endl<<endl;

	int* partNew = new int[numVar];
	double* minValNew =new double[numVar];
	double* maxValNew =new double[numVar];

	double** TR=new double*[numPatterTr];
	double** TS=new double*[numPatterTs];

	int k=0;
	for (int i=0;i<numPatterTr;i++)
	{	k=0;
		TR[i]=new double[numVar];
	    if (i<numPatterTs)
	    	TS[i]=new double[numVar];
		for (int j=0;j<numVar-1;j++)
		{	if (indici[j] != 0)
			{	TR[i][k]=inOutTr[i][j];
				if (i<numPatterTs)
					TS[i][k]=inOutTs[i][j];
				k++;
			}
		}
		TR[i][k]=inOutTr[i][numVar-1];
		if (i<numPatterTs)
			TS[i][k]=inOutTs[i][numVar-1];
	}
	k=0;
	for (int i=0;i<numVar-1;i++)
		if (indici[i] !=0)
		{	partNew[k]=part[i];
			minValNew[k] = minVal[i];
			maxValNew[k] =maxVal[i];
			k++;
		}
	minValNew[k]=minVal[numVar-1];
	maxValNew[k]=maxVal[numVar-1];
	partNew[k]=part[numVar-1];

	for (int i=0;i<numPatterTr;i++)
	{	delete[] inOutTr[i];
		if (i<numPatterTs)
			delete[] inOutTs[i];
	}
	delete[] inOutTr;
	delete[] inOutTs;

	inOutTr=TR;
	inOutTs=TS;

	delete[] part;
	delete[] minVal;
	delete[] maxVal;

	part=partNew;
	minVal=minValNew;
	maxVal=maxValNew;
	numVar=k+1; //RealAtt-1;
	//numVar=maxTerms+1;
}


double HH(int indice1,int indice2, frbs& fis,double** inOutCamp, int numCamp)
	{
		int conta=0;
		double appMedia=0;
		double app1,app2;
		double valH=0;
		//cout<<"variabile "<<indice1<<endl;
		for (int i=0;i<numParts[indice1];i++)
		{	for (int j=0;j<numParts[indice2];j++)
			{	conta=0;
				appMedia=0;
				for (int k=0;k<numCamp;k++)
				{	app1=fis.evalMF(inOutCamp[k][indice1],i+1,indice1,numParts[indice1]);
					app2=fis.evalMF(inOutCamp[k][indice2],j+1,indice2,numParts[indice2]);
					//if (app1 && app2)
					//{	conta++;
						appMedia+=app1*app2;
					//}
				}
			//	cout<<"fuzzy set indice1="<<i<<"indice2 ="<<j<<" "<<appMedia<<endl;
				if (appMedia!=0)
					valH+=(appMedia/(double)numCamp)*log(((double)appMedia/(double)numCamp))/log(2);
			}


		}
		//cout<<"h="<<valH<<endl;
		return (-1)*valH;
	}

	/*double H(int indice,frbs& fis,double** inOutCamp, int numCamp)
	{

		int conta;
		double appMedia,app;

		double valH=0;
		for (int i=0;i<numParts[indice];i++)
		{	conta=0;
			appMedia=0;
			for (int j=0;j<numCamp;j++)
			{	app=fis.evalMF(inOutCamp[j][indice],i+1,indice,numParts[indice]);
				if (app!=0)
				{	conta++;
					appMedia+=app;
				}
			}
			if (conta!=0)
				valH+=(appMedia/(double)numCamp)*log((conta/(double)numCamp))/log(2);

		}

		return (-1)*valH;
	}*/

	double H(int indice,frbs& fis,double** inOutCamp, int numCamp)
	{
		//int conta;
		double appMedia;
		double valH=0;
		//cout<<"variabile "<< indice<<endl;
			for (int i=0;i<numParts[indice];i++)
			{	//conta=0;
				appMedia=0;
				for (int j=0;j<numCamp;j++)
				{	appMedia+=fis.evalMF(inOutCamp[j][indice],i+1,indice,numParts[indice]);
					/*if (app!=0)
					{	conta++;
						appMedia+=app;
					}*/
				}
				//cout<<"fuzzy set "<<i<<' '<<appMedia<<endl;
				if (appMedia!=0)
					valH+=(appMedia/(double)numCamp)*log((appMedia/(double)numCamp))/log(2);

			}
		//	cout<<"H="<<valH<<endl;
			return (-1)*valH;
		}





	void calcolaWMErr()
	{
		int** mat;
		int numRul=0;
		double (*tipoaggr)(double*, int, double);
		tipoaggr = minimo; // function for the AND operator
		bool nocop;
		double output,suma;
		fstream errTr,errTs;
		errTr.open("ErrTr.txt", ios::out | ios::app);
		errTs.open("ErrTs.txt", ios::out | ios::app);

		frbs fis(numVar);
		suma=0;
		fis.createFuzzySet(numParts, 0);
		mat=fis.generateWM(numRul, numParts, inOutTr, numPatterTr);
		for (int i=0;i<numPatterTr;i++)
		{	output = fis.calcOutput(inOutTr[i], *tipoaggr, numParts, mat, numRul, nocop);
			suma += 0.5 * pow(inOutTr[i][numVar - 1] - output, 2.0);
		}
				//cout<<"Tutte le feature:"<<endl<<"R="<<numRul<<" Err "<<suma/(double)numPatterTr<<endl;
				//cout<<"Errori"<<endl;
		errTr<<suma/(double)numPatterTr<<'\t';

		suma=0;
		for (int i=0;i<numPatterTs;i++)
		{	output = fis.calcOutput(inOutTs[i], *tipoaggr, numParts, mat, numRul, nocop);
			suma += 0.5 * pow(inOutTs[i][numVar - 1] - output, 2.0);
		}
		errTs<<suma/(double)numPatterTs<<'\t';
		fis.deleteFuzzySet(numParts);
	}

void calcolaWMrandom(vector<int> indiciF)
{	//frbs fis(numVar);
	//fis.createFuzzySet(numParts, 0);
	int* numRand=RandintDistinct(0,numVar-2,numVar-1);
	double (*tipoaggr)(double*, int, double);
	tipoaggr = minimo; // function for the AND operator
	vector <int> insFeat;
	insFeat.reserve(numVar-1);
	int** mat;
	double **inOutTrOld,**inOutTsOld;
	int* numPartOld;
	double* minOld,*maxOld;
	fstream feat,errTr,errTs,indFeat, nRegole;
	bool nocop;
	double output,suma;

	feat.open("feat.txt", ios::out | ios::app);
	errTr.open("ErrTr.txt", ios::out | ios::app);
	errTs.open("ErrTs.txt", ios::out | ios::app);
	nRegole.open("Regole.txt", ios::out | ios::app);

	int dim=indiciF.size();

	for (int i=0;i<dim;i++)
		feat<<indiciF[i]<<' ';
	feat<<endl;
	int numRul;
	int oldNumVar=numVar;
	inOutTrOld=inOutTr;
	inOutTsOld=inOutTs;
	minOld=minVal;
	maxOld=maxVal;
	numPartOld=numParts;
	for (int k=0;k<dim;k++)
	{	insFeat.push_back(indiciF[k]);
		riduciVariabili(insFeat);
		frbs fis1(numVar);
		fis1.createFuzzySet(numParts, 0);
		mat=fis1.generateWM(numRul, numParts, inOutTr, numPatterTr);
		suma=0;
		if (k==dim-1)
			nRegole<<numRul<<endl;
	   //stampaMatrice(mat,numRul,numVar);

		for (int j=0;j<numPatterTr;j++)
		{	output = fis1.calcOutput(inOutTr[j], *tipoaggr, numParts, mat, numRul, nocop);
			suma += 0.5 * pow(inOutTr[j][numVar - 1] - output, 2.0);
		}
		errTr<<suma/(double)numPatterTr<<'\t';
		suma=0;
		for (int i=0;i<numPatterTs;i++)
		{	output = fis1.calcOutput(inOutTs[i], *tipoaggr, numParts, mat, numRul, nocop);
			suma += 0.5 * pow(inOutTs[i][numVar - 1] - output, 2.0);
		}
		errTs<<suma/(double)numPatterTs<<'\t';

		fis1.deleteFuzzySet(numParts);

		for (int i=0;i<numPatterTr;i++)
		{	delete[] inOutTr[i];
			if (i<numPatterTs)
				delete[] inOutTs[i];
		}
		delete[] inOutTr;
		delete[] inOutTs;
		delete[] numParts;
		delete[] minVal;
		delete[] maxVal;

		for (int i=0;i<numRul;i++)
			delete[] mat[i];
		delete[] mat;
		inOutTr=inOutTrOld;
		inOutTs=inOutTsOld;
		minVal=minOld;
		maxVal=maxOld;
		numParts=numPartOld;
		numVar=oldNumVar;
		//insS.pop_back();
	}

	errTr<<endl;
	errTs<<endl;


}

	int calcolaMI(double** MI, frbs& fis,vector <double>& valMI,vector<double>& indH,int* indici)
	{	int indMax;
		double minH,HHval;
		double Husc=H(numVar-1,fis,inOutTr,numPatterTr); //entropia della variabile di uscita
		int partMax;
		frbs fis1(numVar);
		double maxEnt,ent,H_i;
		for (int i=0;i<numVar-1;i++)
		{	for (int j=2;j<6;j++)
			{	numParts[i]=j;
				fis1.createFuzzySet(numParts, 0);
				H_i=H(i,fis1,inOutTr,numPatterTr);
				if (H_i<Husc)
					minH=H_i;
				else
					minH=Husc;
				HHval=HH(i,numVar-1,fis1,inOutTr,numPatterTr);
				ent=(H_i+Husc-HHval); ///minH;
				//cout<<"entr X_i="<<H_i<<"  entr(XiXf)="<<HHval<<" MI="<<ent<<endl;
				if (j==2 || ent>maxEnt)
				{	maxEnt=ent;
					partMax=j;
				}
				fis1.deleteFuzzySet(numParts);
			}
			numParts[i]=partMax;
			MI[i]=new double[numVar];
			fis1.createFuzzySet(numParts, 0);
			indH[i]=H(i,fis1,inOutTr,numPatterTr); //entropia della i-esima variabile

			if (indH[i]<Husc)
				minH=indH[i];
			else
				minH=Husc;
			HHval=HH(i,numVar-1,fis1,inOutTr,numPatterTr);
			MI[i][i]=(indH[i]+Husc-HHval)/minH;  //Mutual information tra la variabile i e l'uscita
			MI[i][numVar-1]=(2*indH[i]-HH(i,i,fis1,inOutTr,numPatterTr))/indH[i];

			valMI[i]=MI[i][i];
			if (i==0 || MI[i][i]>MI[indMax][indMax])
				indMax=i;
			for (int j=0;j<numVar-1;j++)
				if (i!=j)
					MI[i][j]=-1;
			indici[i]=i;
			fis1.deleteFuzzySet(numParts);
		/*	MI[i]=new double[numVar-1];
			indH[i]=H(i,fis,inOutTr,numPatterTr); //entropia della i-esima variabile
			if (indH[i]<Husc)
				minH=indH[i];
			else
				minH=Husc;
			HHval=HH(i,numVar-1,fis,inOutTr,numPatterTr);
			MI[i][i]=(indH[i]+Husc-HHval)/minH;  //Mutual information tra la variabile i e l'uscita
			valMI[i]=MI[i][i];
			if (i==0 || MI[i][i]>MI[indMax][indMax])
				indMax=i;
			for (int j=0;j<numVar-1;j++)
				if (i!=j)
					MI[i][j]=-1;

			//vetMI[i]=MI[i][i];
			indici[i]=i;*/
		}
		fstream par;
		par.open("partizioni",ios::out|ios::app);
		for (int i=0;i<numVar;i++)
			par<<numParts[i]<<' ';
		par<<endl;
		return indMax;
	}

	double calcolaIns(vector <int> insS,double** MI,vector<double>& MIinsNum,vector<double>& MIinsDen)
	{	double num=0, den=0;
		int volte=0;

		/*for (int i=0;i<insS.size();i++)
			cout<<insS[i]<<' ';
		cout<<endl;

		for (int i=0;i<numVar-1;i++)
		{	for (int j=0;j<numVar;j++)
				cout<<MI[i][j]<<' ';
			cout<<endl;
		}*/


		for (int i=0;i<insS.size();i++)
		{	num+=MI[insS[i]][insS[i]];
			for (int j=i+1;j<insS.size();j++)
			{	volte++;
				//cout<<insS[i]<<' '<<insS[j]<<' '<<endl;
				den+=MI[insS[i]][insS[j]];

			}
			den+=MI[insS[i]][numVar-1];
			volte++;
		}

		num/=insS.size();
		den/=volte;
		MIinsNum.push_back(num);
		MIinsDen.push_back(den);
		return num/den;

	}

	void calcolaIndice(vector <int>& insF,vector <int>& insS,double** MI,vector<double>& indH,frbs& fis,vector <double>& valIndice,vector<double>& MIinsNum,vector<double>& MIinsDen)
	{	vector<int>::iterator it;
		double minV,somma, massimo=-1;
		int indMax;
		double indice=0;
		while(insF.size()>0)
		{	/*for (int k=0;k<insF.size();k++)
				cout<<insF[k]<<endl;*/
			indMax=-1;
			for (int i=0;i<insF.size();i++)
			{	somma=0;
				double MImedia=MI[insF[i]][insF[i]];
				double MIs=0;
				int numeri=0;

				for (int j=0;j<insS.size();j++)
				{
					if (MI[insF[i]][insS[j]]==-1)
					{	if (indH[insF[i]]<indH[insS[j]])
							minV=indH[insF[i]];
						else
							minV=indH[insS[j]];

						MI[insF[i]][insS[j]]=(indH[insF[i]]+indH[insS[j]]-HH(insF[i],insS[j],fis,inOutTr,numPatterTr))/minV;
						MI[insS[j]][insF[i]]=MI[insF[i]][insS[j]];
					}
					somma+=MI[insF[i]][insS[j]];
					MImedia+=MI[insS[j]][insS[j]];
				}

				/*for (int s=0;s<numVar-1;s++)
							{	for (int t=0;t<numVar;t++)
									cout<<MI[s][t]<<'\t';
								cout<<endl;
							}*/

				MIs=0;

				for (int j=0;j<insS.size();j++)
				{	//cout<<insS[j]<<endl;
					for (int k=j+1;k<insS.size();k++)
					{	MIs+=MI[insS[j]][insS[k]];
						numeri++;
					}
					MIs+=MI[insS[j]][numVar-1];
					numeri++;
					MIs+=MI[insS[j]][insF[i]];
					numeri++;
				}
				MIs+=MI[insF[i]][numVar-1];
				numeri++;
				MIs/=(double)numeri;

				MImedia/=((double)insS.size()+1);
				somma=somma/(double)insS.size();  //mutua informazione media insF(i) e l'insieme S
				//cout<<insF[i]<<' '<<MI[insF[i]][insF[i]]<<' '<<somma<<' '<<MI[insF[i]][insF[i]]-somma<<endl;
				indice=MI[insF[i]][insF[i]]-somma;
				//MIs=(MIs+somma)/(numeri+insS.size());

				indice=MImedia/MIs;

				if (indMax==-1 || massimo<indice)
				{	indMax=insF[i];
					massimo=indice;
				}
			}

			/*for (int i=0;i<numVar-1;i++)
			{	for (int j=0;j<numVar-1;j++)
					cout<<MI[i][j]<<' ';
				cout<<endl;
			}*/

			it=find(insF.begin(),insF.end(),indMax);
			insF.erase(it);
			insS.push_back(indMax);
			valIndice.push_back(massimo);
			double ind=calcolaIns(insS,MI,MIinsNum,MIinsDen);


		}
	}

void calcolaErroreWM(vector <int>insS)
{	int** mat;
	int numRul=0;
	double (*tipoaggr)(double*, int, double);
	tipoaggr = minimo; // function for the AND operator
	bool nocop;
	double output,suma;
	fstream errTr,errTs,regole;

	errTr.open("ErrTr.txt", ios::out | ios::app);
	errTs.open("ErrTs.txt", ios::out | ios::app);
	regole.open("Regole.txt",ios::out| ios::app);
	suma=0;
	/*fis.createFuzzySet(numParts, 0);
	mat=fis.generateWM(numRul, numParts, inOutTr, numPatterTr);

	for (int i=0;i<numPatterTr;i++)
	{	output = fis.calcOutput(inOutTr[i], *tipoaggr, numParts, mat, numRul, nocop);
		suma += 0.5 * pow(inOutTr[i][numVar - 1] - output, 2.0);
	}

	errTr<<suma/(double)numPatterTr<<'\t';

	suma=0;
	for (int i=0;i<numPatterTs;i++)
	{	output = fis.calcOutput(inOutTs[i], *tipoaggr, numParts, mat, numRul, nocop);
		suma += 0.5 * pow(inOutTs[i][numVar - 1] - output, 2.0);
	}
	errTs<<suma/(double)numPatterTs<<'\t';
	fis.deleteFuzzySet(numParts);

	for (int i=0;i<numRul;i++)
		delete[] mat[i];
	delete[] mat;
*/
	double **inOutTrOld,**inOutTsOld;
	int* numPartOld;
	double* minOld,*maxOld;
	inOutTrOld=inOutTr;
	inOutTsOld=inOutTs;
	minOld=minVal;
	maxOld=maxVal;
	numPartOld=numParts;
	int oldNumVar=numVar;

	vector <int> insFeat;
	insFeat.reserve(numVar-1);

	for (int k=0;k<insS.size();k++)
	{	insFeat.push_back(insS[k]);
		riduciVariabili(insFeat);

		frbs fis1(numVar);
		fis1.createFuzzySet(numParts, 0);
		mat=fis1.generateWM(numRul, numParts, inOutTr, numPatterTr);
		suma=0;


		for (int j=0;j<numPatterTr;j++)
		{	output = fis1.calcOutput(inOutTr[j], *tipoaggr, numParts, mat, numRul, nocop);
			suma += 0.5 * pow(inOutTr[j][numVar - 1] - output, 2.0);
		}
		errTr<<suma/(double)numPatterTr<<'\t';
		regole<<numRul<<'\t';

		suma=0;
		for (int i=0;i<numPatterTs;i++)
		{	output = fis1.calcOutput(inOutTs[i], *tipoaggr, numParts, mat, numRul, nocop);
			suma += 0.5 * pow(inOutTs[i][numVar - 1] - output, 2.0);
		}
		errTs<<suma/(double)numPatterTs<<'\t';

		fis1.deleteFuzzySet(numParts);

		for (int i=0;i<numPatterTr;i++)
		{	delete[] inOutTr[i];
			if (i<numPatterTs)
				delete[] inOutTs[i];
		}
		delete[] inOutTr;
		delete[] inOutTs;

		delete[] numParts;
		delete[] minVal;
		delete[] maxVal;

		for (int i=0;i<numRul;i++)
			delete[] mat[i];
		delete[] mat;

		inOutTr=inOutTrOld;
		inOutTs=inOutTsOld;
		minVal=minOld;
		maxVal=maxOld;
		numParts=numPartOld;
		numVar=oldNumVar;

	}
	errTr<<endl;
	errTs<<endl;
	regole<<endl;
}


	void seleziona()
	{
		vector <int>insS;
		vector <int> insF(numVar-1);
		vector <double> valIndice;
		vector<double> indH(numVar-1);
		vector <double> valMI(numVar);
		vector<int>::iterator it;

		vector <double>MIinsNum;
		vector <double>MIinsDen;

		int* indici =new int[numVar-1];
		double** MI = new double*[numVar-1];

		fstream val,feat,errTr,errTs,mi,MIinsiemeNum,MIinsiemeDen;

		val.open("valori.txt", ios::out | ios::app);
		mi.open("inMI.txt", ios::out | ios::app);
		feat.open("feat.txt", ios::out | ios::app);
		MIinsiemeNum.open("MIinsNum.txt", ios::out | ios::app);
		MIinsiemeDen.open("MIinsDen.txt", ios::out | ios::app);

		frbs fis(numVar);
		fis.createFuzzySet(numParts, 0);
		for (int i=0;i<numVar-1;i++)
			insF[i]=i;



		insS.reserve(numVar-1);
		valIndice.reserve(numVar-1);
		//valMI.reserve(numVar);

		int indMax=calcolaMI(MI,fis,valMI,indH,indici);

		//fis.deleteFuzzySet(numParts);
		it=find(insF.begin(),insF.end(),indMax);
		insF.erase(it);
		insS.push_back(indMax);

		MIinsNum.push_back(MI[indMax][indMax]);
		MIinsDen.push_back(1);
		valIndice.push_back(MI[indMax][indMax]);

		fis.createFuzzySet(numParts,0);
		calcolaIndice(insF,insS,MI,indH,fis,valIndice,MIinsNum,MIinsDen);
		fis.deleteFuzzySet(numParts);

		for (int i=0;i<insS.size();i++)
		{	feat<<insS[i]<<'\t';
			val<<valIndice[i]<<'\t';
			mi<<valMI[insS[i]]<<'\t';
			MIinsiemeNum<<MIinsNum[i]<<'\t';
			MIinsiemeDen<<MIinsDen[i]<<'\t';
		}
		feat<<endl;
		val<<endl;
		mi<<endl;
		MIinsiemeNum<<endl;
		MIinsiemeDen<<endl;


		calcolaErroreWM(insS);


}


void riduciVariabili(vector<int> indici)
{
	fstream fout;
	fout.open("C45.txt",ios::out | ios::app);
	fout<<"Variabili Selezionate: ";
	for (int i=0;i<indici.size();i++)
		//if (indici[i]==0)
			fout<<indici[i]<<' ';
	fout<<endl<<endl;


	int* partNew = new int[numVar];
	double* minValNew =new double[numVar];
	double* maxValNew =new double[numVar];

	double** TR=new double*[numPatterTr];
	double** TS=new double*[numPatterTs];

	int k=0;
	for (int i=0;i<numPatterTr;i++)
	{	k=0;
		TR[i]=new double[numVar];
	    if (i<numPatterTs)
	    	TS[i]=new double[numVar];
		for (int j=0;j<numVar-1;j++)
		{	if (find(indici.begin(),indici.end(),j)!= indici.end())
			{	TR[i][k]=inOutTr[i][j];
				if (i<numPatterTs)
					TS[i][k]=inOutTs[i][j];
				k++;
			}
		}
		TR[i][k]=inOutTr[i][numVar-1];
		if (i<numPatterTs)
			TS[i][k]=inOutTs[i][numVar-1];
	}
	k=0;
	for (int i=0;i<numVar-1;i++)
		if (find(indici.begin(),indici.end(),i)!= indici.end())
		{	partNew[k]=numParts[i];
			minValNew[k] = minVal[i];
			maxValNew[k] =maxVal[i];
			k++;
		}
	minValNew[k]=minVal[numVar-1];
	maxValNew[k]=maxVal[numVar-1];
	partNew[k]=numParts[numVar-1];

	/*for (int i=0;i<numPatterTr;i++)
	{	delete[] inOutTr[i];
		if (i<numPatterTs)
			delete[] inOutTs[i];
	}*/
	//delete[] inOutTr;
	//delete[] inOutTs;

	inOutTr=TR;
	inOutTs=TS;

	//delete[] part;
//	delete[] minVal;
//	delete[] maxVal;

	numParts=partNew;
	minVal=minValNew;
	maxVal=maxValNew;
	numVar=k+1; //RealAtt-1;
	//numVar=maxTerms+1;
}


void stampaMatrice(int** mat,int r,int c)
{	cout<<endl;
	for (int i=0;i<r;i++)
	{	for (int j=0;j<c;j++)
			cout<<mat[i][j]<<'\t';
		cout<<endl;
	}
	cout<<endl;

}

void stampaVett(int* vett,int n)
{
	cout<<endl;
	for (int i=0;i<n;i++)
		cout<<vett[i]<<'\t';
	cout<<endl;
}


/******************************************************************
*******************************************************************
Used in quicksort to arrange element in ascending order
******************************************************************
*******************************************************************/
int compareDec(const void *_a, const void *_b) {
    //used in quicksort function
        int *a, *b;

        a = (int *) _a;
        b = (int *) _b;

        return (*b - *a);
}


int compareCre(const void *_a, const void *_b) {
    //used in quicksort function
        int *a, *b;

        a = (int *) _a;
        b = (int *) _b;

        return (*a - *b);
}






/******************************************************************
*******************************************************************
Compares two n-dimensional vectors of objective values for minimization
problems
returns 1 if first dominates second, -1 if second dominates first,
and 0 otherwise
******************************************************************
*******************************************************************/
int compare_min(double *first, double *second, int n)
{
  int obj = 0;
  int deflt = 0;
  int current;

  do
  {	if(*first < *second)
		current = 1;
    else
		if(*second < *first)
			current = -1;
      	else
			current = 0;
    if((current)&&(current==-deflt))
		return(0);

    if(current!=0)
		deflt = current;
	obj++;
    *first++;
    *second++;
	}
	while(obj < n);

	return(deflt);
}




unsigned* RandintRule(int* max,int num)
{     unsigned* vect=new unsigned[num];

     for (int i=0;i<num;i++)
     {		  if (i!=num-1)
			vect[i]=Randint(0,max[i]);
		  else
			vect[i]=Randint(1,max[i]);
		//cout<<vect[i]<<'\t';
		  }
		//cout<<endl;
	  return vect;
}

/******************************************************************
*******************************************************************
Returns num distinct random integer between [min,max] allocating
the space for the return vector
	- min minimum value of each generated element
    - max maximum value of each generatedelement
	- num number of elements to be generated
******************************************************************
*******************************************************************/
int* RandintDistinct(int min,int max,int num)
{
     int* vect=new int[num];
     int i,val,app;
     for (i=0;i<num;i++)
         vect[i]=-1;
     i=0;
     while(i<num)
	{	val=Randint(min,max);
        app=find(vect,val,num);
        if(app==0)
		{	vect[i]=val;
            i++;
		}
    }
    return vect;
}

int* RandintDistinctOrd(int min,int max,int num)
{
     int* vect=new int[num];
     int i,val,index,app;
     for (i=0;i<num;i++)
         vect[i]=-1;
     i=0;
     while(i<num)
	{	val=Randint(min,max);
        app=find(vect,val,num);
        if(app==0)
		{	vect[i]=val;
            i++;
		}
    }

	qsort(vect,num, sizeof(int), &compareCre);


    return vect;
}

/******************************************************************
*******************************************************************
Returns a random integer between [min,max]
	- min minimum value of each generated element
    - max maximum value of each generatedelement
******************************************************************
*******************************************************************/
 int Randint (int min, int max)
{ 	int numero=min + (max-min+1) * Rand();
	if (numero>=max)
		numero=max;
	return numero;
}


int RandintSuper(int min, int max)
{
	int numero=min + (max-min+1) * Rand();
	if (numero>=max)
		numero=max;
	return numero;
}






/******************************************************************
*******************************************************************
Generates a random number in [0-1)
******************************************************************
*******************************************************************/
float Rand()
{	Seed= ( (Seed * PRIME) & MASK);
	return (Seed  * 0.4656612875e-9) ;
}


/******************************************************************
*******************************************************************
Counts the number of non-zeros in a vector
******************************************************************
*******************************************************************/
int countNoZero(unsigned*vect,int dim)
{	int cont=0;
    for (int i=0;i<dim;i++)
         if (vect[i]!=0) cont+=1;
	return cont;
}

/******************************************************************
*******************************************************************
Returns the minimun value between a and b
******************************************************************
*******************************************************************/
int minN(int a,int b)
{ return (a < b ? a : b);}


/******************************************************************
*******************************************************************
Writes the training error, the test error and the number of rules
of the most accurate solution
******************************************************************
*******************************************************************/
void WriteResults (char *METHOD, double ec_tra, double ec_tst, int num_rules, int comp)
{	char file[100];
	//FILE *fp;

	sprintf (file,"%scomunR.txt",METHOD);
	fstream fp;
	fp.open(file, ios::out | ios::app);
	if (!fp)
	{	cout<<"Output file open error"<<endl;
		exit(1);
	}
	fp<<num_rules<<endl;
	fp.close();

	sprintf (file,"%scomunTRA.txt",METHOD);

	fp.open(file, ios::out | ios::app);
	if (!fp)
	{	cout<<"Output file open error"<<endl;
		exit(1);
	}
	fp<<ec_tra<<endl;
	fp.close();

	sprintf (file,"%scomunTST.txt",METHOD);
	fp.open(file, ios::out | ios::app);
	if (!fp)
	{	cout<<"Output file open error"<<endl;
		exit(1);
	}
	fp<<ec_tst<<endl;
	fp.close();


	sprintf (file,"%scomunCOM.txt",METHOD);
	fp.open(file, ios::out | ios::app);
	if (!fp)
	{	cout<<"Output file open error"<<endl;
		exit(1);
	}
	fp<<comp<<endl;
	fp.close();



}

/******************************************************************
*******************************************************************
Copies int mat2 the rows of mat1 corresponding to the indices
contained in num
******************************************************************
*******************************************************************/

void reduceDataSet(double** mat1, double** mat2, int* num, int numRed)
{	//printMat(mat1,numPatterTr, numVar, cout);
	for(int i=0;i<numRed;i++ )
		for (int j=0;j<numVar;j++)
			mat2[i][j]=mat1[num[i]][j];

}



/******************************************************************
*******************************************************************
Read the configuration file nomefile1, inizializes the dataset
matrices (training and test), and find the range of each linguistic
variable
******************************************************************
*******************************************************************/
void inizializeVar(char*nomefile1,char*nomefile2,int fold)//,vector <int> gran,double* minimoV,double* maxV)
{
	 readFiles(nomefile1,nomefile2,fold);	// Read the configuration file and inizializes

	 maxVal=new double[numVar];
     minVal=new double[numVar];

      for (int i=0;i<numVar;i++)	// Find the range of each linguistic variable
	  {		maxVal[i]=minVal[i]=inOutTr[0][i];
	   		for (int j=1;j<numPatterTr;j++)
			{	if (j<numPatterTs)
				{	if (inOutTs[j][i]<minVal[i])
						minVal[i]=inOutTs[j][i];
             	 	 if (inOutTs[j][i]>maxVal[i])
             	 		 maxVal[i]=inOutTs[j][i];

				}
            	if (inOutTr[j][i]<minVal[i])
            			minVal[i]=inOutTr[j][i];
            	if (inOutTr[j][i]>maxVal[i])
            		maxVal[i]=inOutTr[j][i];
			}
	  }
      numParts=new int[numVar];
      //seleziono la granularita delle partizioni indici degli obiettivi che regolano l'evoluzione

   /*   if (MAXPARTI[0]==49)
      {	for (int i=0;i<numVar;i++)
      	{  numParts[i]=gran[i];
      	   if (!(minimoV==0 && maxV==0))
      		{	minVal[i]=minimoV[i];
      			maxVal[i]=maxV[i];
      		}
      	}
//      	numParts[numVar-1]=gran[numVar-1];

      	for (int i=0;i<numVar;i++)
      		cout<<numParts[i]<<' ';
      	cout<<endl;
      	return;
      }*/

      if (MAXPARTI[1]!=0) //significa che nel file di configurazione specifico la granularit x ogni singola variabile
      {	  for (int i=0;i<numVar;i++)
      		  numParts[i]=MAXPARTI[i]-'0';
      }
      else //specifico solo quella della variabile di uscita e l massima uguale per tutti
      {		for (int i=0;i<numVar-1;i++)
    	  	  	numParts[i]=MAXPARTO;
      	  numParts[numVar-1]=MAXPARTI[0]-'0';
      }



}


/******************************************************************
*******************************************************************
Read the configuration file nomefile1, inizializes the dataset
matrices (training and test)
******************************************************************
*******************************************************************/
void readFiles(const char* fileInput, const char* filePar,int fold)
{
      FILE* fp;
      char file_tra[300],  input_file[300], output_file[300];
      char numero[2];

       /* Read input data */
      sprintf (input_file,"%s",fileInput);
      if ((fp = fopen(input_file, "r")) == NULL) {
         printf ("Input file error %s",input_file);
         exit(1);
      }
      fscanf (fp,FORMATO_ENT,VAR_ENT);
      fclose (fp);

      /* Read parameters */
      sprintf (output_file,"%s.pgl",filePar);
      if ((fp = fopen (output_file,"w")) == NULL) {
         printf("Output file error %s",output_file);
         exit(1);
      }
      fprintf (fp,FORMATO_SAL,VAR_SAL);
      fclose (fp);

      /* Read training set */
      sprintf(numero,"%i",fold+1);
      strcat(file_tra,numero);
      strcat(file_tra,".tra");
      if ((fp = fopen(file_tra, "r")) ==NULL){
         printf("Error opening file\n");
         exit(1);
      }
      inOutTr=readDataSet(fp,numPatterTr,numVar);

      fclose(fp);
      /* Read test set*/
      strcat(file_test,numero);
      strcat(file_test,".tst");
      if ((fp = fopen(file_test, "r")) ==NULL){
         printf("Error opening file\n");
         exit(1);
      }
      inOutTs=readDataSet(fp,numPatterTs,numVar);
	  Seed = (unsigned long) semilla;
      fclose(fp);
}

/******************************************************************
*******************************************************************
Read the configuration file nomefile1, inizializes the dataset
matrices (training and test)
******************************************************************
*******************************************************************/
double** readDataSet(FILE* fp,int& num, int& nvar){
     float app2;
     double** mat;

     fscanf(fp,"%i",&num);
     fscanf(fp,"%i",&nvar);
     mat=new double* [num];

     for (int i=0;i<num;i++){
         mat[i]=new double[nvar];
         for (int j=0;j<nvar;j++){
             fscanf(fp,"%f",&app2);
             mat[i][j]= (double) app2;
         }
     }
     return mat;
}

/******************************************************************
*******************************************************************
Deallocates the matrix mat
******************************************************************
*******************************************************************/
void deleteMat(int** mat,int num)
{	for (int i=0;i<num;i++)
    	delete[] mat[i];
     delete[] mat;
}



/******************************************************************
*******************************************************************
Calculates the mean value of the elements in vector vett
******************************************************************
*******************************************************************/
double calculateMean(double vett[], int dimVet)
{	double mean=0;
	for (int i=0;i<dimVet;i++)
		mean+=vett[i];
	return mean/dimVet;
}


/******************************************************************
*******************************************************************
Calculates the standard deviation of the elements in vector vett
******************************************************************
*******************************************************************/
double standarDeviation(double* vet, int dimVet){
       double sd=0;
       double mean = calculateMean(vet, dimVet);
       for (int i=0;i < dimVet; i++)
           sd+=(vet[i]-mean)*(vet[i]-mean);
       return sqrt(sd/dimVet);
       }


/******************************************************************
*******************************************************************
calculates the covariance between vet1 and vet2
******************************************************************
*******************************************************************/
double covariance(double* vet1, double* vet2, int dimVet)
{
       double cov=0;
       double mean1 = calculateMean(vet1, dimVet);
       double mean2 = calculateMean(vet2, dimVet);
       for (int i=0;i < dimVet; i++)
           cov+=(vet1[i]-mean1)*(vet2[i]-mean2);
       return (cov/dimVet);
}

/******************************************************************
*******************************************************************
Calculates the squared R coefficient between vet1 and vet2
******************************************************************
*******************************************************************/
double quaredR (double* vet1, double* vet2, int dimVet){
	//calculates the squared R coefficient between vet1 and vet2
       double cov = covariance(vet1,vet2,dimVet);
       double sd1 = standarDeviation(vet1,dimVet);
       double sd2 = standarDeviation(vet2,dimVet);

	   return (cov/(sd1*sd2))*(cov/(sd1*sd2));
       }


/******************************************************************
*******************************************************************
List support functions
******************************************************************
*******************************************************************/

void insert (list *node, int x)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to enter after a NULL pointer, hence exiting \n");
        exit(1);
    }
	temp=new list;

    temp->index = x;
    temp->child = node->child;
    temp->parent = node;
    if (node->child != NULL)
    {
        node->child->parent = temp;
    }
    node->child = temp;

}

/* Delete the node NODE from the list */
list* del (list *node)
{
    list *temp;
    if (node==NULL)
    {
        printf("\n Error!! asked to delete a NULL pointer, hence exiting \n");
        exit(1);
    }
    temp = node->parent;
    temp->child = node->child;
    if (temp->child!=NULL)
    {
        temp->child->parent = temp;
    }
    delete node;
    return (temp);
}







