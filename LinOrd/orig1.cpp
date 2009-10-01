#if 0
// A GRASP implementation for the Linear Ordering Problem
// second formulation, using ratios.
// Bruno Chiarini.

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>

using namespace std;

//--SETTINGS-------------------------------------------
const int    nItLimit = 5000;		//Iteration limit
const double beta  = 0.3;			//cardinality restriction
const double alpha = 0.3;
const double LARGENUM = 1E9;		

bool InputInstance(const int, const char*, int &, double** &);
void ConstructSolution( int*, double**);
void LocalSearch(int* pSol, const int, double** c);
void LocalSearch2k(int* pSol, const int, double** c); //Entire neighborhood
int  compare (const void * a, const void * b);
double getSolValue (int*, double**,bool );

void printM(double** ptr, int nN);				//debugging
void printv(void);								//debugging

//----------------------------------------------------------------------------
struct list {
public:
	double*	value;			//Greedy function
	double* value0;			//Copy of Initial Greedy function
	int*	idx;			//Greedy order
	int*	idx0;			//Copy of Initial greedy order
	int		length;
};
list*	 pList = new list;					//Pointer to ``greedy'' list
//----------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////
// MAIN																		//
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
	int	nIt = 0;							//Iteration Counter
	int nN;									//Size of problem
	double** c;								//cost matrix
	int*	 pSol, *pSolBest;					//current and best Solution
	double   solValue, solValueBest;		//Current and best solution's value
	
	srand ( (unsigned)time(NULL) );					//Initialize random generator
	solValue = solValueBest = 0;

	if ( InputInstance(argc,argv[1],nN,c) == false ) //Exit if input fails.
		return 0;								

	pSol = new int [nN];
	pSolBest = new int [nN];

	//SOLVE
	while (nIt < nItLimit)
	{
		//ConstructGreedyRandomizedSolution(Solution);
 		ConstructSolution(pSol,c);

		//LocalSearch(Solution);
		  //LocalSearch2k(pSol,nN,c);

		//UpdateSolution(Solution);
		solValue = getSolValue(pSol,c,false);  //true is reverse (min)
		if (solValue > solValueBest)
		{
			for (int i=0; i<nN; i++)
				pSolBest[i] = pSol[i];
			solValueBest = solValue;
		}
		nIt++;
		if (nIt%250==0)
		{
			cout << setw(5) << nIt << ":  Current: "
				 << setw(15) << solValue << "  Best MAX: "
				 << setw(15) << solValueBest 
				 << "  Best MIN: "
				 << setw(15) << getSolValue(pSolBest,c,true) << endl;
		}
	}
	//Report Best Solution
	cout << "\nBEST PERMUTATION FOUND:\n";
	for (int i=0; i<nN; i++)
	{
		cout << setw(4) << pSolBest[i];
		if ((i+1)%10==0)
			cout << endl;
	}
	
cout << "\n PROBLEM be75eec.mat\n";
	int pbe75eec[] ={7, 48, 49, 33, 38, 50, 11, 28, 25, 27, 41, 45, 37, 23, 
		             22, 20, 24, 26, 42, 43, 1, 5, 21, 39, 34, 19, 15, 18, 16,
					 17, 14, 32, 8, 9, 10, 35, 36, 30, 46, 29, 31, 13, 12, 3, 
					 44, 2, 40, 4, 6, 47};
	cout << "be75eec.mat Optimum"<< endl;
	cout << sizeof(pbe75eec)/sizeof(int) << endl;
	for (int i=0; i<50; i++)
	{
		pbe75eec[i]--;
		cout << setw(4) << pbe75eec[i];
		if ((i+1)%10==0)
			cout << endl;
	}
	cout << endl << "Value: 264490" << endl;
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// End of MAIN

//////////////////////////////////////////////////////////////////////////////
// FUNCTION: ConstructSolution
void ConstructSolution(int* pSol, double** c)
{
	int s,i,k;
	int nRCL;
	int nN = pList->length;
	double limit;
	int* pRCL;

	for (k=0; k<pList->length; k++)
	{
		//MakeRCL(RCL)
		limit = (pList->value[0] - pList->value[nN-1-k])*alpha;
		i=0;
		nRCL = 1;
		for (i=1; i<nN-k; i++)
		{
			if (pList->value[i] >= limit)
				nRCL++;
			else
				break;
		}
		
		pRCL = new int [nRCL];
		//if  (pList->length - k < nRCL)
		//	nRCL--;
		for (i=0; i<nRCL; i++)
			pRCL[i] = pList->idx[i];

		//SelectElementAtRandom(RCL)
		s = rand()%nRCL;

		//Solution = Solution U {s}
		pSol[k] = pRCL[s];

		//cout << "s    = " << s << endl;    //debugging
		//cout << "nRCL = " << nRCL << endl; //debugging
		//printv();							 //debugging
		//cout << "Sort\n";					 //debugging

		//AdaptGreedyFunction(s)
		pList->value[pSol[k]] = -LARGENUM;				//Send to bottom of list
		for (i=0; i<pSol[k]; i++)						//Adapt sums
		{
			if (pList->value[i] != -LARGENUM)
				pList->value[i] -= ( c[i][ pSol[k] ] - c[ pSol[k] ][i] );
		}
		qsort(pList->idx, pList->length/* - k*/, sizeof(int), compare);	//Sort again

		if (k>=1)
			LocalSearch2k(pSol,k+1,c);

		delete pRCL;
	}
	//Restore list
	for (i=0; i < pList->length; i++)
	{
		pList->value[i] = pList->value0[i];
		pList->idx[i] = pList->idx0[i];
	}
	
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTION: LocalSearch		//Only combinations with the last value added.
void LocalSearch(int* pSol, const int nN, double** c)
{
	double deltaVal = 0;
	int temp=-1;						//Should we consider stopping as soon as
	int i,k;							//we get an improvement.
	int SolBest[2];
	double deltaValBest = 0;
	bool stop = false;
	const int j = nN-1;					//node nN
	for (i=0; i<nN-1; i++)				//permutes with i \in {0,...,nN-1}
	{
		deltaVal = c[ pSol[j] ][ pSol[i] ] - c[ pSol[i] ][ pSol[j] ];
		for (k = i+1; k<j; k++)
		{
			deltaVal += c[ pSol[k] ][ pSol[i] ] - c[ pSol[i] ][ pSol[k] ];
			deltaVal += c[ pSol[j] ][ pSol[k] ] - c[ pSol[k] ][ pSol[j] ];
		}
		if (deltaVal > deltaValBest)	//Record best permutation so far.
		{
			SolBest[0] = i;
			SolBest[1] = j;
			deltaValBest = deltaVal;
		}
	} 

	//If improvement found, do permutation
	if (deltaValBest > 0)
	{
		temp = pSol[ SolBest[0] ];
		pSol[ SolBest[0] ] = pSol[ SolBest[1] ];
		pSol[ SolBest[1] ] = temp;
	}
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTION: LocalSearch2k		//Full neighborhood	
void LocalSearch2k(int* pSol, const int nN, double** c)
{
	double deltaVal = 0;
	int temp=-1;						//Should we consider stopping as soon as
	int i,j,k;							//we get an improvement.
	int SolBest[2];
	double deltaValBest = 0;
	bool stop = false;
	i = 0;
	do									//First node i
	{
		j=i+1;							//permutes with j \in {i+1,...,nN-1}
		do 
		{
			deltaVal = c[ pSol[j] ][ pSol[i] ] - c[ pSol[i] ][ pSol[j] ];
			for (k = i+1; k<j; k++)
			{
				deltaVal += c[ pSol[k] ][ pSol[i] ] - c[ pSol[i] ][ pSol[k] ];
				deltaVal += c[ pSol[j] ][ pSol[k] ] - c[ pSol[k] ][ pSol[j] ];
			}
			if (deltaVal > deltaValBest)	//Record best permutation so far.
			{
				SolBest[0] = i;
				SolBest[1] = j;
				deltaValBest = deltaVal;
			}
			j++;
		} while (j<nN);	   //while no not end.
		i++;
	} while (i<nN-1);
	//If improvement found, do permutation
	if (deltaValBest > 0)
	{
		temp = pSol[ SolBest[0] ];
		pSol[ SolBest[0] ] = pSol[ SolBest[1] ];
		pSol[ SolBest[1] ] = temp;
	}
}
//////////////////////////////////////////////////////////////////////////////
// FUNCTION: InputInstance
bool InputInstance(const int _argc, const char* datafile, 
				   int & nN, double** & c)
{
	if (_argc!=2)
		return false;
	else
	{
		//Read data
		ifstream file(datafile);
		char str[80];
		int  i,k;
		file.getline(str,10000);	  //Ignore first line (title).

		cout << endl;
		cout << "--------------------------------------------\n"
			 << "GRASP for LOP\n"
			 << "--------------------------------------------\n"
			 << "SETTINGS\n"
			 << "  Max #Iterations: " << nItLimit << endl
			 << "  alpha          : " << alpha << endl
			 << "  beta           : " << beta << endl
			 << "PROBLEM" << endl
		     << "  Filename       : " << datafile << endl
			 << "  Title          : " << str << endl;
			 

		file >> str;				  //Read size of problem
		nN = atoi(str);	
		
		cout << "  Size           : " << nN << endl
			 << "--------------------------------------------\n\n";

		c = new double* [nN];		  //Create matrix
		for (k=0; k<nN; k++)		  
			c[k] = new double [nN];	  

		pList->value = new double [nN];	//Store sums
		pList->value0 = new double [nN];//Store copy of sums
		pList->idx	 = new int [nN];	//Store sorted indeces
		pList->idx0  = new int [nN];	//Store copy of sorted indeces
		pList->length = nN;
		
		for (i=0; i<nN; i++)
			pList->idx[i] = i;			//Initialize indeces

		//Read data
		i = 0;
		for (k=0; k<nN*nN; k++)		    
		{
			file >> str;
			c[i][k%nN] = atoi(str);
			if ( (k+1)%nN == 0 ) 
				i++;
		}

		//Collapse matrix and obtain row sums
		for (i=0; i<nN-1; i++)			
		{
			pList->value[i] = 0;			
			for (int j=i+1; j<nN; j++)
				pList->value[i] += ( c[i][j] - c[j][i] );
		}
		pList->value[nN-1] = 0;
		
		//printM(c,nN);  //debugging
		//cout << endl;  //debugging
		//printv(); //debugging
		//Sort indeces
		qsort(pList->idx, nN, sizeof(int), compare);	
		//printv(); //debugging

		//Make copy of indeces and sums
		for (i=0; i<nN; i++)
		{
			pList->value0[i] = pList->value[i];
			pList->idx0[i] = pList->idx[i];
		}
		return true;
	}
}
//////////////////////////////////////////////////////////////////////////////
// AUXILIARY FUNCTIONS
int compare (const void * a, const void * b) //Compare, used to sort.
{
	int a1 = *(int*)a;
	int a2 = *(int*)b;
	if (pList->value[a1] < pList->value[a2])
		return 1;
	else if (pList->value[a1] == pList->value[a2])
		return 0;
	else
		return -1;
}
//----------------------------------------------------------------------------
double getSolValue (int* _pSol, double** _c, bool isReverse)
{
	double val = 0;
	int i,j;
	for (i=0; i < pList->length-1; i++)
	{
		for (j=i+1; j < pList->length; j++)
		{
			if (isReverse)
				val += _c[ _pSol[j] ][ _pSol[i] ];
			else
				val += _c[ _pSol[i] ][ _pSol[j] ];
		}
	}
	return val;
}

#endif