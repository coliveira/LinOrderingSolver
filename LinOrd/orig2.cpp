#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>

using namespace std;

const int    nItLimit = 5000;		//Iteration limit
const double beta  = 0.3;			//cardinality restriction
const double alpha = 1;
const double LARGENUM = 1E9;	

clock_t cpu1,cpu2;	
int		Opt_Val;

void	LocalSearch(int* pSol, const int, double** c);
void	LocalSearch2k(int* pSol, const int, double** c); //Entire neighborhood
int		compare (const void * a, const void * b);

void	printM(double** ptr, int nN);				//debugging
void	printv(void);								//debugging

struct list {
public:
	double*	value;			//Greedy function
	double* value0;			//Copy of Initial Greedy function
	int*	idx;			//Greedy order
	int*	idx0;			//Copy of Initial greedy order
	int		length;

    void  ConstructSolution(int* pSol, double** c);
    bool InputInstance(const int argc, const char* datafile, 
				   int & nN, double** & c);
    double getSolValue (int* _pSol, double** _c, bool isReverse);
    void printv();

};

list	 g_GRASPData;

#include <string>

int main(int argc, char** argv)
{
	srand (1); // (unsigned)time(NULL) ); // initialize random generator
	double currentValue = 0;
    double bestValue = currentValue;

    const char *fileName = argv[1];

    int problem_size = 0;
    double **costs;
    if (g_GRASPData.InputInstance(argc, fileName, problem_size, costs) == false ) 
    {
        printf("error reading instance\n");
        return 1;
    }

    std::string output_file(fileName);
    output_file += "_inbe_10.txt"; 	// alpha+beta = 0.3

    ofstream out_file(output_file.c_str());
	if (out_file.fail())	cout << "Error opening: " << output_file << endl;
	
	int *currentSolution = new int [problem_size];
	int *bestSolution = new int [problem_size];

    double tIt = 0;
    int nIt = 0;
	double averageValue = 0;
	while (nIt < nItLimit)
	{
		cpu1 = clock();
		g_GRASPData.ConstructSolution(currentSolution, costs);
		cpu2 = clock();

		currentValue = g_GRASPData.getSolValue(currentSolution, costs, false);  //true is reverse (min)

		if (currentValue > bestValue) {
			for (int i=0; i<problem_size; i++) bestSolution[i] = currentSolution[i];
			bestValue = currentValue;
		}

		nIt++;
		tIt += cpu2-cpu1;
		averageValue += currentValue;

		out_file << nIt << setw(8) << cpu2-cpu1 << setw(10) << tIt/nIt << setw(15) << (averageValue/nIt) << setw(10) << (Opt_Val-(averageValue/nIt))*100/Opt_Val << setw(15) << bestValue << setw(10) << (Opt_Val-bestValue)*100/Opt_Val << endl;

		if (nIt%250==0)
		{
			cout << setw(5) << nIt << ":  Current: "
				 << setw(15) << currentValue << "  Best MAX: "
				 << setw(15) << bestValue << endl;
		}
	}
	//Report Best Solution
	cout << "\nBEST PERMUTATION FOUND:\n";
	for (int i=0; i<problem_size; i++)
	{
		cout << setw(4) << bestSolution[i];
		if ((i+1)%10==0)
			cout << endl;
	}
	out_file.close();
	return 0;
}

void list::ConstructSolution(int* pSol, double** c)
{
    list *pList = this;
	int s,i,k;
	int nRCL;
	int nN = length;
	double limit;
	int* pRCL;

	for (k=0; k<length; k++)
	{
		limit = (value[0] - value[nN-1-k])*alpha;
		i=0;
		nRCL = 1;
		for (i=1; i<nN-k && value[i] >= limit; i++) 
        {
            nRCL++;
        }
		
		pRCL = new int [nRCL];

		for (i=0; i<nRCL; i++)
			pRCL[i] = idx[i];

		s = rand()%nRCL;

		pSol[k] = pRCL[s];

		value[pSol[k]] = -LARGENUM;				//Send to bottom of list
		for (i=0; i<pSol[k]; i++)						//Adapt sums
		{
			if (value[i] != -LARGENUM)
				value[i] -= ( c[i][ pSol[k] ] - c[ pSol[k] ][i] );
		}
        qsort(idx, length - k, sizeof(int), ::compare);	//Sort again

		if (k>=1)
			LocalSearch2k(pSol,k+1,c);

		delete pRCL;
	}

	for (i=0; i < length; i++)
	{
		value[i] = value0[i];
		idx[i] = idx0[i];
	}
	
}

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

// FUNCTION: InputInstance
bool list::InputInstance(const int argc, const char* datafile, 
				   int & nN, double** & c)
{
    if (argc != 2) 
    {
        printf("usage: %s <instanceFile>\n");
        return false;
    }

    ifstream file(datafile);
    char str[80];
    int  i,k;

    file.getline(str,10000);	  //Ignore first line (title).

    cout << endl;
    cout << "--------------------------------------------\n"
        << "GRASP for LOP\n"
        << "--------------------------------------------\n"
        << "SETTINGS\n"
        << "  Max #Iterations  : " << nItLimit << endl
        << "  alpha            : " << alpha << endl
        << "  beta             : " << beta << endl
        << "PROBLEM" << endl
        << "  Filename         : " << datafile << endl
        << "  Title            : " << str << endl;
    file >> str;				  //Read size of problem
    nN = atoi(str);	
    cout << "  Size             : " << nN << endl
        << "  Optimal Solution : " << Opt_Val << endl
        << "--------------------------------------------\n\n";
    c = new double* [nN];		  //Create matrix
    for (k=0; k<nN; k++)		  
        c[k] = new double [nN];	  

    list *pList = this;
    value = new double [nN];	//Store sums
    value0 = new double [nN];//Store copy of sums
    idx	 = new int [nN];	//Store sorted indeces
    idx0  = new int [nN];	//Store copy of sorted indeces
    length = nN;

    for (i=0; i<nN; i++)
        idx[i] = i;			//Initialize indeces

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
        value[i] = 0;			
        for (int j=i+1; j<nN; j++)
            value[i] += ( c[i][j] - c[j][i] );
    }
    value[nN-1] = 0;

    qsort(idx, nN, sizeof(int), ::compare);	

    //Make copy of indeces and sums
    for (i=0; i<nN; i++)
    {
        value0[i] = value[i];
        idx0[i] = idx[i];
    }
    return true;
}

// AUXILIARY FUNCTIONS
int compare (const void * a, const void * b) //Compare, used to sort.
{
    list *pList = &g_GRASPData;
	int a1 = *(int*)a;
	int a2 = *(int*)b;
	if (pList->value[a1] < pList->value[a2])
		return 1;
	else if (pList->value[a1] == pList->value[a2])
		return 0;
	else
		return -1;
}

double list::getSolValue (int* _pSol, double** _c, bool isReverse)
{
    list *pList = this;
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

void printM(double** ptr, int nN)
{
	cout << endl;
    for (int i=0; i<nN; i++) {
		for (int j=0; j<nN; j++)
			cout << setw(5) << ptr[i][j] << " ";
		cout << endl;
	}
}

void list::printv(void)
{
    list *pList = this;
	cout << endl;
	for (int j=0; j<pList->length; j++)
		cout << setw(2) << pList->idx[j] << setw(8) << pList->value[pList->idx[j]] << endl;
	cout << endl;
}
