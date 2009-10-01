#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string>

using namespace std;

const int    nItLimit = 5000;		//Iteration limit
const double beta  = 0.3;			//cardinality restriction
const double ALPHA = 1;
const double LARGENUM = 1E9;	

int		Opt_Val;

void	LocalSearch(int* pSol, const int, double** c);
void	LocalSearch2k(int* pSol, const int, double** c); //Entire neighborhood
int		compare (const void * a, const void * b);

class GraspData {
	double *greedy_value;
	double *copy_of_values;			//Copy of Initial Greedy function
	int *greedy_order;
	int *copy_of_order;			//Copy of Initial greedy order
    int *RCL;
	int length;
    double **costs;

public:

    void  ConstructSolution(int* pSol);
    bool InputInstance(const int argc, const char* datafile, int& problem_size);
    double getSolValue (int* _pSol, bool isReverse);
    void printv();
    bool AllocData(int problem_size);
    bool ReadInstanceValues(ifstream &file);
    double GreedyValue(int pos) { return this->greedy_value[pos]; }

    GraspData() 
        : greedy_value(0)
        , copy_of_values(0)
        , greedy_order(0)
        , copy_of_order(0)
        , RCL(0)
        , length(0)
        , costs(0) {}
    ~GraspData();

};

GraspData g_GRASPData;


int main(int argc, char** argv)
{
	srand (1); // (unsigned)time(NULL) ); // initialize random generator
	double currentValue = 0;
    double bestValue = currentValue;

    const char *fileName = argv[1];

    int problem_size = 0;
    if (g_GRASPData.InputInstance(argc, fileName, problem_size) == false ) 
    {
        printf("error reading instance\n");
        return 1;
    }

    std::string output_file(fileName);
    output_file += "_inbe_10.txt"; 	// ALPHA+beta = 0.3

    ofstream out_file(output_file.c_str());
	if (out_file.fail())
        cout << "Error opening: " << output_file << endl;	

	int *currentSolution = new int [problem_size];
	int *bestSolution = new int [problem_size];

    double tIt = 0;
    int nIt = 0;
	double averageValue = 0;
	while (nIt < nItLimit)
	{
        clock_t cpu1 = clock();
		g_GRASPData.ConstructSolution(currentSolution);
		clock_t cpu2 = clock();

		currentValue = g_GRASPData.getSolValue(currentSolution, false);  //true is reverse (min)

		if (currentValue > bestValue) {
			for (int i=0; i<problem_size; i++) bestSolution[i] = currentSolution[i];
			bestValue = currentValue;
		}

		nIt++;
		tIt += cpu2-cpu1;
		averageValue += currentValue;


		out_file << nIt; 
        out_file << setw(8);
        out_file << cpu2-cpu1;
        out_file << setw(10);
        out_file << tIt/nIt;
        out_file << setw(15);
        out_file << (averageValue/nIt);
        out_file << setw(10);
        out_file << (Opt_Val-(averageValue/nIt))*100/Opt_Val;
        out_file << setw(15);
        out_file << bestValue;
        out_file << setw(10);
        out_file << (Opt_Val-bestValue)*100/Opt_Val << endl;

		if (nIt%10==0)
		{
            //if (nIt == 250) assert(currentValue == 262758);
            //if (nIt == 250) assert(bestValue == 264455);
            printf("iteration: %d current: %d best %d \n", nIt, (int)currentValue, (int)bestValue);
			//cout << setw(5) << nIt << ":  Current: "
			//	 << setw(15) << currentValue << "  Best MAX: "
			//	 << setw(15) << bestValue << endl;
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

void GraspData::ConstructSolution(int* perm)
{
	double limit;

	for (int k=0, num_remaining=length; k<length; ++k, --num_remaining)
	{
		limit = (greedy_value[0] - greedy_value[num_remaining-1]) * ALPHA;
		
        int rcl_size = 1;
        for (int i=1; i<num_remaining && greedy_value[i] >= limit; i++)  rcl_size++;
		
		for (int i=0; i<rcl_size; i++) RCL[i] = greedy_order[i];

        int elem = perm[k] = RCL[rand() % rcl_size];  // select an element from RCL

		greedy_value[elem] = -LARGENUM;	// send position to bottom of list
		
        for (int i=0; i<elem; i++)       // update sums
			if (greedy_value[i] != -LARGENUM)
				greedy_value[i] -= ( costs[i][elem] - costs[elem][i] );

        qsort(greedy_order, num_remaining, sizeof(int), ::compare);	// update greedy order

		if (k > 0) LocalSearch2k(perm, k+1, costs);  // local improvement
	}

	for (int i=0; i < length; i++) greedy_value[i] = copy_of_values[i];
	for (int i=0; i < length; i++) greedy_order[i] = copy_of_order[i];
}

// LocalSearch. Checks only combinations containing the last value.
void LocalSearch(int* perm, const int nN, double** c)
{
	const int LAST_POS = nN-1;

    int pos;
	double best_delta = 0;
    int j = perm[LAST_POS];
	for (int i=0; i<LAST_POS; i++)	// permutes with i \in {0,...,nN-1}
	{
        double delta = c[j][perm[i]] - c[perm[i]][j];

        for (int k=i+1; k<LAST_POS; ++k) delta += c[perm[k]][perm[i]] - c[perm[i]][perm[k]];
		for (int k=i+1; k<LAST_POS; ++k) delta += c[j][perm[k]]       - c[perm[k]][j];
		
        if (delta > best_delta)	// record best permutation so far
		{
			pos = i;
			best_delta = delta;
		}
	} 

	if (best_delta > 0) // improved solution
	{
		int temp = perm[ pos ];
		perm[ pos ] = perm[ LAST_POS ];
		perm[ LAST_POS ] = temp;
	}
}

// perform search on full neighborhood	
void LocalSearch2k(int* perm, const int nN, double** c)
{
	int pos1, pos2;
	double best_delta = 0;
	bool stop = false;
	for (int i=0; i<nN-1; ++i)
	{
		for (int j=i+1; j<nN; ++j) 
		{
            int pi = perm[i], pj = perm[j];
			
            double delta = c[pj][pi] - c[pi][pj];

            for (int k = i+1; k<j; k++) delta += c[perm[k]][pi] - c[pi][perm[k]];
			for (int k = i+1; k<j; k++) delta += c[pj][perm[k]] - c[perm[k]][pj];

            if (delta > best_delta)	// record best permutation
			{
				pos1 = i;
				pos2 = j;
				best_delta = delta;
			}
		}
	}

	if (best_delta > 0) // solution improved
	{
		int temp = perm[pos1];
		perm[pos1] = perm[ pos2 ];
		perm[pos2] = temp;
	}
}

bool GraspData::AllocData(int problem_size)
{
    greedy_value = new double [problem_size];
    copy_of_values = new double [problem_size];
    greedy_order = new int [problem_size];
    copy_of_order = new int [problem_size];
    RCL = new int [problem_size];            // memory used by RCL
    length = problem_size;

    costs = new double* [problem_size];
    for (int i=0; i<problem_size; i++)		  
        costs[i] = new double [problem_size];	  

    return true;
}

GraspData::~GraspData()
{
    if (greedy_value) delete[] greedy_value;
    if (copy_of_values) delete[] copy_of_values;
    if (greedy_order) delete[] greedy_order;
    if (copy_of_order) delete[] copy_of_order;
    if (RCL) delete[] RCL;
}

bool GraspData::ReadInstanceValues(ifstream &file)
{
    const int buffer_size = 1024*10;
    char str[buffer_size];
    file >> str;				  // read problem size
    int problem_size = atoi(str);	


    if (!AllocData(problem_size)) return false;

    for (int i=0; i<problem_size; i++) greedy_order[i] = i;

    // read cost matrix
    for (int i=0; i<problem_size; ++i)
    {
        for (int j=0; j<problem_size; ++j)
        {
            file >> str;
            costs[i][j] = atoi(str);
        }
    }
    return true;
}

bool GraspData::InputInstance(const int argc, const char* datafile, int & problem_size)
{
    if (argc != 2) 
    {
        printf("usage: %s <instanceFile>\n");
        return false;
    }

    ifstream file(datafile, ifstream::in);
    if (file.fail()) 
    {
        printf("could not open file %s\n", datafile);
        return false;
    }
    const int buffer_size = 1024*10;
    char str[buffer_size];

    file.getline(str, buffer_size);  // ignore first line 

    printf("GRASP for LOP\n");
    printf("Settings:\n");
    printf("  Max #Iterations: %d\n", nItLimit);
    printf("  Alpha: %lf\n", ALPHA);
    printf("  Beta: %lf\n", beta);
    printf("Problem\n");
    printf("  Filename: %s\n", datafile);
    printf("  Title: %s\n", str);

    ReadInstanceValues(file);

    problem_size = length;

    printf("  Size: %d\n", problem_size);


    // collapse matrix and obtain row sums
    for (int i=0; i<problem_size-1; i++) greedy_value[i] = 0;
    for (int i=0; i<problem_size-1; i++)			
        for (int j=i+1; j<problem_size; j++)
            greedy_value[i] += ( costs[i][j] - costs[j][i] );
    
    greedy_value[problem_size-1] = 0;

    qsort(greedy_order, problem_size, sizeof(int), ::compare);	

    // make a copy of initial values
    for (int i=0; i<problem_size; i++) copy_of_values[i] = greedy_value[i];
    for (int i=0; i<problem_size; i++) copy_of_order[i] = greedy_order[i];
    
    return true;
}

int compare (const void * a, const void * b) //Compare, used to sort.
{
    GraspData &data = g_GRASPData;
	int a1 = *(int*)a;
	int a2 = *(int*)b;
	
    if      (data.GreedyValue(a1) <  data.GreedyValue(a2)) return 1;
	else if (data.GreedyValue(a1) == data.GreedyValue(a2)) return 0;
	else return -1;
}

double GraspData::getSolValue(int *perm, bool isReverse)
{
	double val = 0;
	for (int i=0; i < length-1; i++)
	{
		for (int j=i+1; j < length; j++)
		{
			if (isReverse) val += costs[ perm[j] ][ perm[i] ];
			else           val += costs[ perm[i] ][ perm[j] ];
		}
	}
	return val;
}

