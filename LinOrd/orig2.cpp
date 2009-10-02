// (c) 2009, Carlos Oliveira
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
const double ALPHA = 0.3;
const double LARGENUM = 1E9;	

int		Opt_Val;

void LocalSearch(int* pSol, const int, double** c);

// search entire neighborhood
bool LocalSearch2k(int* pSol, const int, double** c, double *delta); 
int compare (const void * a, const void * b);

class GraspData {
	double *greedy_value;
	double *copy_of_values;			//Copy of Initial Greedy function
	int *greedy_order;
	int *copy_of_order;			//Copy of Initial greedy order
	int length;
    double **costs;

public:

    void  ConstructSolution(int* pSol);
    bool InputInstance(const int argc, const char* datafile, int& problem_size);
    double GetSolutionValue (int* permutation);
    void printv();
    bool AllocData(int problem_size);
    bool ReadInstanceValues(ifstream &file);
    double GreedyValue(int pos) { return this->greedy_value[pos]; }
    double **GetCosts() { return costs; }

    GraspData() 
        : greedy_value(0)
        , copy_of_values(0)
        , greedy_order(0)
        , copy_of_order(0)
        , length(0)
        , costs(0) {}
    ~GraspData();

};

GraspData g_GRASPData;


int main(int argc, char** argv)
{
	srand (1); // (unsigned)time(NULL));

    int problem_size = 0;
    const char *fileName = argv[1];
    if (g_GRASPData.InputInstance(argc, fileName, problem_size) == false ) 
    {
        printf("error reading instance %s\n", fileName);
        return 1;
    }

	int *currentSolution = new int [problem_size];
	int *bestSolution = new int [problem_size];

    clock_t total_const_time=0, total_ls_time=0;
    int num_ls_iter = 0;
    double bestValue = 0;
    printf("constructor  (time)   LS   #iter  (time)  \n");
	for (int num_iterations=1; num_iterations < nItLimit; ++num_iterations)
	{
        clock_t ct1 = clock();
		g_GRASPData.ConstructSolution(currentSolution);
		clock_t ct2 = clock();
        total_const_time += ct2-ct1;

        double delta, objective = g_GRASPData.GetSolutionValue(currentSolution);
        printf("%.0lf %ld ", objective, ct2-ct1); fflush(stdout);

        clock_t lst1 = clock();
        int n=0;
        for (; LocalSearch2k(currentSolution, problem_size, g_GRASPData.GetCosts(), &delta); ++n)
            objective += delta;
        clock_t lst2 = clock();
        total_ls_time += lst2 - lst1;
        num_ls_iter += n;
        printf("%.0lf %.2lf %d \n", objective, (double)(lst2-lst1)/n, n); fflush(stdout);

		if (objective > bestValue) {
			for (int i=0; i<problem_size; i++) bestSolution[i] = currentSolution[i];
			bestValue = objective;
		}

		if (num_iterations % 10 == 0)
		{
            if (num_iterations == 10) assert(objective == 3391268);
            if (num_iterations == 10) assert(bestValue == 3391268);
            printf("iter: %d best %d constr_time: %.2lf ls_time: %.2lf\n",
                num_iterations, (int)bestValue,
                (double)total_const_time/num_iterations,
                (double)total_ls_time/num_ls_iter);
		}
	}

	cout << "\nFinal Permutation:\n";
	for (int i=0; i<problem_size; i++)
	{
        printf("%d %s", bestSolution[i], (((i+1) % 10 == 0) ? "\n" : ""));
	}

    delete [] currentSolution;
    delete [] bestSolution;
	return 0;
}

// tests to do:
// 1. define the RCL size based on 
//    a. values (as done currently)
//    b. percentage of maximum size

// other TODOs:
// * check local search
// * in the constructor, check how to find the location where to add the new element
// * in the constructor, implement the local search performed after each item is selected

void GraspData::ConstructSolution(int* perm)
{
	double limit;

	for (int k=0, num_remaining=length; k<length; ++k, --num_remaining)
	{
        double min_value = greedy_value[greedy_order[num_remaining-1]];
		limit = (greedy_value[greedy_order[0]] - min_value) * (1 - ALPHA) + min_value;
		
        int rcl_size = 1;
        for (int i=1; i<num_remaining && greedy_value[greedy_order[i]] >= limit; i++)  rcl_size++;
		
        int elem = perm[k] = greedy_order[rand() % rcl_size];  // select an element from RCL

		greedy_value[elem] = -LARGENUM;
		
        for (int i=0; i<length; ++i) // update attractiveness
			if (greedy_value[i] != -LARGENUM)
				greedy_value[i] -= ( costs[i][elem] - costs[elem][i] );

        qsort(greedy_order, num_remaining, sizeof(int), ::compare);	// update greedy order
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
bool LocalSearch2k(int* perm, const int nN, double** c, double *delta)
{
	int pos1, pos2;
	double &best_delta = *delta;
    best_delta = 0;
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
        return true;
	}
    return false;
}

bool GraspData::AllocData(int problem_size)
{
    greedy_value = new double [problem_size];
    copy_of_values = new double [problem_size];
    greedy_order = new int [problem_size];
    copy_of_order = new int [problem_size];
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


    // calculate attractiveness factors for each position
    for (int i=0; i<problem_size; i++) greedy_value[i] = 0;
    for (int i=0; i<problem_size; i++)			
        for (int j=0; j<problem_size; j++)
            greedy_value[i] += ( costs[i][j] - costs[j][i] );
    

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

double GraspData::GetSolutionValue(int *perm)
{
	double val = 0;
	for (int i=0; i < length-1; ++i)
		for (int j=i+1; j < length; ++j) val += costs[perm[i]][perm[j]];

    return val;
}

