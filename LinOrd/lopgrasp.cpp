// (c) 2009, Carlos Oliveira
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <assert.h>
#include <string>
#include <vector>

using namespace std;

const int    nItLimit = 5000;  // iteration limit
const double ALPHA = 0.3;
const double LARGENUM = 1E9;	
const int ELITE_SIZE = 20;

struct GraspData {
    double *greedy_value;
    double *copy_of_values;
    int *greedy_order;
    int *copy_of_order;
    int length;
    double **costs;
    double **ls_entries; // used by local search
};

void  ConstructSolution(GraspData *d, int* pSol);
bool InputInstance(GraspData *d, const int argc, char **argv, int& problem_size);
double GetSolutionValue (GraspData *d, int* permutation);
void printv(GraspData *d);
bool AllocData(GraspData *d, int problem_size);
bool ReadInstanceValues(GraspData *d, FILE *file);
double GreedyValue(GraspData *d, int pos) { return d->greedy_value[pos]; }
double **GetCosts(GraspData *d) { return d->costs; }
double **GetLSEntries(GraspData *d) { return d->ls_entries; }


GraspData g_GRASPData;

void LocalSearchForLastEntry(int* pSol, const int, double** c);
bool LocalSearch(int* perm, const int n, double** c, double *delta, double **M);
bool LocalSearch(int* perm, const int n, double *delta, GraspData *data);
int compare (const void * a, const void * b); // used by qsort algorithm


struct EliteEntry {
    int size;
    int *perm;
    double cost;
    EliteEntry(int n) : size(n), perm(new int[n]), cost(-LARGENUM) {}
    ~EliteEntry() { delete[] perm; }
    void SetSol(int *p, double c) { for (int i=0; i<size; ++i) perm[i] = p[i]; cost = c;}
    int *GetPerm() { return perm; }
    bool Equal(int *p, double c) {
        if (c!=cost) return false;
        for (int i=0; i<size; ++i) if (perm[i]!=p[i]) return false;
        return true;
    }
};

int instanceError(const char *f) { printf("error reading instance %s\n", f); return 1; }
int usage(const char *p) { printf("usage: %s <instanceFile>\n", p); return 0; }
int failFile(const char *f) { printf("could not open file %s\n", f); return 0; }
int failMem() { printf("error allocating memory\n"); return 1; }

#ifdef NDEBUG
#define DODBG(x) 
#else
#define DODBG(x) x
#endif

void printfRes0(double perc, double objective, double t) {
    DODBG(printf("%.2lf%% %.0lf %ld ", perc, objective, t); fflush(stdout);)
}

void printfRes1(double o, double t1, int n, double perc) {
    DODBG(printf("%.0lf %.2lf %d (imp %.2lf) ", o, t1, n, perc);)
}

void printfRes2(double o, double t2, double perc) {
    DODBG(printf("%.0lf %.2lf (imp %.2lf)", o, t2, perc);)
}

void checkRes(int num_iterations, double objective, double bestValue) {
DODBG(
    if (num_iterations == 10) assert(objective == 3391268);
    if (num_iterations == 10) assert(bestValue == 3391268);
)
}

template<class T>bool allocVec(T *&v, int n) { v = new T[n]; return v != 0; }
#define ALLOC_OR_RET(x, n) if (!allocVec(x, n)) return false;

bool allocMemory(int n, int *&current_solution, int *&best_solution, int *&elite_start,
                 EliteEntry ** &elite_set)
{
    ALLOC_OR_RET(current_solution, n);
    ALLOC_OR_RET(best_solution, n);
    ALLOC_OR_RET(elite_start, n);
    ALLOC_OR_RET(elite_set, ELITE_SIZE);

    for (int i=0; i<ELITE_SIZE; ++i) 
        if (!(elite_set[i] = new EliteEntry(n))) return false;
    return true;
}

// store temporary results, so we don't spend time printing them
struct Results {
    double obj1;
    long time1;
    double obj2;
    double time2;
    int n;
    double perc1;
    double obj3;
    long time3;
    double perc2;
};


void printRes(const vector<Results> &resultList)
{
#ifdef NDEBUG
    FILE *f = fopen("..\\output.txt", "w");
    int it_num=0;
    for (vector<Results>::const_iterator i=resultList.begin(); i!=resultList.end(); i++, it_num++)
    {
        fprintf(f, "%.0lf %ld ", i->obj1, i->time1); fflush(stdout);
        fprintf(f, "%.0lf %.2lf %d (imp %.2lf) ", i->obj2, i->time2, i->n, i->perc1); fflush(stdout);
        if (it_num > ELITE_SIZE)
        {
            if (i->obj3 > i->obj2)
                fprintf(f, "%.0lf %.2lf (imp %.2lf)", i->obj3, (double)i->time3, i->perc2);
        }
        fprintf(f, "\n");
    }
    fclose(f);
#endif
}

void insertSol(EliteEntry **e, int n, int pos, double cost, int *sol) {
    EliteEntry *t = e[n-1];    // reuse last pos    
    for (int i=n-1; i>pos; --i) e[i] = e[i-1]; // make room
    (e[pos] = t)->SetSol(sol, cost); // store at empty pos
}

void insertIntoElite(EliteEntry **e, double cost, int *sol) {
    const int n = ELITE_SIZE;

    // find the first position where a solution can be inserted
    int pos=0;  for (; pos<n && e[pos]->cost > cost;) ++pos; 
    if (pos<n && !e[pos]->Equal(sol, cost)) insertSol(e, n, pos, cost, sol);
}


int main(int argc, char** argv)
{
   if (argc != 2) return usage(argv[0]);

    srand (1); // (unsigned)time(NULL));
    clock_t startTime = clock();
    int problem_size = 0;
    const char *fileName = argv[1];
    
    // input instance
    if (!InputInstance(&g_GRASPData, argc, argv, problem_size)) 
        return instanceError(fileName);

    // alloc memory
    int *current_solution, *best_solution,*elite_start;
    EliteEntry **elite_set;
    if (!allocMemory(problem_size, current_solution, best_solution, elite_start, elite_set)) 
        return failMem();

    vector<Results> resultList;
    clock_t total_const_time=0, total_ls_time=0;
    int num_ls_iter = 0;
    double bestValue = 0;
    printf("constructor  (time)   LS   #iter  (time)  \n");
    for (int num_iterations=1; num_iterations < nItLimit; ++num_iterations)
    {
        // apply constructor
        clock_t ct1 = clock();
        ConstructSolution(&g_GRASPData, current_solution);
        clock_t ct2 = clock();
        total_const_time += ct2-ct1;

        // get cost, store values
        double delta, objective = GetSolutionValue(&g_GRASPData, current_solution);
        Results res;
        res.obj1 = objective; res.time1 = ct2-ct1;
        double percCompleted = 100*(num_iterations/(double)nItLimit);
        printfRes0(percCompleted, objective, ct2-ct1);

        // appy local search
        clock_t lst1 = clock();
        int n=0;
        double orig_objective = objective;
        for (; LocalSearch(current_solution, problem_size, &delta, &g_GRASPData); ++n)
            objective += delta;
        clock_t lst2 = clock();
        total_ls_time += lst2 - lst1;
        num_ls_iter += n;

        // store values
        res.obj2 = objective; res.time2 = (double)(lst2-lst1)/n; res.n = n;
        res.perc1 = 100*(objective-orig_objective)/orig_objective;
        printfRes1(res.obj2, res.time2, n, res.perc1);

        // apply path relinking
        if (num_iterations > ELITE_SIZE) { 
            clock_t prt1 = clock();
            orig_objective = objective;
            int rand_pos = rand() % ELITE_SIZE;
            int *guiding_solution = elite_set[rand_pos]->GetPerm();
            for (int i=0; i<problem_size; ++i) elite_start[i] = current_solution[i];
            for (int i=0; i<problem_size; ++i) {
                int elem = guiding_solution[i];
                int elem_pos = 0;
                for (int j=0; j<problem_size; ++j) if (elite_start[j]==elem) { elem_pos=j; break; }
                // change elite_start to the contain the element on position i
                elite_start[elem_pos] = elite_start[i];
                elite_start[i] = elem;
                // TODO: add local search here
            }
            clock_t prt2 = clock();
            objective = GetSolutionValue(&g_GRASPData, elite_start);

            res.obj3 = objective;
            if (objective > orig_objective) {
                for (int i=0; i<problem_size; ++i) current_solution[i] = elite_start[i];
                res.time3 = prt2-prt1; res.perc2 = 100*(objective-orig_objective)/orig_objective;
                printfRes2(objective, res.time3, res.perc2);
            } else {
                objective = orig_objective;
            }
        }

        DODBG(printf("\n"));

        insertIntoElite(elite_set, objective, current_solution);

        // save solution
        if (objective > bestValue) {
            for (int i=0; i<problem_size; i++) best_solution[i] = current_solution[i];
            bestValue = objective;

            clock_t currentTime = clock();
            printf("iter: %d best %d constr_time: %.2lf ls_time: %.2lf %ld\n",
                num_iterations, (int)bestValue,
                (double)total_const_time/num_iterations,
                (double)total_ls_time/num_ls_iter, currentTime - startTime);
        }

        checkRes(num_iterations, objective, bestValue);
        resultList.push_back(res);
    }

    // print results
    printRes(resultList);
    printf("\nFinal Permutation:\n");
    for (int i=0; i<problem_size; i++)
    {
        printf("%d %s", best_solution[i], (((i+1) % 10 == 0) ? "\n" : ""));
    }

    for (int i=0; i<ELITE_SIZE; ++i) delete elite_set[i];
    delete[] elite_set;

    delete [] elite_start;
    delete [] current_solution;
    delete [] best_solution;
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

void ConstructSolution(GraspData *d, int* perm)
{
    double limit;
    int length = d->length;
    double *greedy_value = d->greedy_value;
    int *greedy_order = d->greedy_order;
    double **costs = d->costs;

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

    for (int i=0; i < length; i++) greedy_value[i] = d->copy_of_values[i];
    for (int i=0; i < length; i++) greedy_order[i] = d->copy_of_order[i];
}

// LocalSearch. Checks only combinations containing the last value.
void LocalSearchForLastEntry(int* perm, const int nN, double** c)
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

bool LocalSearch(int* perm, const int n, double *delta, GraspData *data)
{
    return LocalSearch(perm, n, GetCosts(data), delta, GetLSEntries(data));
}

// perform search on full neighborhood	
bool LocalSearch(int* perm, const int n, double** c, double *delta, double **M)
{
    int pos1, pos2;
    double &best_delta = *delta;
    best_delta = 0;

    for (int k=0; k<n-1; ++k)
    {
        for (int i=0, j=k+1; j<n; ++i, ++j)
        {
            int pi = perm[i], pj = perm[j];
            int ppi = perm[i+1], pmj = perm[j-1];

            double deltaij = c[pj][pi] - c[pi][pj];
            if (j-i > 1) deltaij += M[i][j-1] + M[i+1][j];
            if (j-i > 2) deltaij += - M[i+1][j-1] - (c[pmj][ppi] - c[ppi][pmj]);
            M[i][j] = deltaij;

            if (deltaij > best_delta)  pos1 = i, pos2 = j, best_delta = deltaij;
        }
    }

    static int test = 0;
    if (test == 0) test = 1, assert(best_delta == 35774);
    if (best_delta > 0) // solution improved
    {
        int temp = perm[pos1];
        perm[pos1] = perm[ pos2 ];
        perm[pos2] = temp;
        return true;
    }
    return false;
}

bool AllocData(GraspData *d, int n)
{
    d->length = n;
    ALLOC_OR_RET(d->greedy_value, n);
    ALLOC_OR_RET(d->copy_of_values, n);
    ALLOC_OR_RET(d->greedy_order, n);
    ALLOC_OR_RET(d->copy_of_order, n);

    ALLOC_OR_RET(d->costs, n);
    for (int i=0; i<n; i++) ALLOC_OR_RET(d->costs[i], n);

    ALLOC_OR_RET(d->ls_entries, n);
    for (int i=0; i<n; ++i) ALLOC_OR_RET(d->ls_entries[i], n);

    return true;
}

template<class T>void vecDelete(T *v) { if (v) delete[] v; }
template<class T>void mvecDelete(T *v, int n) { 
    if (v) for (int i=0; i<n; ++i) vecDelete(v[i]);
    vecDelete(v);
}

void deleteGrasp(GraspData *d)
{
    vecDelete(d->greedy_value);
    vecDelete(d->copy_of_values);
    vecDelete(d->greedy_order);
    vecDelete(d->copy_of_order);
    mvecDelete (d->costs, d->length);
    mvecDelete (d->ls_entries, d->length);
}

void readStr(FILE *f, char *s) { fscanf(f, "%s", s); }
void readEntry(FILE *f, char *s, double **c, int i, int j) { readStr(f, s); c[i][j] = atoi(s); }

bool ReadInstanceValues(GraspData *d, FILE *file)
{
    const int buffer_size = 1024*10;
    char str[buffer_size];
    readStr(file, str);
    int n = atoi(str); // read problem size

    if (!AllocData(d, n)) return false;

    for (int i=0; i<n; i++) d->greedy_order[i] = i;

    // read cost matrix
    for (int i=0; i<n; ++i) 
        for (int j=0; j<n; ++j) readEntry(file, str, d->costs, i, j);

    return true;
}

bool InputInstance(GraspData *d, const int argc, char **argv, int& problem_size)
{
    FILE *file = fopen(argv[1], "r");
    if (!file) return failFile(argv[1]);
    
    const int buffer_size = 1024*10;
    char str[buffer_size];

    fgets(str, buffer_size, file);  // ignore first line 

    printf("GRASP for LOP\n");
    printf("Settings:\n");
    printf("  Max #Iterations: %d\n", nItLimit);
    printf("  Alpha: %lf\n", ALPHA);
    printf("Problem\n");
    printf("  Filename: %s\n", argv[1]);
    printf("  Title: %s\n", str);

    ReadInstanceValues(d, file);
    fclose(file);

    problem_size = d->length;
    double *greedy_value = d->greedy_value;

    printf("  Size: %d\n", problem_size);

    // calculate attractiveness factors for each position
    for (int i=0; i<problem_size; i++) greedy_value[i] = 0;
    for (int i=0; i<problem_size; i++)			
        for (int j=0; j<problem_size; j++)
            greedy_value[i] += ( d->costs[i][j] - d->costs[j][i] );

    

    qsort(d->greedy_order, problem_size, sizeof(int), ::compare);	

    // make a copy of initial values
    for (int i=0; i<problem_size; i++) d->copy_of_values[i] = greedy_value[i];
    for (int i=0; i<problem_size; i++) d->copy_of_order[i] = d->greedy_order[i];
    
    return true;
}

int compare (const void * a, const void * b) //Compare, used to sort.
{
    GraspData &data = g_GRASPData;
    int a1 = *(int*)a;
    int a2 = *(int*)b;
    
    if      (GreedyValue(&data, a1) <  GreedyValue(&data, a2)) return 1;
    else if (GreedyValue(&data, a1) == GreedyValue(&data, a2)) return 0;
    else return -1;
}

double GetSolutionValue(GraspData *d, int *perm)
{
    double val = 0;
    for (int i=0; i < d->length-1; ++i)
        for (int j=i+1; j < d->length; ++j) val += d->costs[perm[i]][perm[j]];

    return val;
}

