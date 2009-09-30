// solves the linear ordering problem
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INST "c:\\co\\Personal\\Papers\\LinearOrdering\\instances\\random1\\t1d100.1"

#define NewVec(S, N) (S*)malloc(sizeof(S)*(N));
#define DO(I, N) for (int I=0; I<N; ++I)


FILE *openFile(char *name) {
	FILE *f = fopen(name, "r");
	if (!f) {
		printf("error opening file %s\n", name);
		exit(1);
	}
	return f;
}

int allocRow(int **m, int i, int n) {
	m[i] = NewVec(int, n);
	if (!m[i]) return 0;
	memset(m[i], '\0', sizeof(m[0])*n);
	return 1;
}

int **newMatrix(int n) {
	int **m = NewVec(int*, n);
	if (!m) return 0;
	for (int i=0; i<n; ++i) 
		if (!allocRow(m, i, n)) return 0;
	return m;
}

int readInteger(FILE *f) {
	int n;
	fscanf(f, "%d", &n);
	return n;
}

int readMatrixSize(FILE *f) {
	return readInteger(f);
}

int readMatrix(FILE *f, int **m, int n) {
	DO (i, n) {
		DO (j, n) {
			m[i][j] = readInteger(f);
			if (m[i][j] <= 0) {
				printf("error reading entry %d-%d. Value is %d\n", 
					i, j, m[i][j]);
				return 0;
			}
		}
	}
	return 1;
}

void printMatrix(int **m, int n) {
	DO (i, n) {
		DO (j, n) {
			printf("%d ", m[i][j]);
		}
		printf("\n");
	}
}

int old_main(int argv, char **argc) {
	FILE *f = openFile(INST);
	int n = readMatrixSize(f);
	if (!n) {
		printf("wrong number of elements in the matrix: %d\n", n);
		goto end;
	}
	printf("the value is %d\n", n);
	int **m = newMatrix(n);
	readMatrix(f, m, n);
	printMatrix(m, n);
	printf("matrix created\n");

end:
	getc(stdin);
	return 0;
}