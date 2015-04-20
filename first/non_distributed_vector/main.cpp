#include <mpi.h>
#include <cstdlib>
#include <cstdio>

#define EPSILON 0.00001
#define POSITIVE_TAU  0.00001
#define NEGATIVE_TAU -0.00001
//
#define N 4096

int size, rank;

// creates array that represents part of array
double* create_matrix();

//creates vector
double* create_vector();

// calculates number of rows that will belong 
// to process of given rank
int get_chunk_size(int given_rank);

// calculates vertical position of the first element
int calculate_shift(int given_rank);

// name tells for itself
void print_matrix(double* matrix);
void print_vector(double* vector);

double* multiply_matrix_by_vector(double* matrix, double* vector);
double* subtract_vector_from_vector(double* a, double* b);
double* multiply_scalar_by_vector(double scalar, double* vector);
double squared_norm_of_vector(double* vector);

double* calculate(double* matrix, double* vector);
double measure(double* A, double* x, double* b);

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	double* matrix = create_matrix();
	double* vector = create_vector();

	//print_matrix(matrix);
	double start_time = MPI_Wtime();
	double *result = calculate(matrix, vector);
	double end_time = MPI_Wtime();

	//print_vector(result);
	printf("PROCESS: %d TIME: %f\n", rank, end_time - start_time);

	free(matrix);
	free(vector);
	free(result);

	MPI_Finalize();
	return 0;
}

int get_chunk_size(int given_rank) {
	int basic_chunk = N / size;
	int rest = N % size;
	return basic_chunk + (given_rank < rest ? 1 : 0);
}

int calculate_shift(int given_rank) {
	int result = 0;
	for (int i = 0; i < given_rank; i++) {
		result += get_chunk_size(i);
	}
	return result;
}

double* create_matrix() {
	int chunk_size = get_chunk_size(rank);
	double* chunk = (double*) calloc(chunk_size * N, sizeof(double));
	int shift = calculate_shift(rank);
	for (int i = 0; i < chunk_size; i++) {
		for (int j = 0; j < N; j++) {
			chunk[i * N + j] = (i == (j - shift)) ? 2 : 1;
		}
	}
	return chunk;
}

double* create_vector() {
	double* vector = (double*) calloc(N, sizeof(double));
	for (int i = 0; i < N; i++) {
		vector[i] = N + 1;
	}
	return vector;
}

void print_matrix(double* matrix) {
	int chunk_size = get_chunk_size(rank);
	for (int process = 0; process < size; process++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == process) {
			for (int i = 0; i < chunk_size; i++) {
				for (int j = 0; j < N; j++) {
					printf("%0.0f ", matrix[i * N + j]);
				}
				printf("\n");
			}
		}
	}
}

void print_vector(double* vector) {
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		for (int i = 0; i < N; i++) {
			printf("%0.4f\n", vector[i]);
		}
	}
}

double* multiply_matrix_by_vector(double* matrix, double* vector) {
	int chunk_size = get_chunk_size(rank);
	double* result = (double*) calloc(chunk_size, sizeof(double));
	for (int i = 0; i < chunk_size; i++) {
		for (int j = 0; j < N; j++) {
			result[i] += matrix[N * i + j] * vector[j];
		}
	}
	return result;
}

double* subtract_vector_from_vector(double* a, double* b) {
	double* result = (double*) calloc(N, sizeof(double));
	for (int i = 0; i < N; i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}

double* multiply_scalar_by_vector(double scalar, double* vector) {
	double* result = (double*) calloc(N, sizeof(double));
	for (int i = 0; i < N; i++) {
		result[i] = scalar * vector[i];
	}
	return result;
}

double squared_norm_of_vector(double* vector) {
	double sum = 0;
	for (int i = 0; i < N; i++) {
		sum += vector[i] * vector[i];
	}
	return sum;
}

// sizes
int* count_recvcounts() {
	int* counts = (int*) calloc(size, sizeof(int));
	for (int i = 0; i < size; i++) {
		counts[i] = get_chunk_size(i);
	}
	return counts;
}

// positions
int* count_displs() {
	int* counts = (int*) calloc(size, sizeof(int));
	for(int i = 0; i < size; i++) {
		counts[i] = calculate_shift(i);
	}
	return counts;
}

double* calculate(double* matrix, double* vector) {
	double* result = (double*) calloc(N, sizeof(double));

	int* c_displs = count_displs();
	int* c_recvcounts = count_recvcounts();
	while (1) {
		// Ax
		double* u = multiply_matrix_by_vector(matrix, result);
		// Так как мы умножаем строку на весь вектор в каждом процессе, 
		// то нужно получить весь вектор в каждом процессе
		double* z = (double*) calloc(N, sizeof(double));
		MPI_Allgatherv(u, c_recvcounts[rank], MPI_DOUBLE, z, c_recvcounts, c_displs, MPI_DOUBLE, MPI_COMM_WORLD);

		// Ax - b
		double* y = subtract_vector_from_vector(z, vector);
		//print_vector(result);
		if (squared_norm_of_vector(y)/squared_norm_of_vector(vector) < EPSILON*EPSILON) {
			free(u);
			free(y);
			free(z);
			break;
		}

		// t(Ax - b)
		double* v = multiply_scalar_by_vector(POSITIVE_TAU, y);
		// x - t(Ax - b)
		result = subtract_vector_from_vector(result, v);
	}
	return result;
}
