#include <mpi.h>
#include <cstdio>
#include <cstdlib>

#define EPSILON 0.00001
#define POSITIVE_TAU  0.01
#define NEGATIVE_TAU -0.01

#define TAG_SEND_MATRIX 100
#define TAG_SEND_OFFSET 101
#define TAG_SEND_LENGTH 102
//
#define N 199


int size, rank;

// chunks sizes
int* c_recvcounts;
// chunks positions
int* c_displs;

int get_chunk_size(int given_rank);
int calculate_shift(int given_rank);

double* create_matrix();
double* create_vector();

double* subtract_vector_from_vector(double* a, double* b);
double* multiply_scalar_by_vector(double scalar, double* vector);
double* multiply_matrix_by_vector(double* matrix, double* vector);
double squared_norm_of_vector(double* vector);

double* calculate(double* matrix, double* vector);

int* count_recvcounts();
int* count_displs();

void print_vector(double* vector);
void print_matrix(double* matrix);

int main(int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	c_displs = count_displs();
	c_recvcounts = count_recvcounts();

	double* matrix = create_matrix();
	double* vector = create_vector();
	//print_vector(vector);

	double start_time = MPI_Wtime();
	double *result = calculate(matrix, vector);
	double end_time = MPI_Wtime();

	print_vector(result);
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
	int chunk_size = c_recvcounts[rank];
	double* chunk = (double*) calloc(chunk_size * N, sizeof(double));
	int shift = c_displs[rank];
	for (int i = 0; i < chunk_size; i++) {
		for (int j = 0; j < N; j++) {
			chunk[i * N + j] = (i == (j - shift)) ? 2 : 1;
		}
	}
	return chunk;
}

double* create_vector() {
	int chunk_size = c_recvcounts[rank];
	double* vector = (double*) calloc(chunk_size, sizeof(double));
	for (int i = 0; i < chunk_size; i++) {
		vector[i] = N + 1;
	}
	return vector;
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

void print_vector(double* vector) {
	for (int process = 0; process < size; process++) {
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank == process) {
			for (int i = 0; i < c_recvcounts[rank]; i++) {
				printf("%0.4f\n", vector[i]);
			}
		}
	}
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

double squared_norm_of_vector(double* vector) {
	double sum;
	double chunk_sum = 0;
	for (int i = 0; i < c_recvcounts[rank]; i++) {
		chunk_sum += vector[i] * vector[i];
	}
	MPI_Allreduce(&chunk_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	return sum;
}

double* subtract_vector_from_vector(double* a, double* b) {
	double* result = (double*) calloc(c_recvcounts[rank], sizeof(double));
	for (int i = 0; i < c_recvcounts[rank]; i++) {
		result[i] = a[i] - b[i];
	}
	return result;
}

double* multiply_scalar_by_vector(double scalar, double* vector) {
	double* result = (double*) calloc(c_recvcounts[rank], sizeof(double));
	for (int i = 0; i < c_recvcounts[rank]; i++) {
		result[i] = scalar * vector[i];
	}
	return result;
}

double* multiply_matrix_by_vector(double* matrix, double* vector) {
	int vector_length = c_recvcounts[rank];
	int offset = c_displs[rank];
	int incoming_process_data = 0;

	double* result = (double*) calloc(vector_length, sizeof(double));
	double* v = (double*) calloc(N/size + 1, sizeof(double));
	for (int i = 0; i < vector_length; i++) {
		v[i] = vector[i];
	}
	for (int process = 0; process < size; process++) {
		// index of current part of vector
		incoming_process_data = (rank + process) % size;

		for (int i = 0; i < c_recvcounts[rank]; i++) {
			for (int j = 0; j < c_recvcounts[incoming_process_data]; j++) {
				result[i] += matrix[i * N + j + c_displs[incoming_process_data]] * v[j];
			}
		}
		// switch vector, (vector length, vector offset) between processes
		MPI_Sendrecv_replace(v, N/size + 1, MPI_DOUBLE, (rank+1) % size, TAG_SEND_MATRIX,
			(rank-1) % size, TAG_SEND_MATRIX, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	return result;
}

double* calculate(double* matrix, double* vector) {
	double* result = (double*) calloc(N, sizeof(double));
	// Just a little optimization
	double vector_squared_norm = squared_norm_of_vector(vector);
	while (1) {
		// Ax
		double* z = multiply_matrix_by_vector(matrix, result);
		// Ax - b
		double* y = subtract_vector_from_vector(z, vector);
		double norm = squared_norm_of_vector(y);
		if (norm/vector_squared_norm < EPSILON*EPSILON) {
			free(y);
			free(z);
			break;
		}

		// t(Ax - b)
		double* v = multiply_scalar_by_vector(POSITIVE_TAU, y);
		// x - t(Ax - b)
		result = subtract_vector_from_vector(result, v);
		free(v);
	}
	return result;
}