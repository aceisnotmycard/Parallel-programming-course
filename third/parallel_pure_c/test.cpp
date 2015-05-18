#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

#define in 20
#define jn 20
#define kn 20
#define a 1

int rank, size;
int layer_height;
double* F;
double hx, hy, hz;

#define F(iter, x, y, z) F[iter * ((layer_height+2)*(jn+1)*(kn+1)) + x*(jn+1)*(kn+1) + y*(kn+1) + z]

double Fresh (double, double, double);
double Ro (double, double, double);
void Inic();
int get_chunk_size(int given_rank);

double Fresh (double x, double y, double z) {
	return x + y + z;
}

double Ro (double x, double y, double z) {
	return -a * (x + y + z);
}

void Inic() {
	int starting_row = (rank == 0) ? 1 : 0;
	int ending_row = (rank == size - 1) ? layer_height :  (layer_height + 1);

	int offset = 0;
	for (int i = 0; i < rank; i++) {
		offset += get_chunk_size(i);
	}

	int i, j, k;
	for (i = starting_row; i <= ending_row; i++) {
			for (j = 0; j <= jn; j++) {
				for (k = 0; k <= kn; k++) {
					if ((i != starting_row) && (j != 0) && (k != 0) && (i != ending_row) && (j != jn) && (k != kn)) {
						F(0, i, j, k) = 0;
						F(1, i, j, k) = 0;
					} else {
						F(0, i, j, k) = Fresh((offset + i) * hx, j * hy, k * hz);
						F(1, i, j, k) = Fresh((offset + i) * hx, j * hy, k * hz);
					}
				}
			}
	}

}

int get_chunk_size(int given_rank) {
	int basic_chunk = kn / size;
	int rest = kn % size;
	return basic_chunk + (given_rank < rest ? 1 : 0);
}

int main(int argc, char **argv) {
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double X, Y, Z;
	double max, N, t1, t2;
	double owx, owy, owz, c, e;
	double Fi, Fj, Fk, F1;
	double global_max;
	int i, j, k, mi, mj, mk;
	int R, fl, fcur_it, fl2;
	int it, f, prev_it, cur_it;
	double osdt;
	struct timeval tv1, tv2;
	int global_f;
	layer_height = get_chunk_size(rank);
	F = (double*) malloc(2 * (layer_height + 2) * (jn + 1) * (kn + 1) * sizeof(double));
	int starting_row = (rank == 0) ? 2 : 1;
	int ending_row = (rank == size - 1) ? layer_height : (layer_height + 1);
	MPI_Request reqs1, reqr2, reqs2, reqr1;

	it = 0;
	X = 2.0;
	Y = 2.0;
	Z = 2.0;
	e = 0.00001;
	prev_it = 1;
	cur_it = 0;
	hx = X / in;
	hy = Y / jn;
	hz = Z / kn;
	owx = pow(hx, 2);
	owy = pow(hy, 2);
	owz = pow(hz, 2);
	c = 2 / owx + 2 / owy + 2 / owz + a;

	gettimeofday(&tv1, (struct timezone*) 0);
	Inic();
	do {
		MPI_Irecv(&F(prev_it,0,0,0), jn*kn, MPI_DOUBLE, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, &reqr2);
		MPI_Isend(&F(prev_it,layer_height,0,0), jn*kn, MPI_DOUBLE, (rank + 1) % size, 1, MPI_COMM_WORLD, &reqs1);

		MPI_Isend(&F(prev_it,1,0,0), jn*kn, MPI_DOUBLE, (rank - 1 + size) % size, 2, MPI_COMM_WORLD, &reqr1);
		MPI_Irecv(&F(prev_it,(layer_height + 1),0,0), jn*kn, MPI_DOUBLE, (rank + 1) % size, 2, MPI_COMM_WORLD, &reqs2);

		f = 1;
		prev_it = 1 - prev_it;
		cur_it = 1 - cur_it;

		for (i = starting_row; i < ending_row; i++) {
			for (j = 1; j < jn; j++) {
				for (k = 1; k < kn; k++) {
					Fi = (F(prev_it, (i+1), j, k) + F(prev_it, (i-1), j, k)) / owx;
					Fj = (F(prev_it, i, (j+1), k) + F(prev_it, i, (j-1), k)) / owy;
					Fk = (F(prev_it, i, j, (k+1)) + F(prev_it, i, j, (k-1))) / owz;
					F(cur_it, i, j, k) = (Fi + Fj + Fk - Ro(i * hx, j * hy, k * hz)) / c;
					if (fabs(F(cur_it, i, j, k) - F(prev_it, i, j, k)) > e) {
						f = 0;
					}
				}
			}
		}
		
		if (rank != size - 1) {
			MPI_Wait(&reqr1, MPI_STATUS_IGNORE);
			for (j = 1; j < jn; j++) {
				for (k = 1; k < kn; k++) {
					Fi = (F(prev_it, (layer_height+1), j, k) + F(prev_it, (layer_height-1), j, k)) / owx;
					Fj = (F(prev_it, layer_height, (j+1), k) + F(prev_it, layer_height, (j-1), k)) / owy;
					Fk = (F(prev_it, layer_height, j, (k+1)) + F(prev_it, layer_height, j, (k-1))) / owz;
					F(cur_it, layer_height, j, k) = (Fi + Fj + Fk - Ro(layer_height * hx, j * hy, k * hz)) / c;
					if (fabs(F(cur_it, i, j, k) - F(prev_it, i, j, k)) > e) {
						f = 0;
					}
				}
			}
		}	

		if (rank != 0) {
			MPI_Wait(&reqr2, MPI_STATUS_IGNORE);

			for (j = 1; j < jn; j++) {
				for (k = 1; k < kn; k++) {
					Fi = (F(prev_it, 2, j, k) + F(prev_it, 0, j, k)) / owx;
					Fj = (F(prev_it, 1, (j+1), k) + F(prev_it, 1, (j-1), k)) / owy;
					Fk = (F(prev_it, 1, j, (k+1)) + F(prev_it, 1, j, (k-1))) / owz;
					F(cur_it, 1, j, k) = (Fi + Fj + Fk - Ro(1 * hx, j * hy, k * hz)) / c;
					if (fabs(F(cur_it, i, j, k) - F(prev_it, i, j, k)) > e) {
						f = 0;
					}
				}
			}
		}
		MPI_Allreduce(&f, &global_f, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		it++;
	} while (global_f == 0);

	gettimeofday(&tv2, (struct timezone*) 0);
	osdt = (tv2.tv_sec - tv1.tv_sec) + (double) (tv2.tv_usec - tv1.tv_usec) / 1000000;

	printf("\n rank = %d in = %d iter = %d E = %f T = %fs \n", rank, in, it, e, osdt);

	max = 0.0;
	for (i = starting_row; i < ending_row; i++) {
		for (j = 1; j < jn; j++) {
			for (k = 1; k < kn; k++) {
				if ((F1 = fabs(F(cur_it, i, j, k) - Fresh(i * hx, j * hy, k * hz))) > max) {
					max = F1;
					mi = rank * layer_height + i;
					mj = j;
					mk = k;
				}
			}
		}
	}
	MPI_Allreduce(&max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	printf("\n rank = %d, Max differ = %f in point (%d, %d, %d) \n\n", rank, global_max, mi, mj, mk);

	MPI_Finalize();
}