#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TASKS_LIST_SIZE 20
#define SIZE_MULTIPLIER 555555

// for MPI
int rank, size;
MPI_Comm comm;
// logical variable
// 1 if worker is done and waiting for more
// 0 if it is busy
int worker_is_ready;
// array with sizes of tasks
int* tasks_list;

// prepares data in 0 process
// TODO: remove argument
void prepare_tasks(int rank);

void* worker();

void listener();

int main(int argc, char** argv) {
	int provided;

	// prepare MPI
	comm = MPI_COMM_WORLD;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	//prepare data
	tasks_list = (int*) malloc(TASKS_LIST_SIZE * sizeof(int));
	prepare_tasks(rank);

	//prepare working thread
	pthread_t working_thread;
	pthread_attr_t attrs;

	if (pthread_attr_init(&attrs) != 0) {
		perror("Cannot init thread attrs");
		exit(1);
	}
	if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE) != 0) {
		perror("Cannot set thread joinable");
		exit(1);
	}
	pthread_create(&working_thread, &attrs, worker, NULL);
	pthread_attr_destroy(&attrs);

	// start listening for messages from another processes
	listener();

	// ok, going home
	pthread_join(working_thread, NULL);
	MPI_Finalize();
	return 0;
}

void prepare_tasks(int rank) {
	int i;
	if (rank == 0) {
		for (i = 0; i < TASKS_LIST_SIZE; i++) {
			tasks_list[i] = (SIZE_MULTIPLIER * i) % size;
		}
	}
	MPI_Bcast(tasks_list, TASKS_LIST_SIZE, MPI_INT, 0, comm);
}

void* worker() {
	return NULL;
}

void listener() {
	// logical array
	// if i-th element is set to 1 then such task is calculated
	int* tasks_status = (int*) malloc(TASKS_LIST_SIZE * sizeof(int));
	while (1) {

	}
}
