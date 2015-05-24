#include <mpi.h>
#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define TASKS_LIST_SIZE 20
#define SIZE_MULTIPLIER 1234000
#define N 1000

// for MPI
int rank, size;
MPI_Comm comm;
// logical variable
// 1 if worker is done and waiting for more
// 0 if it is busy
int iter_counter;
int is_worker_ready;
// mutex
pthread_mutex_t mutex;

void* task_parser();
void task_getter();
int get_task_from_queue(int iteration);

// returns -1 if there is no free tasks
// int find_neartest_free_task(int* tasks_status);

int main(int argc, char** argv) {
	int provided;

	// prepare MPI
	comm = MPI_COMM_WORLD;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_size(comm, &size);
	MPI_Comm_rank(comm, &rank);

	//prepare data
	//tasks_list = (int*) malloc(TASKS_LIST_SIZE * sizeof(int));
	//prepare_tasks(rank);

	//prepare working thread
	pthread_t working_thread;
	pthread_attr_t attrs;

	if (pthread_attr_init(&attrs) != 0) {
		perror("Cannot init thread attrs");
		abort();
	}
	if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE) != 0) {
		perror("Cannot set thread joinable");
		abort();
	}
	if (pthread_mutex_init(&mutex, NULL) != 0) {
		perror("Cannot init mutex");
		abort();	
	}
	pthread_create(&working_thread, &attrs, task_parser, NULL);
	pthread_attr_destroy(&attrs);

	// start listening for messages from another processes
	task_getter();

	// ok, going home
	if (pthread_join(working_thread, NULL) != 0) {
		perror("Cannot join thread");
		abort();
	}
	MPI_Finalize();
	return 0;
}

void task_getter() {
	int* states = (int*) malloc(sizeof(int) * size);
	int worker_state;
	int current_state;
	int everyone_is_completed = 0;
	int sum;
	int task;
	int counter;
	current_task = rank;
	iter_counter = 1;
	while(1) {
		
		// check is everyone complete
		current_state = (current_task == -1) ? 1 : 0;
		MPI_Allreduce(&current_state, &everyone_is_completed, 1, MPI_INT, MPI_PROD, comm);
		if (everyone_is_completed) {
			printf("bye-bye\n");
			break;
		}

		//check how many processes complete
		pthread_mutex_lock(&mutex);
		worker_state = is_worker_ready;
		pthread_mutex_unlock(&mutex);

		MPI_Allgather(&worker_state, 1, MPI_INT, states, 1, MPI_INT, comm);

		//get new task
		counter = 0;
		for (int i = 0; i < size; i++) {
			if (states[i] == 1) {
				counter++;
				if (rank == i) {
					task = get_task_from_queue(iter_counter+counter);
					printf("%d\n", iter_counter + counter);
					counter = 1;
					break;
				}
			}
		}
		// get number of iterations
		MPI_Allreduce(&counter, &sum, 1, MPI_INT, MPI_SUM, comm);
		iter_counter += sum;

		if (worker_state) {
			pthread_mutex_lock(&mutex);
			current_task = task;
			is_worker_ready = 0;
			pthread_mutex_unlock(&mutex);
		}
	}
}

void* task_parser() {
	int task;
	double result;
	is_worker_ready = 0;
	int iwr;
	while (task != -1) {
		pthread_mutex_lock(&mutex);
		iwr = is_worker_ready;
		pthread_mutex_unlock(&mutex);
		if (!iwr) {
			pthread_mutex_lock(&mutex);
			task = current_task;
			pthread_mutex_unlock(&mutex);
			for (int i = 0; i < SIZE_MULTIPLIER * task; i++) {
				result+=sqrt(i);
			}
			pthread_mutex_lock(&mutex);
			is_worker_ready = 1;
			pthread_mutex_unlock(&mutex);
		}
	}
	return NULL;
}


int get_task_from_queue(int iteration) {
	return iteration < N ? abs((rank - (iteration%size))) : -1;
}