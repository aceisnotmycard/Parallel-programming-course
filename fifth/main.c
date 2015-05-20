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
int is_worker_ready;
// mutex
pthread_mutex_t mutex;
// array with sizes of tasks
int* tasks_list;
// current element of tasks_list 
// that should be processed by worker
// if it is equals to -1 than work is done 
int current_task;

// prepares data in 0 process
// TODO: remove argument
void prepare_tasks(int rank);

void* worker();

void listener();

// returns -1 if there is no free tasks
int find_neartest_free_task(int* tasks_status);

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
	pthread_create(&working_thread, &attrs, worker, NULL);
	pthread_attr_destroy(&attrs);

	// start listening for messages from another processes
	listener();

	// ok, going home
	if (pthread_join(working_thread, NULL) != 0) {
		perror("Cannot join thread");
		abort();
	}
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
	current_task = rank;
	MPI_Bcast(tasks_list, TASKS_LIST_SIZE, MPI_INT, 0, comm);
}

void* worker() {
	int i;
	int task;
	task = current_task;
	pthread_mutex_lock(&mutex);
	is_worker_ready = 0;
	pthread_mutex_unlock(&mutex);

	while (task != -1) {
		int task_size = tasks_list[task];
		for (i = 0; i < task_size; i++) {
			sqrt(i);
		}
		printf("Worker at %d process calculated value for %d task\n", rank, task);

		pthread_mutex_lock(&mutex);
		is_worker_ready = 0;
		pthread_mutex_unlock(&mutex);

		while(task == current_task) {}
		task = current_task;
	}

	return NULL;
}

void listener() {
	int process;
	int worker_state;
	//int free_task;
	// logical array
	// if i-th element is set to 1 then such task is calculated
	int* tasks_status = (int*) malloc(TASKS_LIST_SIZE * sizeof(int));

	// current task of i-th process
	int* current_tasks = (int*) malloc(size * sizeof(int));
	tasks_status[rank] = 1;

	// if current_state is set to 1 than work is completed
	int current_state;

	int everyone_is_completed;
	// main loop
	while (1) {
		// Exchanging information with other processes
		current_state = (current_task == -1) ? 1 : 0;
		MPI_Allreduce(&current_state, &everyone_is_completed, 1, MPI_INT, MPI_PROD, comm);
		if (everyone_is_completed) {
			break;
		}
		MPI_Allgather(&current_task, 1, MPI_INT, current_tasks, 1, MPI_INT, comm);
		for (process = 0; process < size; process++) {
			tasks_status[current_tasks[process]] = 1;
		}
		pthread_mutex_lock(&mutex);
		worker_state = is_worker_ready;
		pthread_mutex_unlock(&mutex);

		// if worker is ready then find available task
		if (worker_state == 0) {
			printf("Worker at %d process taking %d task\n", rank, current_task);
			current_task = find_neartest_free_task(tasks_status);
			pthread_mutex_lock(&mutex);
				is_worker_ready = 1;
			pthread_mutex_unlock(&mutex);
		}
	}
}

// returns -1 if there is no free tasks
int find_neartest_free_task(int* tasks_status) {
	int i;
	for (i = 0; i < TASKS_LIST_SIZE; i++) {
		if(tasks_status[i] == 0) {
			return i;
		}
	}
	return -1;
}