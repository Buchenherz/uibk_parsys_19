CC=chpl

SRC=monte_carlo.chpl
binaries=monte_pi

all: monte

monte: $(SEQ_SRC)
	@$(CC) $(SRC) -o monte_pi 

.PHONEY: clean
clean:
	@rm -f $(binaries)

.PHONEY: run
run: all
	@./monte_pi
	
.PHONEY: debug
debug: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=4
	@./monte_pi --debug=true --number_of_points=10000000 --num_tasks=4
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 1thread
1thread: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=1
	@./monte_pi --num_tasks=1
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 2thread
2thread: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=2
	@./monte_pi --num_tasks=2
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 4thread
4thread: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=4
	@./monte_pi --num_tasks=4
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 8thread
8thread: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=8
	@./monte_pi --num_tasks=8 --debug=true
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8