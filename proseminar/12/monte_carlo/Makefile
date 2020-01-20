CC=chpl

SRC=monte_carlo.chpl
binaries=monte_pi

all: monte

monte: monte_carlo.chpl
	@$(CC) $(SRC) -o monte_pi

fast-monte: fast_monte_carlo.chpl
	@$(CC) fast_monte_carlo.chpl -o monte_pi

.PHONEY: clean
clean:
	@rm -f $(binaries)

.PHONEY: run
run: all
	@./monte_pi
	
.PHONEY: debug
debug: monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./monte_pi --debug=true --number_of_points=100000000 --num_tasks=4

.PHONEY: 1thread
1thread: monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=1 ./monte_pi --num_tasks=1 --debug=true

.PHONEY: 2thread
2thread: monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=2 ./monte_pi --num_tasks=2 --debug=true

.PHONEY: 4thread
4thread: monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./monte_pi --num_tasks=4 --debug=true

.PHONEY: 8thread
8thread: monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8 ./monte_pi --num_tasks=8 --debug=true

.PHONEY: fast-debug
fast-debug: fast-monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./monte_pi

.PHONEY: fast-1thread
fast-1thread: fast-monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=1 ./monte_pi

.PHONEY: fast-2thread
fast-2thread: fast-monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=2 ./monte_pi

.PHONEY: fast-4thread
fast-4thread: fast-monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./monte_pi

.PHONEY: fast-8thread
fast-8thread: fast-monte
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8 ./monte_pi