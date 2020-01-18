CC=chpl

SRC1=monte_carlo.chpl
SRC2=matrix_mul.chpl
binaries1=monte_pi matrix_mul

all: monte matrix

monte: $(SEQ_SRC)
	@$(CC) $(SRC1) -o monte_pi 

matrix: $(SEQ_SRC)
	@$(CC) $(SRC2) -o matrix_mul

.PHONEY: clean
clean:
	@rm -f $(binaries)

.PHONEY: run
run: all
	@./monte_pi
	@./matrix_mul
	
.PHONEY: debugMonte
debugMonte: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=4
	@./monte_pi --debug=true --number_of_points=10000000 --num_tasks=4
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: debugMatrix
debugMatrix: matrix
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=4
	@./matrix_mul --debug=true --rows_and_columns=2552 --num_tasks=4
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 1threadMonte
1threadMonte: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=1
	@./monte_pi --num_tasks=1
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8
.PHONEY: 1threadMatrix
1threadMatrix: matrix
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=1
	@./matrix_mul --num_tasks=1
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 2threadMonte
2threadMonte: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=2
	@./monte_pi --num_tasks=2
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8
.PHONEY: 2threadMatrix
2threadMatrix: matrix
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=2
	@./matrix_mul --num_tasks=2
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 4threadMonte
4threadMonte: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=4
	@./monte_pi --num_tasks=4
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8
.PHONEY: 4threadMatrix
4threadMatrix: matrix
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=4
	@./matrix_mul --num_tasks=4
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8

.PHONEY: 8threadMonte
8threadMonte: monte
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=8
	@./monte_pi --num_tasks=8
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8
.PHONEY: 8threadMatrix
8threadMatrix: matrix
	@export CHPL_RT_NUM_THREADS_PER_LOCALE=8
	@./matrix_mul --num_tasks=8
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8