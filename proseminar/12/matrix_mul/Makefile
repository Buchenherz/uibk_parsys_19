CC=chpl

SRC=matrix_mul_singlenode.chpl
binaries1=matrix_mul_singlenode

all:  matrix

matrix: $(SEQ_SRC)
	@$(CC) $(SRC) -o matrix_mul_singlenode

.PHONEY: clean
clean:
	@rm -f $(binaries)

.PHONEY: run
run: all
	@./matrix_mul_singlenode

.PHONEY: debug
debug: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./matrix_mul_singlenode --debug=true --rows_and_columns=4

.PHONEY: 1thread
1thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=1 ./matrix_mul_singlenode

.PHONEY: 2thread
2thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=2 ./matrix_mul_singlenode

.PHONEY: 4thread
4thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./matrix_mul_singlenode

.PHONEY: 8thread
8thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8 ./matrix_mul_singlenode