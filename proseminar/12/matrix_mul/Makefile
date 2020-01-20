CC=chpl

SRC=matrix_mul.chpl
binaries1=matrix_mul

all:  matrix

matrix: $(SEQ_SRC)
	@$(CC) $(SRC) -o matrix_mul

.PHONEY: clean
clean:
	@rm -f $(binaries)

.PHONEY: run
run: all
	@./matrix_mul

.PHONEY: debug
debug: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./matrix_mul --debug=true --rows_and_columns=4

.PHONEY: 1thread
1thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=1 ./matrix_mul

.PHONEY: 2thread
2thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=2 ./matrix_mul

.PHONEY: 4thread
4thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=4 ./matrix_mul

.PHONEY: 8thread
8thread: matrix
	@CHPL_RT_NUM_THREADS_PER_LOCALE=8 ./matrix_mul