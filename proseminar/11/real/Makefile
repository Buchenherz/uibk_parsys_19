CC=gcc
CC_FLAGS_OPTIMIZED= -g -O3 -lm -march=native -mtune=native -ftree-vectorize -fopenmp
CC_FLAGS_NON_OPTIMIZED= -g -lm -fopenmp
CC_FLAGS_PROFILING_CLION= -fno-omit-frame-pointer -g -O0 -lm -march=native -mtune=native -ftree-vectorize -fopenmp

.PHONEY: all
all: real-optimized-seq real-non-optimized-seq

SRC=$(patsubst %.c,%.o,$(wildcard *.c))
binaries=real_optimized_seq real_non_optimized_seq real

real-optimized-seq: $(SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real_optimized_seq $(SRC)

real-non-optimized-seq: $(SRC)
	@$(CC) $(CC_FLAGS_NON_OPTIMIZED) -o real_non_optimized_seq $(SRC)

real-optimized-seq-clion: $(SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real $(SRC)


.PHONEY: clean
clean:
	@rm -f $(binaries) *.0 *.bin

.PHONEY: run
run: all
	@echo "Sequential:"
	@./real_optimized_seq.bin
	@echo
