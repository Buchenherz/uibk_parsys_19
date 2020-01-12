CC=gcc
CC_MAC=gcc-9
CC_FLAGS_OPTIMIZED= -g -O3 -march=native -mtune=native -ftree-vectorize -lm -fopenmp
CC_FLAGS_NON_OPTIMIZED= -O3 -g -lm -fopenmp
CC_FLAGS_PROFILING_CLION= -fno-omit-frame-pointer -g -O0  -march=native -mtune=native -ftree-vectorize -lm -fopenmp

.PHONEY: all
all: real-optimized-seq real-non-optimized-seq real-optimized-seq-clion

SRC=$(wildcard *.c)
binaries=real_optimized_seq real_non_optimized_seq real

real-optimized-seq: $(SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real_optimized_seq $(SRC)

real-non-optimized-seq: $(SRC)
	@$(CC) $(CC_FLAGS_NON_OPTIMIZED) -o real_non_optimized_seq $(SRC)

real-optimized-seq-clion: $(SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real $(SRC)

mac: $(SRC)
	@$(CC_MAC) $(CC_FLAGS_OPTIMIZED) -o real $(SRC)

.PHONEY: clean
clean:
	@rm -f $(binaries) *.0 *.bin *.o

.PHONEY: run
run: all
	@echo "Sequential:"
	@./real_optimized_seq.bin
	@echo
