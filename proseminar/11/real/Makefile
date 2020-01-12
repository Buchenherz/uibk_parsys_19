CC=gcc
CC_MAC=gcc-9
CC_FLAGS_OPTIMIZED= -g -O3 -march=native -mtune=native -ftree-vectorize -lm -fopenmp
CC_FLAGS_NON_OPTIMIZED= -O3 -g -lm -fopenmp
CC_FLAGS_PROFILING_CLION= -fno-omit-frame-pointer -g -O0 -march=native -mtune=native -ftree-vectorize -lm -fopenmp

.PHONEY: all
all: real-optimized real-non-optimized 
clion: real-clion

SRC=$(wildcard *.c)
binaries=real_optimized real_non_optimized real

real-optimized: $(SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real $(SRC)

real-non-optimized: $(SRC)
	@$(CC) $(CC_FLAGS_NON_OPTIMIZED) -o real_non_optimized $(SRC)

real-clion: $(SRC)
	@$(CC_MAC) $(CC_FLAGS_PROFILING_CLION) -o real $(SRC)

mac: $(SRC)
	@$(CC_MAC) $(CC_FLAGS_OPTIMIZED) -o real $(SRC)

.PHONEY: clean
clean:
	@rm -f $(binaries) *.0 *.bin *.o

.PHONEY: run
run: all
	@echo "Sequential:"
	@./real
	@echo
