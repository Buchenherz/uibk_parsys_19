CC=gcc
CC_MAC=gcc-9
CC_FLAGS_OPTIMIZED= -g -O3 -Wall -Wextra -march=native -mtune=native -ftree-vectorize -lm -fopenmp
CC_FLAGS_NON_OPTIMIZED= -O0 -Wall -Wextra -g -lm -fopenmp
CC_FLAGS_PROFILING_CLION= -fno-omit-frame-pointer -g -O0 -Wall -Wextra -march=native -mtune=native -ftree-vectorize -lprofiler -lm -fopenmp
# you could also use -ftree-vectorizer-verbose=1 or 2 to get verbose output of vectorization
# 02 might give better performance
.PHONEY: all
all: real-optimized real-non-optimized seq
clion: real-clion

SRC=c_timers.c print_results.c randdp.c real.c wtime.c
SEQ_SRC=c_timers.c print_results.c randdp.c real_seq.c wtime.c
binaries=real_optimized real_non_optimized real

seq: $(SEQ_SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real_seq $(SEQ_SRC)

seq-mac: $(SEQ_SRC)
	@$(CC_MAC) $(CC_FLAGS_OPTIMIZED) -o real_seq $(SEQ_SRC)

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
