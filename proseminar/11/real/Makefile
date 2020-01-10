CC=gcc-9
CC_FLAGS_OPTIMIZED= -O3 -lm -march=native -mtune=native -ftree-vectorize -fopenmp
CC_FLAGS_NON_OPTIMIZED= -O3 -lm -fopenmp

.PHONEY: all
all: real-optimized real-non-optimized

SRC=$(patsubst %.c,%.o,$(wildcard *.c))

real-optimized: $(SRC)
	@$(CC) $(CC_FLAGS_OPTIMIZED) -o real_optimized.bin $(SRC) 

real-non-optimized: $(SRC)
	@$(CC) $(CC_FLAGS_NON_OPTIMIZED) -o real_non_optimized.bin $(SRC) 

.PHONEY: clean
clean:
	@rm *.bin
	@rm *.o

.PHONEY: run
run: all
	@echo "Sequential:"
	@./real
	@echo
