CC=gcc
CC_FLAGS= -g -O3 -lm -lprofiler

.PHONEY: all
all: real

SRC=$(patsubst %.c,%.o,$(wildcard *.c))

real: $(SRC)
	@$(CC) $(CC_FLAGS) -o real.bin $(SRC) 

.PHONEY: clean
clean:
	@rm *.bin
	@rm *.o

.PHONEY: run
run: all
	@echo "Sequential:"
	@./real
	@echo
