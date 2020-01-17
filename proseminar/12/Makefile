CC=chpl

SRC=monte_carlo.chpl
binaries=monte_pi

all: monte

monte: $(SEQ_SRC)
	@$(CC) $(SRC) -o monte_pi 

.PHONEY: clean
clean:
	@rm -f $(binaries)

.PHONEY: run
run: all
	@./monte_pi
	
.PHONEY: debug
debug: monte
	@./monte_pi --debug=true
