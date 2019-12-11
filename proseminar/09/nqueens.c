#include <arpa/inet.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// #region https://jameshfisher.com/2017/02/23/printing-bits/
void print_byte_as_bits(char val) {
    for (int i = 7; 0 <= i; i--) {
        printf("%c", (val & (1 << i)) ? '1' : '0');
    }
}

void print_bits(char *ty, char *val, unsigned char *bytes, size_t num_bytes) {
    printf("(%*s) %*s = [ ", 15, ty, 16, val);
    for (size_t i = 0; i < num_bytes; i++) {
        print_byte_as_bits(bytes[i]);
        printf(" ");
    }
    printf("]\n");
}

#define SHOW(T, V)                                          \
    do {                                                    \
        T x = V;                                            \
        print_bits(#T, #V, (unsigned char *)&x, sizeof(x)); \
    } while (0)
// #endregion
#define MAXN 31

int nqueens(int n) {
    int q0, q1;
    int cols[MAXN], diagl[MAXN], diagr[MAXN],
        posibs[MAXN];  // Our backtracking 'stack'
    int num = 0;
    //
    // The top level is two fors, to save one bit of symmetry in the enumeration
    // by forcing second queen to be AFTER the first queen.
    //
    for (q0 = 0; q0 < n - 2; q0++) {
        for (q1 = q0 + 2; q1 < n; q1++) {
            int bit0 = 1 << q0;
            SHOW(int, bit0);
            int bit1 = 1 << q1;
            SHOW(int, bit1);
            SHOW(int, bit0 | bit1);
            SHOW(int, (-1 << n));
            int d = 0;  // d is our depth in the backtrack stack
            cols[0] =
                bit0 | bit1 | (-1 << n);  // The -1 here is used to fill all
                                          // 'coloumn' bits after n ...
            SHOW(int, cols[0]);
            diagl[0] = (bit0 << 1 | bit1) << 1;
            SHOW(int, diagl[0]);
            diagr[0] = (bit0 >> 1 | bit1) >> 1;
            SHOW(int, diagr[0]);
            //  The variable posib contains the bitmask of possibilities we
            //  still have to try in a given row ...
            int posib = ~(cols[0] | diagl[0] | diagr[0]);

            while (d >= 0) {
                while (posib) {
                    int bit = posib & -posib;  // The standard trick for getting
                                               // the rightmost bit in the mask
                    int ncols = cols[d] | bit;
                    int ndiagl = (diagl[d] | bit) << 1;
                    int ndiagr = (diagr[d] | bit) >> 1;
                    int nposib = ~(ncols | ndiagl | ndiagr);
                    posib ^= bit;  // Eliminate the tried possibility.

                    // The following is the main additional trick here, as
                    // recognizing solution can not be done using stack level
                    // (d), since we save the depth+backtrack time at the end of
                    // the enumeration loop. However by noticing all coloumns
                    // are filled (comparison to -1) we know a solution was
                    // reached ... Notice also that avoiding an if on the
                    // ncols==-1 comparison is more efficient!
                    num += ncols == -1;

                    if (nposib) {
                        if (posib) {  // This if saves stack depth + backtrack
                                      // operations when we passed the last
                                      // possibility in a row.
                            posibs[d++] = posib;  // Go lower in stack ..
                        }
                        cols[d] = ncols;
                        diagl[d] = ndiagl;
                        diagr[d] = ndiagr;
                        posib = nposib;
                    }
                }
                posib = posibs[--d];  // backtrack ...
            }
        }
    }
    return num * 2;
}

int main(int ac, char **av) {
    if (ac != 2) {
        printf("usage: nq n\n");
        return 1;
    }
    int n = atoi(av[1]);
    if (n < 1 || n > MAXN) {
        printf("n must be between 2 and 31!\n");
    }
    printf("Number of solution for %d is %d\n", n, nqueens(n));
}