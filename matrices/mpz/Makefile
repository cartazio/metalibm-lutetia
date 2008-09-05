CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = -lmpfr -lgmp  
OBJS      = mpz-matrix.o

all: mpz-test-matrix
 
mpz-test-matrix: mpz-test-matrix.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

mpz-matrix.o: mpz-matrix.c mpz-matrix.h mpz-matrix-types.h

clean:
	@rm -f mpz-test-matrix *.o
	@rm -f *~ *% #*#

