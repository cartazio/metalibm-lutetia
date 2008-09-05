CC        = gcc
CFLAGS    = -fPIC -Wall -g -static
LDLIB     = -lmpfr -lgmp  
OBJS      = mpfr-matrix.o ../../vectors/mpfr/mpfr-vector.o

all: mpfr-test-matrix
 
mpfr-test-matrix: mpfr-test-matrix.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

mpfr-matrix.o: mpfr-matrix.c mpfr-matrix.h mpfr-matrix-types.h

clean:
	@rm -f mpfr-test-matrix *.o
	@rm -f *~ *% #*#	