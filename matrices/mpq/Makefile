CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = -lmpfr -lgmp  
OBJS      = mpq-matrix.o ../../vectors/mpq/mpq-vector.o

all: mpq-test-matrix
 
mpq-test-matrix: mpq-test-matrix.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

mpq-matrix.o: mpq-matrix.c mpq-matrix.h mpq-matrix-types.h

clean:
	@rm -f mpq-test-matrix *.o
	@rm -f *~ *% #*#

