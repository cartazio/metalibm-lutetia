CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = -lmpfr -lgmp  
OBJS      = mpfr-vector.o

all: mpfr-test-vector
 
mpfr-test-vector: mpfr-test-vector.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

mpfr-vector.o: mpfr-vector.c mpfr-vector.h mpfr-vector-types.h

clean:
	@rm -f mpfr-test-vector *.o
	@rm -f *~ *% #*#
