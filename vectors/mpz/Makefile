CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = -lmpfr -lgmp 
OBJS      = mpz-vector.o

all: mpz-test-vector
 
mpz-test-vector: mpz-test-vector.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

mpz-vector.o: mpz-vector.c mpz-vector.h mpz-vector-types.h

clean:
	@rm -f mpz-test-vector *.o
	@rm -f *~ *% #*#
