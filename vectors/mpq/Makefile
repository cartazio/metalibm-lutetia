CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = -lmpfr -lgmp 
OBJS      = mpq-vector.o

all: mpq-test-vector
 
mpq-test-vector: mpq-test-vector.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

mpq-vector.o: mpq-vector.c mpq-vector.h mpq-vector-types.h

clean:
	@rm -f mpq-test-vector *.o
	@rm -f *~ *% #*#
