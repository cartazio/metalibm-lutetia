CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = 
OBJS      = si-vector.o

all: si-test-vector
 
si-test-vector: si-test-vector.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

si-vector.o: si-vector.c si-vector.h si-vector-types.h

clean:
	@rm -f si-test-vector *.o
	@rm -f *~ *% #*#
