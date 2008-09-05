CC        = gcc
CFLAGS    = -fPIC -Wall -g -static
LDLIB     = 
OBJS      = ui-vector.o

all: ui-test-vector
 
ui-test-vector: ui-test-vector.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

ui-vector.o: ui-vector.c ui-vector.h ui-vector-types.h

clean:
	@rm -f ui-test-vector *.o
	@rm -f *~ *% #*#
