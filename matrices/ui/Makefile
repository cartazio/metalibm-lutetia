CC        = gcc
CFLAGS    = -fPIC -Wall -g -static 
LDLIB     = 
OBJS      = ui-matrix.o \
            ../../vectors/ui/ui-vector.o

all: ui-test-matrix
 
ui-test-matrix: ui-test-matrix.c $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIB) -o $@ 

ui-matrix.o: ui-matrix.c ui-matrix.h ui-matrix-types.h

clean:
	@rm -f ui-test-matrix *.o
	@rm -f *~ *% #*#
