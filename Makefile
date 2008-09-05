#
# Copyright 2008 by 
# 
# Laboratoire de l'Informatique du Parallélisme, 
# UMR CNRS - ENS Lyon - UCB Lyon 1 - INRIA 5668
#
# Compile script for Sollya vandercoeff external procedure
#
# This software is based on scientific research made by
# Serge Torres and Nicolas Brisebarre.
#
# Contributors: 
#        Serge Torres (ENS Lyon) -- serge.torres@ens-lyon.fr
#        Christoph Lauter (ENS Lyon) -- christoph.lauter@ens-lyon.fr
#
# This file is an integrated part of the metalibm library developed by the 
# Arenaire project at Ecole Normale Superieure de Lyon
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation; either version 2 of the License, or 
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
#

# This makefile is a very provisonal version that
# will be sperceeded by an autotools generated version.
#
# TODO Makefile.in
#
DATE=`date +%Y-%m-%d`
# Clear all the suffixes
.SUFFIXES:

CXX       = g++
CC        = gcc
CFLAGS    = -fPIC -Wall -g -DNDEBUG 
CXXFLAGS  = -fPIC -Wall 
LDLIBS    = -lmpfr -lgmp -lm -lstdc++ 

INCS      = utils.h \
            ./vectors/mpfr/mpfr-vector.h \
            ./vectors/mpq/mpq-vector.h \
            ./vectors/mpz/mpz-vector.h \
            ./vectors/si/si-vector.h \
            ./vectors/ui/ui-vector.h \
            ./matrices/mpfr/mpfr-matrix.h \
            ./matrices/mpq/mpq-matrix.h \
            ./matrices/mpz/mpz-matrix.h \
            ./matrices/ui/ui-matrix.h \
            ./misc-utils/std-exit-errors.h \
	    ./libVanderCoeffsSparse.h \
	    ./mpfr++.h

SOURCES   = vanderCoeffsSparseAbsErr.c \
            vanderCoeffsSparseRelErr.c \
            libVanderCoeffsSparse.c \
            vandercoeffs.c \
            utils.h utils.c \
            mpfr++.cpp


OBJS      = utils.o \
            libVanderCoeffsSparse.o \
            ./vectors/mpfr/mpfr-vector.o \
            ./vectors/mpq/mpq-vector.o \
            ./vectors/mpz/mpz-vector.o \
            ./vectors/si/si-vector.o \
            ./vectors/ui/ui-vector.o \
            ./matrices/mpfr/mpfr-matrix.o \
            ./matrices/mpq/mpq-matrix.o \
            ./matrices/mpz/mpz-matrix.o \
            ./matrices/ui/ui-matrix.o \
            mpfr++.o


# ---------------------------------------------------------------------


all: vandercoeff 

./vectors/mpfr/mpfr-vector.o: vectors/mpfr/mpfr-vector.h \
                              vectors/mpfr/mpfr-vector.c 
	$(MAKE) -C vectors/mpfr

./vectors/mpq/mpq-vector.o: vectors/mpq/mpq-vector.h \
                            vectors/mpq/mpq-vector.c 
	$(MAKE) -C vectors/mpq

./vectors/mpz/mpz-vector.o: vectors/mpz/mpz-vector.h \
                            vectors/mpz/mpz-vector.c 
	$(MAKE) -C vectors/mpz

./vectors/si/si-vector.o: vectors/si/si-vector.h \
                          vectors/si/si-vector.c 
	$(MAKE) -C vectors/si

./vectors/ui/ui-vector.o: vectors/ui/ui-vector.h \
                          vectors/ui/ui-vector.c 
	$(MAKE) -C vectors/ui

./matrices/mpfr/mpfr-matrix.o: matrices/mpfr/mpfr-matrix.h \
                               matrices/mpfr/mpfr-matrix.c 
	$(MAKE) -C matrices/mpfr

./matrices/mpq/mpq-matrix.o: matrices/mpq/mpq-matrix.h \
                             matrices/mpq/mpq-matrix.c 
	$(MAKE) -C matrices/mpq


./matrices/mpz/mpz-matrix.o: matrices/mpz/mpz-matrix.h \
                             matrices/mpz/mpz-matrix.c 
	$(MAKE) -C matrices/mpz

./matrices/ui/ui-matrix.o: matrices/ui/ui-matrix.h \
                           matrices/ui/ui-matrix.c 
	$(MAKE) -C matrices/ui

vandercoeff.o: vandercoeff.c $(INCS) $(OBJS)
	$(CC) $(CFLAGS) -c $< -o $@

vandercoeff : vandercoeff.o libVanderCoeffsSparse.o $(OBJS)
	$(CC) $(CFLAGS) -shared $(OBJS) vandercoeff.o  -L. $(LDLIBS) -o $@

mpfr++.o: mpfr++.cpp mpfr++.h 
	$(CXX) $(CFLAGS) -c -o mpfr++.o mpfr++.cpp

%.o: %.c %.h $(INCS)
	$(CC) $(CFLAGS) -c -o $@ $< 

% : %.o $(OBJS)
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

# ---------------------------------------------------------------------

.PHONY: clean scratch

clean:
	@rm -f *.o a.out
	@rm -f *~ *% #*#
	@rm -rf html
	@$(MAKE) -C vectors/mpfr clean   > /dev/null
	@$(MAKE) -C vectors/mpq clean    > /dev/null
	@$(MAKE) -C vectors/mpz clean    > /dev/null
	@$(MAKE) -C vectors/si clean     > /dev/null
	@$(MAKE) -C vectors/ui clean     > /dev/null
	@$(MAKE) -C matrices/mpfr clean  > /dev/null
	@$(MAKE) -C matrices/mpq clean   > /dev/null
	@$(MAKE) -C matrices/mpz clean   > /dev/null
	@$(MAKE) -C matrices/ui clean    > /dev/null

scratch: clean
	rm -f $(SCRATCH_ALL)

