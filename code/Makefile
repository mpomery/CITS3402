#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -Wall -Werror -O3 -fopenmp

#-----------------------------------------------------------------------------

all : clean $(BIN) test

$(BIN) :
	$(CC) $(FLAGS) -o breather breather.c

clean:
	$(RM) $(DIRT)
	$(RM) $(BIN)

clobber:
	$(RM) $(DIRT)
	$(RM) $(BIN)
	$(RM) *.dat

test:
	./breather
	diff acce.dat orig/acce.dat
	diff cmass.dat orig/cmass.dat
	diff ke.dat orig/ke.dat
	diff pe.dat orig/pe.dat
	diff restart.dat orig/restart.dat
	diff strsh.dat orig/strsh.dat
	diff toten.dat orig/toten.dat
	diff velsh.dat orig/velsh.dat
