#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -O3 -fopenmp

#-----------------------------------------------------------------------------

all : clean $(BIN) test

original:
	cd orig;\
	ls;\
	$(CC) $(FLAGS) -o breather breather.c;\
	./breather

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
	diff acce.dat orig/acce.dat | wc -l
	diff cmass.dat orig/cmass.dat | wc -l
	diff ke.dat orig/ke.dat | wc -l
	diff pe.dat orig/pe.dat | wc -l
	diff restart.dat orig/restart.dat | wc -l
	diff strsh.dat orig/strsh.dat | wc -l
	diff toten.dat orig/toten.dat | wc -l
	diff velsh.dat orig/velsh.dat | wc -l
