#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -Wall -Werror -fopenmp
FLAGS2 = -std=c99
COUNT = | wc -l

#-----------------------------------------------------------------------------

all : clean $(BIN) original new test

new: clean $(BIN)

original:
	cd orig;\
	ls;\
	$(CC) $(FLAGS2) -o breather breather.c;\
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
	diff acce.dat orig/acce.dat $(COUNT)
	diff cmass.dat orig/cmass.dat $(COUNT)
	diff ke.dat orig/ke.dat $(COUNT)
	diff pe.dat orig/pe.dat $(COUNT)
	diff restart.dat orig/restart.dat $(COUNT)
	diff strsh.dat orig/strsh.dat $(COUNT)
	diff toten.dat orig/toten.dat $(COUNT)
	diff velsh.dat orig/velsh.dat $(COUNT)
