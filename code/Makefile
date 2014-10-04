#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -Wall -Werror -fopenmp
FLAGS2 = -std=c99
COUNT = | wc -l

#-----------------------------------------------------------------------------

all : clean $(BIN) new original test

new: clean $(BIN)

original:
	cd orig;\
	$(CC) $(FLAGS2) -o breather breather.c;\

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
	time ./breather
	time ./orig/breather
	diff acce.dat orig/acce.dat $(COUNT)
	acce.dat | wc -l
	orig/acce.dat | wc -l
	diff cmass.dat orig/cmass.dat $(COUNT)
	cmass.dat | wc -l
	orig/cmass.dat | wc -l
	diff ke.dat orig/ke.dat $(COUNT)
	ke.dat | wc -l
	orig/ke.dat | wc -l
	diff pe.dat orig/pe.dat $(COUNT)
	pe.dat | wc -l
	orig/pe.dat | wc -l
	diff restart.dat orig/restart.dat $(COUNT)
	restart.dat | wc -l
	orig/restart.dat | wc -l
	diff strsh.dat orig/strsh.dat $(COUNT)
	strsh.dat | wc -l
	orig/strsh.dat | wc -l
	diff toten.dat orig/toten.dat $(COUNT)
	toten.dat | wc -l
	orig/toten.dat | wc -l
	diff velsh.dat orig/velsh.dat $(COUNT)
	velsh.dat | wc -l
	orig/velsh.dat | wc -l
