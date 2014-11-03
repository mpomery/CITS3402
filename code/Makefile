#! /usr/bin/make

CC=/usr/lib64/openmpi/bin/mpicc
MPIEXEC=/usr/lib64/openmpi/bin/mpiexec

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -Wall -Werror -O3 -fopenmp
FLAGS2 = -std=c99 -O3
COUNT = | wc -l

#-----------------------------------------------------------------------------

all : clean parallel mpitest

original:
	cd orig;\
	$(CC) $(FLAGS2) -o breather breather.c;\

parallel:
	$(CC) $(FLAGS) -o breather breather.c

mpitest:
	time $(MPIEXEC) -hostfile hosts.txt -np 1 ./breather

clean:
	$(RM) $(DIRT)
	$(RM) $(BIN)

clobber:
	$(RM) $(DIRT)
	$(RM) $(BIN)
	$(RM) *.dat

testparallel:
	./breather

test:
	time ./breather
	time ./orig/breather
	diff acce.dat orig/acce.dat $(COUNT)
	cat acce.dat | wc -l
	cat orig/acce.dat | wc -l
	diff cmass.dat orig/cmass.dat $(COUNT)
	cat cmass.dat | wc -l
	cat orig/cmass.dat | wc -l
	diff ke.dat orig/ke.dat $(COUNT)
	cat ke.dat | wc -l
	cat orig/ke.dat | wc -l
	diff pe.dat orig/pe.dat $(COUNT)
	cat pe.dat | wc -l
	cat orig/pe.dat | wc -l
	diff restart.dat orig/restart.dat $(COUNT)
	cat restart.dat | wc -l
	cat orig/restart.dat | wc -l
	diff strsh.dat orig/strsh.dat $(COUNT)
	cat strsh.dat | wc -l
	cat orig/strsh.dat | wc -l
	diff toten.dat orig/toten.dat $(COUNT)
	cat toten.dat | wc -l
	cat orig/toten.dat | wc -l
	diff velsh.dat orig/velsh.dat $(COUNT)
	cat velsh.dat | wc -l
	cat orig/velsh.dat | wc -l
