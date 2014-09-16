#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -Wall -Werror -fpermissive -O3

#-----------------------------------------------------------------------------

all : $(BIN)

$(BIN) :
	$(CC) $(FLAGS) -o breather breather.o

clean:
	$(RM) $(DIRT)

clobber:
	$(RM) $(DIRT)
	$(RM) $(BIN)
