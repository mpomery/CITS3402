#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -std=c99 -Wall -Werror -O3

#-----------------------------------------------------------------------------

all : $(BIN)

$(BIN) :
	$(CC) $(FLAGS) -o breather breather.c

clean:
	$(RM) $(DIRT)

clobber:
	$(RM) $(DIRT)
	$(RM) $(BIN)
