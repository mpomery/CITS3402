#! /usr/bin/make

CC=gcc

BIN = breather
OBJECTS = breather.o
DIRT = *.o *.i *~ */*~ *.log *.exe
FLAGS = -Wall -fpermissive -O3

#-----------------------------------------------------------------------------

all : $(BIN)

$(BIN) : $(OBJECTS) $(HEADERS)
	$(CC) $(FLAGS) -o breather breather.o

clean:
	$(RM) $(DIRT)

clobber:
	$(RM) $(DIRT)
	$(RM) $(BIN)