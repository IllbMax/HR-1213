# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS =  -pedantic -Wall -Wextra -O2 -ggdb -g
LFLAGS = $(CFLAGS)

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: timempi timempi2

timempi: timempi.o Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS)

timempi2: timempi2.o Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS)

clean:
	$(RM) *.o *~

timempi.o: timempi.c Makefile

timempi2.o: timempi2.c Makefile

