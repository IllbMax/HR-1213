# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -pthread -pedantic -Wall -Wextra -O2 -ggdb -g
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-posix.o askparams.o displaymatrix.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-posix

partdiff-posix: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	$(RM) *.o *~

partdiff-posix.o: partdiff-posix.c Makefile

askparams.o: askparams.c Makefile

displaymatrix.o: displaymatrix.c Makefile
