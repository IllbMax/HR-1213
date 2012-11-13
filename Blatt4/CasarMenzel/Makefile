# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O2 -ggdb -fopenmp
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff-seq.o askparams.o displaymatrix.o
# Objektvariablen für das Aufgabenblatt
OBJS_OMP = partdiff-openmp.o askparams.o displaymatrix.o
OBJS_OMP_S = partdiff-openmp_spalten.o askparams.o displaymatrix.o
OBJS_OMP_E = partdiff-openmp_element.o askparams.o displaymatrix.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-seq partdiff-openmp

partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

# Target für Aufgabe 1
partdiff-openmp: $(OBJS_OMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS_OMP) $(LIBS)

# Targets für Aufgabe 2
partdiff-openmp-element: $(OBJS_OMP_E) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS_OMP_E) $(LIBS)

partdiff-openmp-spalten: $(OBJS_OMP_S) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS_OMP_S) $(LIBS)

partdiff-openmp-zeilen: $(OBJS_OMP) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS_OMP) $(LIBS)

clean:
	$(RM) *.o *~

partdiff-seq.o: partdiff-seq.c Makefile

partdiff-openmp.o: partdiff-openmp.c Makefile
partdiff-openmp_spalten.o: partdiff-openmp_element.c Makefile
partdiff-openmp_element.o: partdiff-openmp_spalten.c Makefile

askparams.o: askparams.c Makefile

displaymatrix.o: displaymatrix.c Makefile
