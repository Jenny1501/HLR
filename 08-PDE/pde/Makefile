# Common definitions
CC = mpicc

# Compiler flags, paths and libraries
CFLAGS = -std=c99 -pedantic -Wall -Wextra -O0 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS1 = partdiff-par.o askparams.o

# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: partdiff-par

partdiff-par: $(OBJS1) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS1) $(LIBS)

clean:
	$(RM) *.o *~

partdiff-par.o: partdiff-par.c Makefile

askparams.o: askparams.c Makefile
