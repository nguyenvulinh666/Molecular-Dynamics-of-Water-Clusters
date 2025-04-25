SHELL = /bin/sh
EXE = SKWPRG.X
# set FC to f77 for unix and g77 for cygwin
FC = g77
FFLAGS = -c
# set CC to cc for unix and gcc for cygwin 
CC = gcc
CFLAGS = -c
# set LD to f77 for unix and g77 for cygwin 
LD = g77
LDFLAGS = -o $(EXE)
OBJ = skwprg.o dencor.o start.o forden.o correl.o denfft.o fft.o invert.o

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ)

o.c:
	$(CC) $(CFLAGS) $<

o.f:
	$(FC) $(FFLAGS) $<

invert.o: ../FORT/invert.f
	$(FC) $(FFLAGS) ../FORT/invert.f

clean: $(OBJ)
	rm $(OBJ)

