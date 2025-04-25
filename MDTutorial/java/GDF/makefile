SHELL = /bin/sh
EXE = GDFPRG.X
# set FC to f77 for unix and g77 for cygwin
FC = g77
FFLAGS = -c
# set CC to cc for unix and gcc for cygwin 
CC = gcc
CFLAGS = -c
# set LD to f77 for unix and g77 for cygwin 
LD = g77
LDFLAGS = -o $(EXE)
OBJ = gdfprg.o gdiff.o hread.o uread.o invert.o dcell.o

$(EXE): $(OBJ)
	$(FC) $(LDFLAGS) $(OBJ)

o.c:
	$(CC) $(CFLAGS) $<

o.f:
	$(FC) $(FFLAGS) $<

uread.o: ../FORT/uread.f
	$(FC) $(FFLAGS) ../FORT/uread.f

dcell.o: ../FORT/dcell.f
	$(FC) $(FFLAGS) ../FORT/dcell.f

hread.o: ../FORT/hread.f
	$(FC) $(FFLAGS) ../FORT/hread.f

invert.o: ../FORT/invert.f
	$(FC) $(FFLAGS) ../FORT/invert.f

clean: $(OBJ)
	rm $(OBJ)

