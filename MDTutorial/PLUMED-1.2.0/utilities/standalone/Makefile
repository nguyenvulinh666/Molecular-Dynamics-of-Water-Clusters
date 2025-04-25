ffiles = $(wildcard *.f90)
cfiles = $(wildcard *.c)
OBJ= $(ffiles:.f90=.o) $(cfiles:.c=.o) 


ifeq ($(arch),gnu)
    F90 = g95 -O3 -fno-second-underscore 
    CC  = gcc -O3 -DSTANDALONE
    LINK = g95
endif

ifeq ($(arch),dbg)
    F90 = g95 -g -fno-second-underscore 
    CC  = gcc -g -DSTANDALONE
    LINK = g95
endif

ifeq ($(arch),intel)
    F90 = ifort -O3 
    CC  = icc -O3 -DSTANDALONE
    LINK = ifort 
endif

ifeq ($(arch),'')
  error:
        @echo 'ERROR: Please specify the right architecture (see Makefile for choice)'
endif


EXEC = plumed_standalone 

all:$(EXEC)

gnu:
	@make arch=gnu
dbg:
	@make arch=dbg
intel:
	@make arch=intel

$(EXEC): $(OBJ)
	$(LINK) $(CFlAGS) $(OBJ) -o $(EXEC)  $(LIBS)

%.o: %.c
	$(CC) $(CFLAGS) -c  -o $@ $*.c

%.o: %.f90
	$(F90) $(PFLAGS) -c  -o $@ $^

clean:
	@rm -rf *.o  
