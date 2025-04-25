ffiles = $(wildcard *.f90)
cfiles = $(wildcard *.c)
OBJ= $(ffiles:.f90=.o) $(cfiles:.c=.o) 


ifeq ($(arch),gfortran)
    F90 = gfortran -O3 -fno-second-underscore 
    CC  = gcc -O3 -DDRIVER 
    LINK = gfortran 
endif

ifeq ($(arch),g95)
    F90 = g95 -O3 -fno-second-underscore     
    CC  = gcc -O3 -DDRIVER     
    LINK = g95
endif

ifeq ($(arch),gfortran_64)
    F90 = gfortran -m64 -O3 -fno-second-underscore     
    CC  = gcc -O3 -DDRIVER -m64     
    LINK = gfortran -m64
endif

ifeq ($(arch),intel)
    F90 = ifort -O3 
    CC  = icc -O3 -DDRIVER
    LINK = ifort 
endif

ifeq ($(arch),'')
  error:
        @echo 'ERROR: Please specify the right architecture (see Makefile for choice)'
endif


EXEC = driver 

all:$(EXEC)

g95:
	@make arch=g95

gfortran:
	@make arch=gfortran

gfortran_64:
	@make arch=gfortran_64

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
