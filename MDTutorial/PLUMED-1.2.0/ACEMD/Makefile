CPPFLAGS=-I./
SRC= *.c

all: $(SRC) 
	gcc  $(CPPFLAGS) -O3 -DACEMD -fPIC --shared -o plugin.so $(SRC)
