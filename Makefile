#Rafael Marinho
#fdtd Makefile

.PHONY: all clean

CORE = $(basename $(FILE))
OUTPUT = fdtd

all:	fdtd

fdtd:	$(CORE).o
	g++ $(CORE).o -o $(OUTPUT)

$(CORE).o:	
	g++ $(CORE).cpp -larmadillo -lstdc++ -lm -lhdf5 -o $(CORE).o

clean:
	rm -rf *.o 
	rm -rf $(OUTPUT)

