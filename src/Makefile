CC=g++
CFLAGS=-I. -std=c++0x -O3 -ftree-vectorize -msse2

OBJ=alignment_incidence_matrix.o sample_allelic_expression.o alignment_import.o emase.o

default: emase-zero

emase-zero: $(OBJ)
	$(CC) -o emase-zero $(OBJ) -lz

emase-dump: emase_dump.o
	$(CC) -o emase-dump emase_dump.o -lz 

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	$(RM) emase-zero emase-dump *.o
