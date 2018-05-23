CC=g++
CFLAGS=-I. -std=c++0x -O3 -ftree-vectorize -msse2

OBJ=alignment_incidence_matrix.o sample_allelic_expression.o alignment_import.o emase.o

default: emase-zero

emase-zero: $(OBJ)
	$(CC) -lz -o emase-zero $(OBJ)

emase-dump: emase_bin.o
	$(CC) -lz -o emase-dump emase_bin.o

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	$(RM) emase-zero *.o
