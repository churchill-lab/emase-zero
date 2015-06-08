CC=g++
CFLAGS=-I. -std=c++0x -O3 -ftree-vectorize -msse2
LDFLAGS=`python-config --ldflags`

OBJ=alignment_incidence_matrix.o python_interface.o sample_allelic_expression.o emase.o kallisto_import.o

default: emase2


emase2: $(OBJ) 
	$(CC) -o emase2 $(OBJ) $(LDFLAGS)

# the sed command here strips out a flag that doesn't make sense for C++ and may
# generate a warning on some systems.
# python-config --cflags may output -Wstrict-prototypes as a flag, which g++ may
# warn about because it is not a valid flag for C++ but it is for C
python_interface.o: python_interface.cpp
	$(CC) $(CFLAGS) `python-config --cflags | sed 's/-Wstrict-prototypes//g'` -c python_interface.cpp

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	$(RM) emase2 *.o
