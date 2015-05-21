CC=g++ 
CFLAGS=-I. -std=c++0x -O3 -ftree-vectorize -msse2  
LDFLAGS=`python-config --ldflags`

OBJ=alignment_incidence_matrix.o python_interface.o sample_allelic_expression.o populase.o

default: emase2


emase2: alignment_incidence_matrix.o python_interface.o sample_allelic_expression.o populase.o
	$(CC) -o emase2 $(OBJ) $(LDFLAGS)

python_interface.o: python_interface.cpp 
	$(CC) $(CFLAGS) `python-config --cflags | sed 's/-Wstrict-prototypes//g'` -c python_interface.cpp

%.o: %.cpp
	$(CC) -c -o $@ $< $(CFLAGS) 

clean:
	$(RM) emase2 *.o
