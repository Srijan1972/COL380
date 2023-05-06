CC=mpic++
FLAGS=-std=c++17 -O3 -g
TEST=5
# FLAGS=-std=c++11 -O2 -g -fopenmp

sources=hybridtruss.cpp
# sources=test.cpp
objects=$(sources:.cpp=.o)

a2:$(objects)
	$(CC) $^ $(FLAGS) -o $@

run:a2
	mpirun -n 6 ./a2 --inputpath=test$(TEST)/test-input-$(TEST).gra --headerpath=test$(TEST)/test-header-$(TEST).dat --outputpath=test$(TEST)/test-output-$(TEST).txt --verbose=0 --startk=1 --endk=8

%.o: %.cpp
	$(CC) $(FLAGS) -c $<

clean:
	rm *.o a2