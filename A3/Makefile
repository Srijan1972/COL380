CC=mpic++
FLAGS=-std=c++17 -O3 -g -fopenmp
TEST=1

sources=hybridtruss.cpp
objects=$(sources:.cpp=.o)

a3:$(objects)
	$(CC) $^ $(FLAGS) -o $@

run:a3
	mpirun -n 6 ./a3 --taskid=2 --inputpath=test$(TEST)/test-input-$(TEST).gra --headerpath=test$(TEST)/test-header-$(TEST).dat --outputpath=test$(TEST)/test-output-$(TEST).txt --verbose=0 --startk=2 --endk=2 --p=2

%.o: %.cpp
	$(CC) $(FLAGS) -c $<

clean:
	rm *.o a3
