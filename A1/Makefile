CC=g++
FLAGS=-std=c++11 -O3 -fopenmp -L. -lrary
# FLAGS=-std=c++11 -O2 -g -fopenmp

sources=main.cpp multiply.cpp
objects=$(sources:.cpp=.o)

exec:$(objects)
	$(CC) $^ $(FLAGS) -o $@

run:exec
	./exec input out

%.o: %.cpp library.hpp multiply.hpp
	$(CC) $(FLAGS) -c $<

clean:
	rm *.o exec
