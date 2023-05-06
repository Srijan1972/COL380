CC=nvcc
FLAGS=-arch=sm_35 -std=c++17 -O3 -g
TEST=2
N=100
M=4
KA=200
KB=200

sources=main.cu
objects=$(sources:.cu=.o)

exec:$(objects)
	$(CC) $^ $(FLAGS) -o $@

gen:
	python gen-input.py -n $(N) -m $(M) -z $(KA) -o Atest$(TEST)
	python gen-input.py -n $(N) -m $(M) -z $(KB) -o Btest$(TEST)
	python gen-output.py -f Atest$(TEST) Btest$(TEST) -e 2 -o Ctest$(TEST)_match

run:exec
	./exec Atest$(TEST) Btest$(TEST) Ctest$(TEST)
	# python checker.py -f Ctest$(TEST)_match Ctest$(TEST) -e 4

%.o: %.cu
	$(CC) $(FLAGS) -c $<

clean:
	rm *.o exec