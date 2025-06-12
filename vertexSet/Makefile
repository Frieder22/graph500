GCC = mpicc
CFLAGS = -O3 -I.

target = main
SOURCES = vertexSet.c vertexSetIterator.c
HEADERS = vertexSet.h bitmap_custom.h vertexSetIterator.h
 
build = ../build

main: main.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/$(target).o $(target).c $(SOURCES)

run: main 
	mpirun -np 4 $(build)/$(target).o

test: test.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/test.o test.c $(SOURCES)

test_it: test
	mpirun -np 1 ./$(build)/test.o 0

test_red_exactHalfing: test
	mpirun -np 4 ./$(build)/test.o 1

test_red_approxHalfing: test
	mpirun -np 4 ./$(build)/test.o 2

perf: performance.c
	$(GCC) $(CFLAGS) -o $(build)/perf.o performance.c $(SOURCES)

perf_red_exactHalfing: perf
	mpirun -np 4 ./$(build)/perf.o 0

perf_red_approxHalfing: perf
	mpirun -np 4 ./$(build)/perf.o 1

clean:
	rm -rf ../build