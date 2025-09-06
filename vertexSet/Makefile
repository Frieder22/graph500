GCC = mpicc
CFLAGS = -O3 -I.

target = main
SOURCES = vertexSet.c vertexSetIterator.c
HEADERS = vertexSet.h bitmap_custom.h vertexSetIterator.h testStructure.h
 
build = ../build

main: main.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/$(target).o $(target).c $(SOURCES)

run: main 
	mpirun -np 4 $(build)/$(target).o

test: test.c $(SOURCES) $(HEADERS)
	$(GCC) $(CFLAGS) -o $(build)/test.o test.c $(SOURCES)

test_it: test
	mpirun -np 1 ./$(build)/test.o 0

test_it_shifted: test
	mpirun -np 3 ./$(build)/test.o 0

test_red_exactHalfing: test
	mpirun -np 4 ./$(build)/test.o 1

test_red_approxHalfing: test
	mpirun -np 3 ./$(build)/test.o 2

test_red_ringComm: test
	mpirun -np 4 ./$(build)/test.o 3

test_red_dense: test
	mpirun -np 4 ./$(build)/test.o 4

test_vertexset: test
	mpirun -np 1 ./$(build)/test.o 5

perf: performance.c
	$(GCC) $(CFLAGS) -o $(build)/perf.o performance.c $(SOURCES)

perf_red_exactHalfing: perf
	mpirun -np 4 ./$(build)/perf.o 0

perf_red_approxHalfing: perf
	mpirun -np 4 ./$(build)/perf.o 1

perf_red_ringComm: perf
	mpirun -np 4 ./$(build)/perf.o 2

perf_red_dense: perf
	mpirun -np 4 ./$(build)/perf.o 3

clean:
	rm -rf ../build