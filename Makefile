CFLAGS = -std=c++11 -pedantic -Wall -I../msieve/include
LDFLAGS = -L../msieve -lgmpxx -lgmp -lmsieve -lz -lecm -lpthread -lpng -ldl

all: exp

exp: exp.cpp experiment.cpp experiment.h exp_container.cpp exp_container.h fretriever.cpp fretriever.h helpers.cpp helpers.h heap.h factorize.cpp factorize.h
	g++ $(CFLAGS) -march=native -O3 -o exp exp.cpp experiment.cpp exp_container.cpp fretriever.cpp helpers.cpp factorize.cpp $(LDFLAGS)

exp_debug: exp.cpp experiment.cpp experiment.h exp_container.cpp exp_container.h fretriever.cpp fretriever.h helpers.cpp helpers.h heap.h factorize.cpp factorize.h
	g++ $(CFLAGS) -g -pg -o exp exp.cpp experiment.cpp exp_container.cpp fretriever.cpp helpers.cpp factorize.cpp $(LDFLAGS)

clean: 
	rm exp
