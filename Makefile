all: exp

old_util: util.c heap.c heap.h
	g++ -march=native -Wall -Imsieve/msieve-1.50/include -Lmsieve/msieve-1.50 -O3 -o util util.c heap.c -lm -lgmp -lmsieve -lz -lecm -lpthread

exp: exp.cpp experiment.cpp experiment.h exp_container.cpp exp_container.h fretriever.cpp fretriever.h helpers.cpp helpers.h heap.h factorize.cpp factorize.h
	g++ -std=c++11 -pedantic -march=native -Wall -Imsieve/include -Lmsieve -O3 -o exp exp.cpp experiment.cpp exp_container.cpp fretriever.cpp helpers.cpp factorize.cpp -lgmpxx -lgmp -lmsieve -lz -lecm -lpthread -lpng -ldl

exp_debug: exp.cpp experiment.cpp experiment.h exp_container.cpp exp_container.h fretriever.cpp fretriever.h helpers.cpp helpers.h heap.h factorize.cpp factorize.h
	g++ -std=c++11 -pedantic -Wall -Imsieve/include -Lmsieve -g -pg -o exp exp.cpp experiment.cpp exp_container.cpp fretriever.cpp helpers.cpp factorize.cpp -lgmpxx -lgmp -lmsieve -lz -lecm -lpthread -lpng -ldl

clean: 
	rm exp
