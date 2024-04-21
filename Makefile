CFLAGS= -Wall -Wextra -g -march=native -fopenmp -O3

BINS=find_lbox heuristic_implem L32x3 L32x4

all: $(BINS)

heuristic_implem: heuristic_implem.c
	clang -Wall -Wextra -g -march=native -O3 $< -o $@

clean:
	rm -f *.o $(BINS)

.PHONY: all clean
