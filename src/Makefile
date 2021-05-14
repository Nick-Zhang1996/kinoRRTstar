CC=g++
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
#CFLAGS=-I. -I$(ROOT_DIR)/../eigen -lm -g -pg
CFLAGS=-I. -lm 
DEPS = common.h

main: main.o tree.o world.o quad_optimal_control.o kino_rrt_star.o
	$(CC) -o $@ $^ $(CFLAGS)

test: test_quad_optimal_control test_tree

test_quad_optimal_control: test_quad_optimal_control.o quad_optimal_control.o quad_oc_tests.o
	$(CC) -o $@ $^ $(CFLAGS)

eigen_test: eigen.o
	$(CC) -o $@ $^ $(CFLAGS)

test_tree: test_tree.o tree.o
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.cpp %.h $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

clean: 
	rm -f *.o
.PHONY: clean
.PHONY: test
