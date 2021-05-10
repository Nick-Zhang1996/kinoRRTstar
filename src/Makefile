CC=g++
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
CFLAGS=-I. -lm -g -pg -fpic
DEPS = common.h

INCLUDE_PYTHON=/usr/include/python3.6

all: main kinoRRT.so

main: main.o tree.o world.o quad_optimal_control.o kino_rrt_star.o
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.cpp %.h $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.o: main.cpp world.h tree.h common.h kino_rrt_star.h
	$(CC) -c -o $@ $< $(CFLAGS)


boost_wrapper.o: boost_wrapper.cpp common.h world.h kino_rrt_star.h
	$(CC) -c -o $@ $< $(CFLAGS) -I$(INCLUDE_PYTHON) -lboost_python3

kinoRRT.so: boost_wrapper.o world.o kino_rrt_star.o quad_optimal_control.o tree.o
	$(CC) -shared -o $@ $^ $(CFLAGS) -I$(INCLUDE_PYTHON) -lboost_python3

##### tests -------------

test: test_quad_optimal_control test_tree

test_quad_optimal_control: test_quad_optimal_control.o quad_optimal_control.o quad_oc_tests.o
	$(CC) -o $@ $^ $(CFLAGS)

eigen_test: eigen.o
	$(CC) -o $@ $^ $(CFLAGS)

test_tree: test_tree.o tree.o
	$(CC) -o $@ $^ $(CFLAGS)

clean: 
	rm *.o
	rm *.so
.PHONY: clean
.PHONY: test
