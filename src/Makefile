CC=g++
ROOT_DIR:=$(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
CFLAGS=-I. -lm -fpic -Wall -g
DEPS = common.h
LIBRARY_PATH = usr/lib/x86_64-linux-gnu:$LIBRARY_PATH

INCLUDE_PYTHON=/usr/include/python3.8

all: main kinoRRT.so

main: main.o tree.o world.o quad_optimal_control.o kino_rrt_star.o
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.cpp %.h $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.o: main.cpp world.h tree.h common.h kino_rrt_star.h
	$(CC) -c -o $@ $< $(CFLAGS)


boost_wrapper.o: boost_wrapper.cpp common.h world.h kino_rrt_star.h
	$(CC) -c -o $@ $< $(CFLAGS) -I$(INCLUDE_PYTHON) -lboost_python38

kinoRRT.so: boost_wrapper.o world.o kino_rrt_star.o quad_optimal_control.o tree.o
	$(CC) -shared -o $@ $^ $(CFLAGS) -I$(INCLUDE_PYTHON) -lboost_python38

##### tests -------------

test: test_quad_optimal_control test_tree

test_quad_optimal_control: test_quad_optimal_control.o quad_optimal_control.o quad_oc_tests.o
	$(CC) -o $@ $^ $(CFLAGS)

eigen_test: eigen.o
	$(CC) -o $@ $^ $(CFLAGS)

test_tree: test_tree.o tree.o
	$(CC) -o $@ $^ $(CFLAGS)

clean: 
	rm -f *.o
	rm -f *.so
.PHONY: clean
.PHONY: test
