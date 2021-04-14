CC=g++
CFLAGS=-I. -lm
DEPS = 

test: test_quad_optimal_control
%.o: %.cpp %.h $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)
test_quad_optimal_control: test_quad_optimal_control.o quad_optimal_control.o
	$(CC) -o $@ $^ $(CFLAGS)
clean: 
	rm -f *.o
.PHONY: clean
.PHONY: test
