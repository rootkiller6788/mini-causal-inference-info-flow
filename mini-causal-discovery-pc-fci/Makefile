CC=gcc
CFLAGS=-std=c11 -Wall -Wextra -O2 -Iinclude
LIB=libcdf.a
SRCS=$(wildcard src/*.c)
OBJS=$(SRCS:.c=.o)
all: $(LIB)
$(LIB): $(OBJS)
	ar rcs $@ $^
src/%.o: src/%.c include/*.h
	$(CC) $(CFLAGS) -c $< -o $@
test: $(LIB)
	for t in tests/test_*.c; do $(CC) $(CFLAGS) $$t -L. -lcdf -lm -o $$t.out && ./$$t.out; rm -f $$t.out; done

demo: $(LIB)
	for d in demos/demo_*.c; do $(CC) $(CFLAGS) $$d -L. -lcdf -lm -o $$d.out && echo "--- Running $$d ---" && ./$$d.out && echo ""; rm -f $$d.out; done

bench: $(LIB)
	for b in benches/bench_*.c; do $(CC) $(CFLAGS) $$b -L. -lcdf -lm -o $$b.out && echo "--- Running $$b ---" && ./$$b.out && echo ""; rm -f $$b.out; done

clean:
	rm -f src/*.o $(LIB) tests/*.out demos/*.out benches/*.out
