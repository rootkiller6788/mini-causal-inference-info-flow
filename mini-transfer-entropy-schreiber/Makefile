CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude -Wno-misleading-indentation -Wno-unused-function
AR = ar
TARGET = libte.a
SRCS = $(wildcard src/*.c)
OBJS = $(SRCS:.c=.o)
TEST_SRCS = $(wildcard tests/test_*.c)
TEST_BINS = $(TEST_SRCS:.c=)

.PHONY: all clean test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(AR) rcs $@ $^

src/%.o: src/%.c
	$(CC) $(CFLAGS) -c -o $@ $<

test: $(TARGET) $(TEST_BINS)
	@for t in $(TEST_BINS); do echo Running $$t...; ./$$t && echo PASS: $$t || echo FAIL: $$t; done

tests/%: tests/%.c $(TARGET)
	$(CC) $(CFLAGS) -o $@ $< -L. -lte -lm

clean:
	rm -f src/*.o $(TARGET) tests/*.o $(TEST_BINS)

