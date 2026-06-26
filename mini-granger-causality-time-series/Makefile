CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libgranger.a
# Library sources (no main() function)
LIB_SRCS = src/ts_core.c src/granger_test.c src/conditional_granger.c \
           src/transfer_entropy.c src/spectral_granger.c src/fevd.c \
           src/nonlinear_granger.c src/causality_graph.c \
           src/time_varying_granger.c src/granger_bootstrap.c \
           src/partial_granger.c
LIB_OBJS = $(patsubst src/%.c,build/%.o,$(LIB_SRCS))
# Application sources (have main(); compiled standalone)
APP_SRCS = src/granger_app_economics.c src/granger_app_neuroscience.c
APP_BINS = $(patsubst src/%.c,build/%,$(APP_SRCS))
TEST_SRC = tests/test_granger.c
TEST_BIN = build/test_granger
.PHONY: all builddir clean test examples apps demo
all: builddir $(LIB)
builddir:; @mkdir -p build
$(LIB): $(LIB_OBJS); $(AR) rcs $@ $^
build/%.o: src/%.c; $(CC) $(CFLAGS) -c $< -o $@
test: all; $(CC) $(CFLAGS) -o $(TEST_BIN) $(TEST_SRC) $(LIB_OBJS) $(LDFLAGS) && ./$(TEST_BIN)
clean:; rm -rf build
EXAMPLE_SRCS = $(wildcard examples/*.c)
EXAMPLE_BINS = $(patsubst examples/%.c,build/%,$(EXAMPLE_SRCS))
examples: all $(EXAMPLE_BINS); @echo "All examples built."
$(EXAMPLE_BINS): build/%: examples/%.c $(LIB_OBJS); $(CC) $(CFLAGS) -o $@ $< $(LIB_OBJS) $(LDFLAGS)
apps: all $(APP_BINS); @echo "All applications built."
$(APP_BINS): build/%: src/%.c $(LIB_OBJS); $(CC) $(CFLAGS) -o $@ $< $(LIB_OBJS) $(LDFLAGS)
demo: examples apps; @for e in $(EXAMPLE_BINS) $(APP_BINS); do echo "=== $$(basename $$e) ==="; ./$$e; done
