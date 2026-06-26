CC = gcc
AR = ar
CFLAGS = -std=c11 -Wall -Wextra -O2 -Iinclude
LDFLAGS = -lm
LIB = build/libcausal.a

# Source files for the library (exclude _app*.c and examples)
LIB_SRCS = src/dag_representation.c src/scm_model.c src/causal_identification.c src/adjustment_sets.c src/effect_estimators.c src/sensitivity_analysis.c src/mediation_analysis.c
OBJS = $(patsubst src/%.c,build/%.o,$(LIB_SRCS))

# Application binaries
APP_SRCS = $(wildcard src/*_app*.c)
APP_BINS = $(patsubst src/%.c,build/%,$(APP_SRCS))

# Test binaries
TESTS = $(wildcard tests/test_*.c)
TBINS = $(patsubst tests/%.c,build/%,$(TESTS))

all: builddir $(LIB) apps

builddir:
	@mkdir -p build

apps: $(APP_BINS)

$(LIB): $(OBJS)
	$(AR) rcs $@ $^

build/%.o: src/%.c
	$(CC) $(CFLAGS) -c $< -o $@

build/%_app1: src/%_app1.c builddir $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcausal $(LDFLAGS)

build/%_app2: src/%_app2.c builddir $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcausal $(LDFLAGS)

test: builddir $(LIB) $(TBINS)
	@passed=0; failed=0; for t in $(TBINS); do echo "Running $$t..."; if ./$$t; then passed=$$((passed+1)); else failed=$$((failed+1)); fi; done; echo "=== $$passed passed, $$failed failed ==="; if [ $$failed -gt 0 ]; then exit 1; fi

build/test_%: tests/test_%.c builddir $(LIB)
	$(CC) $(CFLAGS) -o $@ $< -Lbuild -lcausal $(LDFLAGS)

clean:
	rm -rf build
