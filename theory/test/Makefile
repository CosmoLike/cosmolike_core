CXX=gcc-6

# DEFINE DEBUG TESTS
DEBUG := 1

# DEFINE GSL LIBRARY
GSL   := 1

# DEFINE NAG LIBRARY
#NAG  := 1

LDFLAGS  =
CXXFLAGS += -fmax-errors=1

ifdef DEBUG
CXXFLAGS += -O0 -g -DDEBUG -W -Wall -Wextra -Wcomment -pedantic 
else
CXXFLAGS += -O3 -march=native
endif

ifdef GSL
CXXFLAGS += -DHAVE_INLINE -DGSLLIB -I ${GSL_INC} 
LDFLAGS  += -L ${GSL_LIB} -lgsl -lgslcblas -lm
endif

ifdef NAG
CXXFLAGS += -I ${NAG_INC} -DNAGLIB
LDFLAGS += ${NAG_LIB}/libnagc_mkl.so -lpthread
endif 

all: test

main.o: main.c
	$(CXX) -c $(CXXFLAGS) main.c

test: main.o
	$(CXX) main.o $(LDFLAGS) -o $@

clean:
	rm -rf test *.o
  
