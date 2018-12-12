CC=gcc
CFLAGS=-Wall -pedantic -std=c++11 -O3
EXEC = main_orig
DEPS = 
RM = rm
SOURCES = $(shell ls *.cpp)
OBJS = $(SOURCES:%.cpp=%.o)
OBJS = cec19_func.o
LIBS = -lm

main_orig: $(OBJS) main_orig.o
	$(CC) $(CFLAGS) -o main_orig main_orig.o $(OBJS) $(LIBS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) *.o $(EXEC)
