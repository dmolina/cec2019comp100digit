CC=g++
CFLAGS=-Wall -pedantic -std=c++11 -O3
EXEC = main
DEPS = 
RM = rm
SOURCES = $(shell ls *.cpp)
OBJS = $(SOURCES:%.cpp=%.o)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJS)

%.o: %.cpp $(DEPS)
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) *.o $(EXEC)
