SRC = main.cpp
OBJ = $(SRC:.cpp = .o)
PROGRAM = MonteCarlo
CC = g++
CFLAGS = -std=c++11 -g
all:$(PROGRAM)

$(PROGRAM):$(OBJ)
	$(CC) $(CFLAGS) -o $(PROGRAM) $(OBJ)

SUFFIXES: .cpp

cpp.o:
	$(CC) $(CFLAGS) -c $<
clean:
	-rm -rf *.o $(PROGRAM)
