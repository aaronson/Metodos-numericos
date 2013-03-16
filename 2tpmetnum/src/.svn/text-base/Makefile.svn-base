CC=g++
CFLAGS=-g 
SOURCES= main.cpp

all: main  

%: %.cpp Makefile
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f main
