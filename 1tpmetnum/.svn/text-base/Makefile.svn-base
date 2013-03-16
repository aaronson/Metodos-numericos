CC=g++
CFLAGS=-g 
SOURCES= main.cpp mainpara2a.cpp

all: main mainpara2a 

%: %.cpp Makefile
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f main mainpara2a
