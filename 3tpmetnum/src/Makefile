CC=g++
CFLAGS=-g 
SOURCES= uniforme.cpp chord.cpp centrip.cpp

all: uniforme centrip chord  

%: %.cpp Makefile
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f uniforme centrip chord
