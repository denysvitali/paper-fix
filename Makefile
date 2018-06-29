CFLAGS = `pkg-config --cflags opencv`
LIBS = `pkg-config --libs opencv`

main:
	g++ -Wall $(CFLAGS) $(LIBS) main.cpp -o main.out

.PHONY: main

run: main
	./main.out
