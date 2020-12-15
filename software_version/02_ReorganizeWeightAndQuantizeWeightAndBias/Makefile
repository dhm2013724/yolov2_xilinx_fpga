
#DEFS = -D REORG_TEST
#all: test

test:
	g++ -O3 -w -o test_layers main.cpp -I . -lm -D REORG_TEST
gen:
	g++ -O3 -w -o test_layers main.cpp -I . -lm -D REORG_GEN
quanti:
	g++ -O3 -w -o test_layers main.cpp -I . -lm -D REORG_TEST -D QUANTI
clean:
	rm test_layers

