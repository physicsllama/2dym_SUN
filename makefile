CPPOPT = g++ -O3 -std=c++11
CPPDBG = g++ -g -std=c++11

all: plot

2dym: 2dym.o 2dym-analytic.o
	g++ -O3 2dym.o 2dym-analytic.o -o 2dym

2dym.o: 2dym.cpp 2dym.h 2dym-analytic.h
	$(CPPOPT) -c 2dym.cpp

2dym-analytic.o: 2dym-analytic.cpp 2dym-analytic.h
	$(CPPOPT) -c 2dym-analytic.cpp

fermion-config: fermion-config.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c fermion-config.cpp
	g++ fermion-config.o 2dym-analytic.o 2dym.o -o fermion-config

fermion-config.pdf: fermion-config
	./fermion-config 5.0 2.0 50
	./fermion-config 9.5 2.0 50
	./fermion-config 15.0 2.0 50
	gnuplot rho.plt

heuristic-entropy: heuristic-entropy.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c heuristic-entropy.cpp
	g++ heuristic-entropy.o 2dym-analytic.o 2dym.o -o heuristic-entropy

heuristic-entropy.dat: heuristic-entropy
	./heuristic-entropy > heuristic-entropy.dat

heuristic-entropy.pdf: heuristic-entropy.dat heuristic-entropy.plt
	gnuplot heuristic-entropy.plt

free-energy: free-energy.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c free-energy.cpp
	g++ free-energy.o 2dym-analytic.o 2dym.o -o free-energy

free-energy.pdf: free-energy
	./free-energy
	gnuplot free-energy.plt

first-law: first-law.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c first-law.cpp
	g++ first-law.o 2dym-analytic.o 2dym.o -o first-law

brute-force: brute-force.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c brute-force.cpp
	g++ brute-force.o 2dym-analytic.o 2dym.o -o brute-force

variance-and-subtract: variance-and-subtract.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c variance-and-subtract.cpp
	g++ variance-and-subtract.o 2dym-analytic.o 2dym.o -o variance-and-subtract

energies: energies.cpp 2dym.o 2dym-analytic.o
	$(CPPOPT) -c energies.cpp
	g++ energies.o 2dym-analytic.o 2dym.o -o energies

first-law.dat: first-law
	./first-law

brute-force.dat: brute-force
	./brute-force

variance-and-subtract.dat: variance-and-subtract
	./variance-and-subtract

S-shannon-exp.pdf: S-shannon-exp.plt
	gnuplot S-shannon-exp.plt

S-shannon-brute.pdf: S-shannon-brute.plt
	gnuplot S-shannon-brute.plt

S-shannon-fixedN.pdf : S-shannon-fixedN.plt
	gnuplot S-shannon-fixedN.plt

S-boltzmann.pdf : S-boltzmann.plt
	gnuplot S-boltzmann.plt

rho.pdf: 2dym.dat 2dym-analytic.dat
	gnuplot rho.plt

clean:
	rm ./*.o
