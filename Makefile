
all: main-Replisim2-0.8

main-Replisim2-0.8: main-Replisim2-0.8.C
	g++ -fopenmp -o main-Replisim2-0.8 main-Replisim2-0.8.C

.DONE:

