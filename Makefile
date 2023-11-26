CC = mpicc
CFLAGS = -Wall -lm
EXEC = heatsink_para

NB_PROC ?= 1
DIM ?= 3

all: $(EXEC)

$(EXEC): heatsink_v1.o
	$(CC) -o $(EXEC) heatsink_v1.o $(CFLAGS)

heatsink_v1.o: heatsink_v1.c
	$(CC) -c heatsink_v1.c $(CFLAGS)

run: $(EXEC)
	mpiexec --n $(NB_PROC) --hostfile $$OAR_NODEFILE $(EXEC) $(DIM)

clean:
	rm -f $(EXEC) heatsink_v1.o