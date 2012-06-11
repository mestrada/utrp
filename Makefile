OBJS = main.o   matrices.o solution.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)     

all: $(OBJS)  
	$(CC)   $(LFLAGS)   $(OBJS) -o  urtp

clean:
	\rm *.o