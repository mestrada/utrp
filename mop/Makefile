#OBJS = main.o   matrices.o solution.o Graph.o
OBJS = main.o matrix.o solution.o Graph.o
CC = g++ -std=c++11 -stdlib=libc++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)     

all: $(OBJS)  
	$(CC)   $(LFLAGS)   $(OBJS) -o  urtp

clean:
	\rm *.o
