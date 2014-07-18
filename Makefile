CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#  -Ofast can randomly cause nan values.
#CFLAGS = -std=gnu++11 -g -Wall
CFLAGS = -std=gnu++11 -Ofast -Wall
#CFLAGS = -std=gnu++11 -O3 -Wall

LIBS = -lfftw3_threads -lfftw3 -lm -lpthread -lboost_filesystem -lboost_iostreams -lboost_system

SRCS = ExpTimeDiff.cpp etdrk4b.cpp model_g.cpp

OBJS = $(SRCS:.cpp=.o)

MAIN = model_g

.PHONY: depend clean

all: $(MAIN)
	@echo $(MAIN) has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS) $(LIBS)

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $^

# DO NOT DELETE THIS LINE -- make depend needs it
