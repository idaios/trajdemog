CC=gcc
CFLAGS=
DEPS=freqints.h
LIBS=-lm

%.o: %.c  $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: trajdemog stepftn2

trajdemog: popsize.o binomial.o rand1.o trajdemog.o
	gcc -o trajdemog popsize.o binomial.o rand1.o trajdemog.o $(LIBS)

stepftn2: stepftn2.o
	gcc -o stepftn2 stepftn2.o $(LIBS)

clean:
	rm -f *.o trajdemog stepftn2





