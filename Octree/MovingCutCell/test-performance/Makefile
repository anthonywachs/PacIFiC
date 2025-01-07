all: 	starting.tst \
	starting-bmr.tst \
	starting-uniform.tst \
	starting-moving.tst \
	starting-moving-bmr.tst \
	starting-moving-uniform.tst

include $(BASILISK)/Makefile.defs

# the default CFLAGS are set in $(BASILISK)/config
CFLAGS += -O2 -DMTRACE=2 -DTRACE=2 -DDISPLAY=1

# The following tests have special dependencies/compilation requirements

starting-bmr.tst: CFLAGS += -DBMR=1
starting-bmr.c: starting.c
	ln -s starting.c starting-bmr.c

starting-uniform.tst: CFLAGS += -DUNIFORM=1
starting-uniform.c: starting.c
	ln -s starting.c starting-uniform.c

starting-moving-bmr.tst: CFLAGS += -DBMR=1
starting-moving-bmr.c: starting-moving.c
	ln -s starting-moving.c starting-moving-bmr.c

starting-moving-uniform.tst: CFLAGS += -DUNIFORM=1
starting-moving-uniform.c: starting-moving.c
	ln -s starting-moving.c starting-moving-uniform.c


