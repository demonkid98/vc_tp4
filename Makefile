#######################################
# Makefile PBM                        #
#                                     #
# E.B.                                #
#######################################


PROG = gradient interest_point

all : $(PROG)

# Variables for file compilation
CC        =  gcc
CFLAGS    =  -g -Wall
CPPFLAGS  =  -DDEBUG
LDFLAGS   =  -g -lm


gradient: gradient.o Util.o
interest_point: interest_point.o Util.o

clean :
	@rm -f *.o

cleanall : clean
	@rm -f $(PROG)
