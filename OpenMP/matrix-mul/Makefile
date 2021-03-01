################################################################################
#
#  Universidade Federal Rural do Rio de Janeiro
#  por: Marcelo Zamith - mzamith@ufrrj.br
#  Multiplicação de matrizes
#
################################################################################
EXEFILE      = mat-mul
CPUCC     = g++
CPPFLAGS  = -mavx -std=c++11 -m64 -fpermissive -fopenmp
DEFS      = $(APPLICATION)  $(VERSION)
INCLUDES  = -I. 
LIBDIR    = -L/usr/lib
LIBS     =  -lm
LINK     =  $(LIBDIR) $(LIBS)
DEFS     =  

ifeq ($(TARGET), release)
  DEFS      +=
  CPPFLAGS  +=  -O3
else
  DEFS     += -D_DEBUG
  CPPFLAGS += -g -Wall -O0
endif

CPU_COMPILE = $(CPUCC) $(DEFS) $(INCLUDES) $(CPPFLAGS)


all:	main
	$(CPU_COMPILE)	main.o $(LINK) -o $(EXEFILE) $(PAPILIBS)

main:
	$(CPU_COMPILE) -c main.cpp


clean:
	rm *.o $(EXEFILE)
