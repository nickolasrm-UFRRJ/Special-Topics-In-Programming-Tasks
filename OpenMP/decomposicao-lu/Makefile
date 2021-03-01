CPUC      = gcc
CFLAGS    = -O3 --std=c11 -m64 -fopenmp 	 #     -Wunused-variable #-Wno-conversion-null -Wdelete-non-virtual-dtor
DEFS      =
INCLUDES  =
LIBDIR   = -L/usr/lib
LIBS     = -lm
LINK     =  $(LIBDIR) $(LIBS)
CPU_COMPILE_C  = $(CPUC) $(DEFS) $(INCLUDES) $(CFLAGS)

all:
	$(CPU_COMPILE_C) main-lu.c  $(LINK)   -o lu

backup:
	$(CPU_COMPILE_C) backup.c  $(LINK)   -o lu-backup

clean:
	rm lu
