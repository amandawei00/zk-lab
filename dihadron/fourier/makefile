# PC running Linux
LINUX_ON_PC = 1

# Something else
#LINUX_ON_PC = 0

# Use Fortran calling conventions
LINK_F = 1

# Use C calling conventions
#LINK_F = 0

#############################################################################

GCCFLAGS = -Wall
G77FLAGS = -funderscoring -fno-second-underscore -Wall
OPTFLAGS = -O3 -ffast-math
#OPTFLAGS = -O1
GSL =
OBJECTS_C = example_c.o fourier.o
OBJECTS_F = example_f.o libFcompat.o fourier.o

LIBS_C = -lgsl -lgslcblas -lm
LIBS_F = -lgsl -lgslcblas -lm

############################################################################

ifeq ($(LINUX_ON_PC),1)
	PC = -D X86
else
	PC = 
endif

ifeq ($(LINK_F),1)
	FF = -D FORTRAN
else
	FF = 
endif

############################################################################

example_f: $(OBJECTS_F)
	gfortran $(GCCFLAGS) -o example_f $(OBJECTS_F) $(LIBS_F)

example_f.o: example.f
	gfortran $(G77FLAGS) $(OPTFLAGS) -c -o example_f.o example.f

############## for calling the library from C ############################# 

example_c: $(OBJECTS_C)
	gcc $(GCCFLAGS) -o example_c $(OBJECTS_C) $(LIBS_C)

example_c.o: example.c
	gcc $(GCCFLAGS) $(OPTFLAGS) -c -o example_c.o example.c

##########################################################################

libFcompat.o: libFcompat.c
	gcc $(GCCFLAGS) $(FF) -c -o libFcompat.o libFcompat.c

fourier.o: fourier.c
	gcc $(GCCFLAGS) $(OPTFLAGS) $(GSL) $(PC) $(FF) -c -o fourier.o fourier.c 
###########################################################################

clean:
	rm -f *.o example_f example_c



