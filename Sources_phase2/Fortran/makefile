include Make.inc

# Machines salle TP
OBJS  = main.o 
MODS  = m_subspace_iter.o

all: lib main

main:  $(MODS) $(OBJS)
	$(FC) -o $@ $^  $(LIBS)

lib:
	(cd matgen; make)

mexfile:
	mex -Dmex -largeArrayDims mex_subspace_iter_sv.c m_subspace_iter.F90 $(PLIBS)


main.o: matgen/libaux.a

clean:
	(rm -f *.o *.mod main; cd matgen; make clean)

%.o: %.f90
	$(FC) -c $<  

%.o: %.F90
	$(FC) -c $<  

%.o: %.f
	$(FC) -c $<  

