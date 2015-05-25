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
	/usr/local/matlab/bin/mex -c mex_subspace_iter_sv.c
	/usr/local/matlab/bin/mex -Dmex -largeArrayDim -c mex_subspace_iter_sv.o m_subspace_iter.F90 $(PLIBS)


main.o: matgen/libaux.a

clean:
	(rm -f *.o *.mod main; cd matgen; make clean)

%.o: %.f90
	$(FC) -c $<  

%.o: %.F90
	$(FC) -c $<  

%.o: %.f
	$(FC) -c $<  

