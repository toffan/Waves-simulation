# FC    = gfortran -g -fcheck=bounds,array-temps -fbacktrace
 FC    = gfortran-4.4 -O

# Machines salle TP
LIBS  = matgen/libaux.a  /mnt/n7fs/ens/tp_abuttari/acml-4.4.0/gfortran64/lib/libacml.a -lm
PLIBS  = matgen/libaux.a /mnt/n7fs/ens/tp_abuttari/acml-4.4.0/gfortran64_mp/lib/libacml_mp.a
# LIBS  = matgen/libaux.a $(LBLAS)
# PLIBS = matgen/libaux.a $(LBLAS) 

%.o: %.f90
	$(FC) -c $<  

%.o: %.F90
	$(FC) -c $<  

%.o: %.f
	$(FC) -c $<  

