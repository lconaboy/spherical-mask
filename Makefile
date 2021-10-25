spherical_mask:
	gfortran spherical_mask.f95 -O2 -o sm

debug:
	gfortran spherical_mask.f95 -fbounds-check -ffpe-trap='invalid,zero' -fcheck='all' -Wall -o sm

.PHONY : clean
clean :
	-rm -f sm
