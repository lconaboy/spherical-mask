Description
===========

This program combines ginnungagap files to produce ramses-readable
(i.e. in the grafic format) zoom intial conditions (ICs) files. This 
is based heavily on the routines in mesh.hh in MUSIC (Hahn & Abel, 2011). 


Usage
=====

You can easily run by using the provided Makefile (type `make`) and
run.sh. If these don't work, the process is as follows:

    gfortran cubic_mask.f95 -o cubic_mask    
    mkdir level_***
    ./cubic_mask paramfile.param

First compile cubic_mask.f95, then make the relevant directories, one
for each level of refinment from levelmin to levelmax (the program
will write the IC files into these). An example .param file is
included; to avoid it being overwritten when changes are pulled into
the repo, copy it, rename it to cubic_mask.param and fill in the
relevant fields to use it. The following fields need to be specified
in the namelist:

    - levelmin: base level of refinement
    - levelmax: highest level of refinement
    - ix:       x-corner of refinement mask
    - iy:       y-corner of refinement mask
    - iz:       z-corner of refinement mask
    - nx:       x-extent of refinement mask
    - ny:       y-extent of refinement mask
    - nz:       z-extent of refinement mask
    - npad:     number of padding cells
    - path:     path to ginnungagap files

Note that ix, iy, iz, nx, ny and nz should all be in units of
finest-level grid cells.

Most of a ramses namelist will be written, but you will need to add to
the OUTPUT_PARAMS and HYDRO_PARAMS blocks. You may also need to adjust
some other parameters, e.g. ngridmax and npartmax, or the full path to
the initfiles.

Acknowledgements
================

Parts originally written by Sergey Pilipenko. Parts of the refinement
mask generation were based on the MUSIC code, described in Hahn & Abel
(2011).
