#Include commands common to every Makefile.am
include $(top_srcdir)/config/makefile_templates/demo_drivers

# Name of executable
check_PROGRAMS=airy_cantilever airy_cantilever2_adapt \
airy_cantilever2_noadapt spherical_solid spherical_solid_compressed \
spherical_solid_imposedDisp \
solid_contact spherical_solid_w_contact \
spherical_solid_w_contact_static \
spherical_solid_w_contact_adaptive

#-------------------------------------------------

# Sources for executable
solid_contact_SOURCES = solid_contact.cc 

# Required libraries:
# $(FLIBS) is included in case the solver involves fortran sources.
solid_contact_LDADD = -L@libdir@  -lsolid -lconstitutive -lgeneric \
				$(EXTERNAL_LIBS) $(FLIBS)


#-------------------------------------------------

#-----------------------------------------------------------------------

# Sources for executable
airy_cantilever_SOURCES = airy_cantilever.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------

# Sources for executable
airy_cantilever2_adapt_SOURCES = airy_cantilever2.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever2_adapt_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

airy_cantilever2_adapt_CXXFLAGS=-DREFINE

#-----------------------------------------------------------------------

# Sources for executable
airy_cantilever2_noadapt_SOURCES = airy_cantilever2.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
airy_cantilever2_noadapt_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------

# Sources for executable
spherical_solid_SOURCES = spherical_solid.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
spherical_solid_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------

# Sources for executable
spherical_solid_compressed_SOURCES = spherical_solid_compressed.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
spherical_solid_compressed_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------


#-----------------------------------------------------------------------

# Sources for executable
spherical_solid_imposedDisp_SOURCES = spherical_solid_imposedDisp.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
spherical_solid_imposedDisp_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

# Sources for executable
spherical_solid_w_contact_SOURCES = spherical_solid_w_contact.cc

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
spherical_solid_w_contact_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)

#-----------------------------------------------------------------------

#-----------------------------------------------------------------------

# Sources for executable
spherical_solid_w_contact_static_SOURCES = spherical_solid_w_contact.cc

spherical_solid_w_contact_static_LDFLAGS =  -static -static-libgcc -static-libstdc++ -static-libgfortran

# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
spherical_solid_w_contact_static_LDADD = -L@libdir@ -l:libsolid.a -l:libconstitutive.a \
                        -l:libgeneric.a -l:liboomph_arpack.a -l:liboomph_superlu_4.3.a\
			-l:liboomph_lapack.a -l:liboomph_flapack.a -l:liboomph_blas.a \
			-l:libgfortran.a \
			-l:libquadmath.a  \
			$(EXTERNAL_LIBS) $(FLIBS)
#-l:libm.a -l:libc.a  

#spherical_solid_w_contact_static_CFLAGS=-static -static-libgcc -static-libstdc++ -static-libgfortran
spherical_solid_w_contact_static_CXXFLAGS=-DHTCONDOR
#-----------------------------------------------------------------------


#-----------------------------------------------------------------------

# Sources for executable
spherical_solid_w_contact_adaptive_SOURCES = spherical_solid_w_contact.cc



# Required libraries:
# $(FLIBS) is included in case the  solver involves fortran sources.
spherical_solid_w_contact_adaptive_LDADD = -L@libdir@ -lsolid -lconstitutive \
                        -lgeneric $(EXTERNAL_LIBS) $(FLIBS)
                        
spherical_solid_w_contact_adaptive_CXXFLAGS=-DREFINE

#-----------------------------------------------------------------------

