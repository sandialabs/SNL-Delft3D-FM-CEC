#-------------------------------------------------------------------------------
#   Makefile for delftflow
#
#   Irv.Elshoff@deltares.nl
#   Adri.Mourits@deltares.nl
#   Fedor.Baart@deltares.nl
#   11 jan 2007
#
#   Copyright (C)  Stichting Deltares, 2011-2013.
#-------------------------------------------------------------------------------

include $(top_srcdir)/common.am
noinst_LTLIBRARIES = libDataGDP.la

#-------------------------------------------------------------------------------
#Includes
FINCLUDES += \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/include

FMODULES += \
        -I$(top_builddir)/utils_gpl/morphology/packages/morphology_data/src \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/src/basics \
	-I$(top_builddir)/engines_gpl/flow2d3d/packages/data/src/parallel_mpi \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common/src \
	-I$(top_builddir)/utils_lgpl/deltares_common/packages/deltares_common_mpi/src \
	-I$(top_builddir)/utils_lgpl/ec_module/packages/ec_module/src \
	-I$(top_builddir)/utils_gpl/morphology/packages/morphology_kernel/src

#-------------------------------------------------------------------------------
#Sources
libDataGDP_la_SOURCES = \
         gdp_for_dll.f90 \
         clradv2d.f90 \
         clrbedformpar.f90 \
         clrdredge.f90 \
         clrflwpar.f90 \
         clrmassbal.f90 \
         clrtrachy.f90 \
         clrwaqpar.f90 \
         gdp_alloc.f90 \
         gdp_alloc_arrays.f90 \
         gdp_dealloc.f90 \
         initadv2d.f90 \
         initarrays.f90 \
         initbcdat.f90 \
         initbedformpar.f90 \
         initbubble.f90 \
         initcline.f90 \
         initcoup.f90 \
         initculver.f90 \
         initdischarge.f90 \
         initveg3d.f90 \
         initdredge.f90 \
         initf0isf1.f90 \
         initflwpar.f90 \
         initfmtbcc.f90 \
         initfmtbct.f90 \
         initfmtdis.f90 \
         initfourier.f90 \
         initheat.f90 \
         initincbc.f90 \
         initincbcc.f90 \
         initincwav.f90 \
         initkeywtd.f90 \
         initmassbal.f90 \
         initpostpr.f90 \
         initprocs.f90 \
         initrestart.f90 \
         initrtc.f90 \
         initsafe.f90 \
         initscour.f90 \
         initsdu.f90 \
         initsnel.f90 \
         initsobek.f90 \
         initstations.f90 \
         inittimers.f90 \
         inittrachy.f90 \
         inittrisol.f90 \
         initu_ppr.f90 \
         initupdbcc.f90 \
         initupdbct.f90 \
         initupddis.f90 \
         initwaqpar.f90 \
         initwrirst.f90 \
         initwrline.f90 \
         initz_initcg.f90 \
         initzmodel.f90 \
         iofiles_alloc.f90 \
         iofiles_dealloc.f90


# flow2d3dinclude_HEADERS=*.mod
#-------------------------------------------------------------------------------
#Libraries

