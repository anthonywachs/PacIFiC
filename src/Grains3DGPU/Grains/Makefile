# Compilation & linking flags
CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
NVCC = $(PACIFIC_NVCC)

include ../mk/nvcc_link_flags.mk

NVCC_ARCHFLAGS += -arch=sm_${PACIFIC_NVCC_ARCH}
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS) $(PACIFIC_CURAND_CFLAGS)
LDLIBS   += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_CURAND_LFLAGS) $(PACIFIC_CUDART_LFLAGS)
NVCC_LDLIBS := $(call to_nvcc_ldflags,$(LDLIBS))

# Target directories
OBJDIR := $(LIBGRAINS_GPU_BUILDDIR)/obj
LIBDIR := $(LIBGRAINS_GPU_BUILDDIR)/lib
SHAREDIR := $(LIBGRAINS_GPU_BUILDDIR)/share/Grains3DGPU
INCDIR := $(LIBGRAINS_GPU_BUILDDIR)/include

LIBSO := libGrains3DGPU.so
LIBSO_WP := $(LIBDIR)/$(LIBSO)
OBJLIB_MANIFEST := $(LIBDIR)/libGrains3DGPU_objects.mk
RPATH_FLAGS := -Xlinker -rpath -Xlinker $(LIBDIR)

ifeq ($(strip $(SHAREDIR)),)
SHAREDIR := $(LIBDIR)/../share
endif

ALLINC = 

# Source and include directories
SRCDIRBASE := Base/src
INCDIRBASE := Base/include
INCDIRS += $(INCDIRBASE)
ALLINC += $(wildcard $(INCDIRBASE)/*.hh)

SRCDIRCOLLISION := CollisionDetection/src
INCDIRCOLLISION := CollisionDetection/include
INCDIRS += $(INCDIRCOLLISION)
ALLINC += $(wildcard $(INCDIRCOLLISION)/*.hh)

SRCDIRCOMPO := Component/src
INCDIRCOMPO := Component/include
INCDIRS += $(INCDIRCOMPO)
ALLINC += $(wildcard $(INCDIRCOMPO)/*.hh)

SRCDIRCONTACT := ContactForce/src
INCDIRCONTACT := ContactForce/include
INCDIRS += $(INCDIRCONTACT)
ALLINC += $(wildcard $(INCDIRCONTACT)/*.hh)

SRCDIRGEOM := Geometry/src
INCDIRGEOM := Geometry/include
INCDIRS += $(INCDIRGEOM)
ALLINC += $(wildcard $(INCDIRGEOM)/*.hh)

SRCDIRPP := PostProcessing/src
INCDIRPP := PostProcessing/include
INCDIRS += $(INCDIRPP)
ALLINC += $(wildcard $(INCDIRPP)/*.hh)

SRCDIRTIME := TimeIntegration/src
INCDIRTIME := TimeIntegration/include
INCDIRS += $(INCDIRTIME)
ALLINC += $(wildcard $(INCDIRTIME)/*.hh)

SRCDIRUTIL := Utility/src
INCDIRUTIL := Utility/include
INCDIRS += $(INCDIRUTIL)
ALLINC += $(wildcard $(INCDIRUTIL)/*.hh)

SRCDIRXML := XML3/src
INCDIRXML := XML3/include
INCDIRS += $(INCDIRXML)
ALLINC += $(wildcard $(INCDIRXML)/*.hh)


# Source files per package
# _WP_ means "with path"
SRCFILES_WP_BASE := $(wildcard $(SRCDIRBASE)/*.cpp)
SRCFILESBASE := $(subst $(SRCDIRBASE)/,,$(SRCFILES_WP_BASE))

SRCFILES_WP_COLLISION := $(wildcard $(SRCDIRCOLLISION)/*.cpp)
SRCFILESCOLLISION := $(subst $(SRCDIRCOLLISION)/,,$(SRCFILES_WP_COLLISION))

SRCFILES_WP_COMPO := $(wildcard $(SRCDIRCOMPO)/*.cpp)
SRCFILESCOMPO := $(subst $(SRCDIRCOMPO)/,,$(SRCFILES_WP_COMPO))

SRCFILES_WP_CONTACT := $(wildcard $(SRCDIRCONTACT)/*.cpp)
SRCFILESCONTACT := $(subst $(SRCDIRCONTACT)/,,$(SRCFILES_WP_CONTACT))

SRCFILES_WP_GEOM := $(wildcard $(SRCDIRGEOM)/*.cpp)
SRCFILESGEOM := $(subst $(SRCDIRGEOM)/,,$(SRCFILES_WP_GEOM))

SRCFILES_WP_PP := $(wildcard $(SRCDIRPP)/*.cpp)
SRCFILESPP := $(subst $(SRCDIRPP)/,,$(SRCFILES_WP_PP))

SRCFILES_WP_TIME := $(wildcard $(SRCDIRTIME)/*.cpp)
SRCFILESTIME := $(subst $(SRCDIRTIME)/,,$(SRCFILES_WP_TIME))

SRCFILES_WP_UTIL := $(wildcard $(SRCDIRUTIL)/*.cpp)
SRCFILESUTIL := $(subst $(SRCDIRUTIL)/,,$(SRCFILES_WP_UTIL))

SRCFILES_WP_XML := $(wildcard $(SRCDIRXML)/*.cpp)
SRCFILESXML := $(subst $(SRCDIRXML)/,,$(SRCFILES_WP_XML))


# Substitutions for objects
# Here we use the base file names combined to VPATH
OBJS_BASE := $(SRCFILESBASE:%.cpp=$(OBJDIR)/%.o)
OBJS_COLLISION := $(SRCFILESCOLLISION:%.cpp=$(OBJDIR)/%.o)
OBJS_COMPO := $(SRCFILESCOMPO:%.cpp=$(OBJDIR)/%.o)
OBJS_CONTACT := $(SRCFILESCONTACT:%.cpp=$(OBJDIR)/%.o)
OBJS_GEOM := $(SRCFILESGEOM:%.cpp=$(OBJDIR)/%.o)
OBJS_PP := $(SRCFILESPP:%.cpp=$(OBJDIR)/%.o)
OBJS_TIME := $(SRCFILESTIME:%.cpp=$(OBJDIR)/%.o)
OBJS_XML := $(SRCFILESXML:%.cpp=$(OBJDIR)/%.o)
OBJS_UTIL := $(SRCFILESUTIL:%.cpp=$(OBJDIR)/%.o)
OBJS := $(OBJS_BASE) $(OBJS_COLLISION) $(OBJS_COMPO) $(OBJS_CONTACT) $(OBJS_GEOM) $(OBJS_PP) $(OBJS_TIME) $(OBJS_XML) $(OBJS_UTIL) 
DEPS := $(OBJS:.o=.d)


# Defining VPATH for automatic search of source files in the 
# different source directories 
# This enables to have a single rule/target for all source files in 
# various directories and all objects in a single directory
VPATH := $(SRCDIRBASE):$(SRCDIRCOLLISION):$(SRCDIRCOMPO):$(SRCDIRCONTACT):$(SRCDIRGEOM):$(SRCDIRPP):$(SRCDIRTIME):$(SRCDIRXML):$(SRCDIRUTIL)
CPPFLAGS += $(addprefix -I,$(INCDIRS))
NVCCFLAGS += -O3 -g -Xcompiler -fPIC -x cu -rdc=true --extended-lambda --expt-relaxed-constexpr


#------------------------------------------------------------------------------
# Production rules
.SILENT:

all    : MESS objlib $(LIBSO_WP) dtd tools CPINC

objlib : $(OBJS) $(OBJLIB_MANIFEST)

$(OBJLIB_MANIFEST): $(OBJS) | $(LIBDIR)
	printf 'LIBGRAINS_GPU_OBJECTS :=' > $@
	printf ' %s' $(OBJS) >> $@
	printf '\n' >> $@

OBJ    : GEOM MGEOM XML MXML COLLISION MCOLLISION COMPO MCOMPO CONTACT MCONTACT BASE MBASE TIME MTIME PP MPP UTIL MUTIL
	@echo -e 'Compilation is complete\n'

LINK    : objlib $(LIBSO_WP)

$(LIBSO_WP) : objlib | $(LIBDIR)
	$(NVCC) $(NVCC_ARCHFLAGS) -shared -Xlinker -soname -Xlinker $(LIBSO) $(RPATH_FLAGS) -o $(LIBSO_WP) $(LIBGRAINS_GPU_OBJECTS) $(NVCC_LDLIBS)
	@echo -e '--------------------------------\n''Library libGrainsGPU.so generated !\n''--------------------------------'	

BASE   : $(OBJS_BASE)

COLLISION  : $(OBJS_COLLISION)

COMPO  : $(OBJS_COMPO)

CONTACT  : $(OBJS_CONTACT)

GEOM   : $(OBJS_GEOM)

PP     : $(OBJS_PP)

TIME   : $(OBJS_TIME)
	
XML    : $(OBJS_XML)

UTIL   : $(OBJS_UTIL)

MBASE   : 
	@echo -e '*** Package Base compiled\n'

MCOLLISION  : 
	@echo -e '*** Package CollisionDetection compiled\n'

MCOMPO  : 
	@echo -e '*** Package Component compiled\n'

MCONTACT  : 
	@echo -e '*** Package ContactForce compiled\n'

MGEOM   :
	@echo -e '*** Package Geometry compiled\n'

MPP     : 
	@echo -e '*** Package PostProcessing compiled\n'

MTIME   : 
	@echo -e '*** Package TimeIntegration compiled\n'
	
MXML    : 
	@echo -e '*** Package XML compiled\n'
	
MUTIL   : 
	@echo -e '*** Package Util compiled\n'
	
MESS   :
	@echo -e '************************************************************\n''Grains3DGPU compiled with' $(NVCC) $(CPPFLAGS) $(NVCCFLAGS)'\n************************************************************'

$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(NVCC) $(NVCC_ARCHFLAGS) $(CPPFLAGS) $(NVCCFLAGS) -MMD -MP -MF $(@:.o=.d) -MT $@ -c $< -o $@

$(LIBDIR):
	@mkdir -p $@
	
$(INCDIR):
	@mkdir -p $@
	
CPINC  : | $(INCDIR)
	cp $(ALLINC) $(INCDIR)/	

-include $(DEPS)

clean:
	rm -rf $(OBJDIR)/* \
	rm -f $(LIBSO_WP) \
	rm -f $(OBJLIB_MANIFEST) \
	rm -rf $(SHAREDIR)/* \
	rm -rf $(INCDIR)/*	

dtd:
	@mkdir -p $(SHAREDIR)/Dtd
	@cp -f \
		Dtd/Contact.dtd \
		Dtd/Convex.dtd \
		Dtd/Grains2D_InFluid.dtd \
		Dtd/Grains2D.dtd \
		Dtd/Grains3D_InFluid.dtd \
		Dtd/Grains3D.dtd \
		$(SHAREDIR)/Dtd/

tools:
	@mkdir -p $(SHAREDIR)/Tools/ExecScripts
	@cp -rf \
	  Tools/ExecScripts/*.exec \
		$(SHAREDIR)/Tools/ExecScripts/

-include $(OBJLIB_MANIFEST)
