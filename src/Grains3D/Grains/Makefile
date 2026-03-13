# Compilation & linking flags
CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
CPPFLAGS += $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS)
LDLIBS   += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS)

# Target directories
OBJDIR := $(LIBGRAINS_OBJDIR)
LIBDIR := $(LIBGRAINS_LIBDIR)
SHAREDIR := $(LIBGRAINS_SHAREDIR)
INCDIR := $(LIBGRAINS_INCDIR)

LIBSO := libGrains3D.so
LIBSO_WP := $(LIBDIR)/$(LIBSO)

ifeq ($(strip $(SHAREDIR)),)
SHAREDIR := $(LIBDIR)/../share
endif

ALLINC = 

# Source and include directories
SRCDIRBASE := Base/src
INCDIRBASE := Base/include
INCDIRS += $(INCDIRBASE)
ALLINC += $(wildcard $(INCDIRBASE)/*.hh)

SRCDIRCOMPO := Component/src
INCDIRCOMPO := Component/include
INCDIRS += $(INCDIRCOMPO)
ALLINC += $(wildcard $(INCDIRCOMPO)/*.hh)

SRCDIRFLUID := Fluid/src
INCDIRFLUID := Fluid/include
INCDIRS += $(INCDIRFLUID)
ALLINC += $(wildcard $(INCDIRFLUID)/*.hh)

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

SRCDIRBAS := InterfaceCBasilisk/src
INCDIRBAS := InterfaceCBasilisk/include
INCDIRS += $(INCDIRBAS)
ALLINC += $(wildcard $(INCDIRBAS)/*.h)


# Source files per package
# _WP_ means "with path"
SRCFILES_WP_BASE := $(wildcard $(SRCDIRBASE)/*.cpp)
SRCFILESBASE := $(subst $(SRCDIRBASE)/,,$(SRCFILES_WP_BASE))

SRCFILES_WP_COMPO := $(wildcard $(SRCDIRCOMPO)/*.cpp)
SRCFILESCOMPO := $(subst $(SRCDIRCOMPO)/,,$(SRCFILES_WP_COMPO))

SRCFILES_WP_FLUID := $(wildcard $(SRCDIRFLUID)/*.cpp)
SRCFILESFLUID := $(subst $(SRCDIRFLUID)/,,$(SRCFILES_WP_FLUID))

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

SRCFILES_WP_BAS := $(wildcard $(SRCDIRBAS)/*.cpp)
SRCFILESBAS := $(subst $(SRCDIRBAS)/,,$(SRCFILES_WP_BAS))


# Substitutions for objects
# Here we use the base file names combined to VPATH
OBJS_BASE := $(SRCFILESBASE:%.cpp=$(OBJDIR)/%.o)
OBJS_COMPO := $(SRCFILESCOMPO:%.cpp=$(OBJDIR)/%.o)
OBJS_FLUID := $(SRCFILESFLUID:%.cpp=$(OBJDIR)/%.o)
OBJS_GEOM := $(SRCFILESGEOM:%.cpp=$(OBJDIR)/%.o)
OBJS_PP := $(SRCFILESPP:%.cpp=$(OBJDIR)/%.o)
OBJS_TIME := $(SRCFILESTIME:%.cpp=$(OBJDIR)/%.o)
OBJS_XML := $(SRCFILESXML:%.cpp=$(OBJDIR)/%.o)
OBJS_UTIL := $(SRCFILESUTIL:%.cpp=$(OBJDIR)/%.o)
OBJS_BAS := $(SRCFILESBAS:%.cpp=$(OBJDIR)/%.o)
OBJS := $(OBJS_BASE) $(OBJS_COMPO) $(OBJS_FLUID) $(OBJS_GEOM) $(OBJS_PP) $(OBJS_TIME) $(OBJS_XML) $(OBJS_UTIL) $(OBJS_BAS)
DEPS := $(OBJS:.o=.d)


# Defining VPATH for automatic search of source files in the 
# different source directories 
# This enables to have a single rule/target for all source files in 
# various directories and all objects in a single directory
VPATH := $(SRCDIRBASE):$(SRCDIRCOMPO):$(SRCDIRFLUID):$(SRCDIRGEOM):$(SRCDIRPP):$(SRCDIRTIME):$(SRCDIRXML):$(SRCDIRUTIL):$(SRCDIRBAS)


CPPFLAGS += $(addprefix -I,$(INCDIRS))

CXXFLAGS += -O3 -g -MMD -MP -fPIC


#------------------------------------------------------------------------------
# Production rules
.SILENT:

all    : MESS $(LIBSO_WP) dtd tools CPINC

OBJ    : GEOM MGEOM XML MXML COMPO MCOMPO BASE MBASE FLUID MFLUID BAS MBAS TIME MTIME PP MPP UTIL MUTIL
	@echo -e 'Compilation is complete\n'

LINK    : $(LIBSO_WP)

$(LIBSO_WP) : $(OBJS) | $(LIBDIR)
	$(CXX) -shared -Wl,-soname,$(LIBSO) -o $(LIBSO_WP) $(OBJS) $(LDFLAGS) $(LDLIBS)
	@echo -e '--------------------------------\n''Library libGrains.so generated !\n''--------------------------------'	

BASE   : $(OBJS_BASE)

COMPO  : $(OBJS_COMPO)

FLUID  : $(OBJS_FLUID)

GEOM   : $(OBJS_GEOM)

PP     : $(OBJS_PP)

TIME   : $(OBJS_TIME)
	
XML    : $(OBJS_XML)

UTIL   : $(OBJS_UTIL)

BAS    : $(OBJS_BAS)

MBASE   : 
	@echo -e '*** Package Base compiled\n'

MCOMPO  : 
	@echo -e '*** Package Component compiled\n'

MFLUID  : 
	@echo -e '*** Package Fluid compiled\n'

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
	
MBAS    : 
	@echo -e '*** Package Basilisk compiled\n'	
		
MESS   :
	@echo -e '************************************************************\n''Grains3D compiled with' $(CXX) $(CPPFLAGS) $(CXXFLAGS)'\n************************************************************'	

$(OBJDIR)/%.o : %.cpp
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

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
