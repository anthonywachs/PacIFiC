CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
CPPFLAGS += $(PACIFIC_ZLIB_CFLAGS) $(PACIFIC_PETSC_CFLAGS)
LDLIBS   += $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_PETSC_LFLAGS)

OBJDIR := $(MAC_OBJDIR)
LIBDIR := $(MAC_LIBDIR)
INCDIR := $(MAC_INCDIR)

LIBSO := libMAC.so
LIBSO_WP := $(LIBDIR)/$(LIBSO)

ALLINC = 

# Source and include directories
SRCDIRARRV := ArrayVector/src
INCDIRARRV := ArrayVector/include
INCDIRS += $(INCDIRARRV)
ALLINC += $(wildcard $(INCDIRARRV)/*.hh)
ALLINC += $(wildcard $(INCDIRARRV)/*.icc)

SRCDIRBASE := Base/src
INCDIRBASE := Base/include
INCDIRS += $(INCDIRBASE)
ALLINC += $(wildcard $(INCDIRBASE)/*.hh)
ALLINC += $(wildcard $(INCDIRBASE)/*.icc)

SRCDIRMPI := ExternalAPI/MPI/src
INCDIRMPI := ExternalAPI/MPI/include
INCDIRS += $(INCDIRMPI)
ALLINC += $(wildcard $(INCDIRMPI)/*.hh)
ALLINC += $(wildcard $(INCDIRMPI)/*.icc)

SRCDIRPETSC := ExternalAPI/PETSc_3.24.3/src
INCDIRPETSC := ExternalAPI/PETSc_3.24.3/include
INCDIRS += $(INCDIRPETSC)
ALLINC += $(wildcard $(INCDIRPETSC)/*.hh)
ALLINC += $(wildcard $(INCDIRPETSC)/*.icc)

SRCDIRFV := FiniteVolume/src
INCDIRFV := FiniteVolume/include
INCDIRS += $(INCDIRFV)
ALLINC += $(wildcard $(INCDIRFV)/*.hh)
ALLINC += $(wildcard $(INCDIRFV)/*.icc)

SRCDIRLA := LinearAlgebra/src
INCDIRLA := LinearAlgebra/include
INCDIRS += $(INCDIRLA)
ALLINC += $(wildcard $(INCDIRLA)/*.hh)
ALLINC += $(wildcard $(INCDIRLA)/*.icc)

# Source files per package
# _WP_ means "with path"
SRCFILES_WP_ARRV := $(wildcard $(SRCDIRARRV)/*.cc)
SRCFILESARRV := $(subst $(SRCDIRARRV)/,,$(SRCFILES_WP_ARRV))

SRCFILES_WP_BASE := $(wildcard $(SRCDIRBASE)/*.cc)
SRCFILESBASE := $(subst $(SRCDIRBASE)/,,$(SRCFILES_WP_BASE))

SRCFILES_WP_MPI := $(wildcard $(SRCDIRMPI)/*.cc)
SRCFILESMPI := $(subst $(SRCDIRMPI)/,,$(SRCFILES_WP_MPI))

SRCFILES_WP_PETSC := $(wildcard $(SRCDIRPETSC)/*.cc)
SRCFILESPETSC := $(subst $(SRCDIRPETSC)/,,$(SRCFILES_WP_PETSC))

SRCFILES_WP_FV := $(wildcard $(SRCDIRFV)/*.cc)
SRCFILESFV := $(subst $(SRCDIRFV)/,,$(SRCFILES_WP_FV))

SRCFILES_WP_LA := $(wildcard $(SRCDIRLA)/*.cc)
SRCFILESLA := $(subst $(SRCDIRLA)/,,$(SRCFILES_WP_LA))

# Substitutions for objects
# Here we use the base file names combined to VPATH
OBJS_ARRV := $(SRCFILESARRV:%.cc=$(OBJDIR)/%.o)
OBJS_BASE := $(SRCFILESBASE:%.cc=$(OBJDIR)/%.o)
OBJS_MPI := $(SRCFILESMPI:%.cc=$(OBJDIR)/%.o)
OBJS_PETSC := $(SRCFILESPETSC:%.cc=$(OBJDIR)/%.o)
OBJS_FV := $(SRCFILESFV:%.cc=$(OBJDIR)/%.o)
OBJS_LA := $(SRCFILESLA:%.cc=$(OBJDIR)/%.o)
OBJS := $(OBJS_ARRV) $(OBJS_BASE) $(OBJS_MPI) $(OBJS_PETSC) $(OBJS_FV) $(OBJS_LA)
DEPS := $(OBJS:.o=.d)

VPATH := $(SRCDIRARRV):$(SRCDIRBASE):$(SRCDIRMPI):$(SRCDIRPETSC):$(SRCDIRFV):$(SRCDIRLA)

CPPFLAGS += $(addprefix -I,$(INCDIRS))
CPPFLAGS += -DMPIRUN=\"$(PACIFIC_MPIEXEC)\"
CPPFLAGS += -DLEVEL=0
CXXFLAGS += -O3 -g -MMD -MP -fPIC

#------------------------------------------------------------------------------
# Production rules
.SILENT:

all    : MESS $(LIBSO_WP) CPINC

OBJ    : ARRV MARRV BASE MBASE MPI MMPI PETSC MPETSC FV MFV LA MLA
	@echo -e 'Compilation is complete\n'

LINK : $(LIBSO_WP)

$(LIBSO_WP) : $(OBJS) | $(LIBDIR)
	$(CXX) -shared -Wl,-soname,$(LIBSO) -o $(LIBSO_WP) $(OBJS) $(LDFLAGS) $(LDLIBS)
	@echo -e '--------------------------------\n''Library libMAC.so generated !\n''--------------------------------'	

ARRV : $(OBJS_ARRV)

BASE  : $(OBJS_BASE)

MPI  : $(OBJS_MPI)

PETSC   : $(OBJS_PETSC)

FV     : $(OBJS_FV)

LA   : $(OBJS_LA)
	
MARRV     : 
	@echo -e '*** Package ArrayVector compiled\n'

MBASE   : 
	@echo -e '*** Package Base compiled\n'

MMPI    : 
	@echo -e '*** Package ExternalAPI/MPI compiled\n'
	
MPETSC    : 
	@echo -e '*** Package ExternalAPI/PETSc compiled\n'
	
MFV   : 
	@echo -e '*** Package FiniteVolume compiled\n'
	
MLA    : 
	@echo -e '*** Package LinearAlgebra compiled\n'	
		
MESS   :
	@echo -e '************************************************************\n''MAC compiled with' $(CXX) $(CPPFLAGS) $(CXXFLAGS)'\n************************************************************'	

$(OBJDIR)/%.o : %.cc
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(LIBDIR):
	@mkdir -p $@
	
$(INCDIR):
	@mkdir -p $@
	
CPINC  : $(INCDIR)
	cp $(ALLINC) $(INCDIR)/	

-include $(DEPS)

clean:
	rm -rf $(OBJDIR)/*.o \
	rm -f $(LIBSO_WP) \
	rm -rf $(INCDIR)/*.hh \
	rm -rf $(INCDIR)/*.icc
 
