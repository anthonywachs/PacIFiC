CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_ZLIB_CFLAGS) $(PACIFIC_PETSC_CFLAGS) $(PACIFIC_HDF5_CFLAGS)
LDLIBS   += $(PACIFIC_ZLIB_LFLAGS) $(PACIFIC_PETSC_LFLAGS) $(PACIFIC_HDF5_LFLAGS)

OBJDIR := $(MAC_OBJDIR)
LIBDIR := $(MAC_LIBDIR)
INCDIR := $(MAC_INCDIR)

LIBSO := libMAC.so
LIBSO_WP = $(LIBDIR)/$(LIBSO)

LIBSO0 := libMAC0.so
LIBSO0_WP := $(LIBDIR)/$(LIBSO0)

LIBSO2 := libMAC2.so
LIBSO2_WP := $(LIBDIR)/$(LIBSO2)


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

MAC_PETSC_API_VERSION ?= 3.25.3
SRCDIRPETSC := ExternalAPI/PETSc_$(MAC_PETSC_API_VERSION)/src
INCDIRPETSC := ExternalAPI/PETSc_$(MAC_PETSC_API_VERSION)/include
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
OBJS_ARRV_0 := $(SRCFILESARRV:%.cc=$(OBJDIR)/%.0.o)
OBJS_BASE_0 := $(SRCFILESBASE:%.cc=$(OBJDIR)/%.0.o)
OBJS_MPI_0 := $(SRCFILESMPI:%.cc=$(OBJDIR)/%.0.o)
OBJS_PETSC_0 := $(SRCFILESPETSC:%.cc=$(OBJDIR)/%.0.o)
OBJS_FV_0 := $(SRCFILESFV:%.cc=$(OBJDIR)/%.0.o)
OBJS_LA_0 := $(SRCFILESLA:%.cc=$(OBJDIR)/%.0.o)
OBJS_0 := $(OBJS_ARRV_0) $(OBJS_BASE_0) $(OBJS_MPI_0) $(OBJS_PETSC_0) $(OBJS_FV_0) $(OBJS_LA_0)

OBJS_ARRV_2 := $(SRCFILESARRV:%.cc=$(OBJDIR)/%.2.o)
OBJS_BASE_2 := $(SRCFILESBASE:%.cc=$(OBJDIR)/%.2.o)
OBJS_MPI_2 := $(SRCFILESMPI:%.cc=$(OBJDIR)/%.2.o)
OBJS_PETSC_2 := $(SRCFILESPETSC:%.cc=$(OBJDIR)/%.2.o)
OBJS_FV_2 := $(SRCFILESFV:%.cc=$(OBJDIR)/%.2.o)
OBJS_LA_2 := $(SRCFILESLA:%.cc=$(OBJDIR)/%.2.o)
OBJS_2 := $(OBJS_ARRV_2) $(OBJS_BASE_2) $(OBJS_MPI_2) $(OBJS_PETSC_2) $(OBJS_FV_2) $(OBJS_LA_2)

DEPS := $(OBJS_0:.o=.d) $(OBJS_2:.o=.d)

VPATH := $(SRCDIRARRV):$(SRCDIRBASE):$(SRCDIRMPI):$(SRCDIRPETSC):$(SRCDIRFV):$(SRCDIRLA)

CPPFLAGS += $(addprefix -I,$(INCDIRS))
CPPFLAGS += -DMPIRUN=\"$(PACIFIC_MPIEXEC)\"
CXXFLAGS += -O3 -g -MMD -MP -fPIC

CPPFLAGS_LEVEL0 = -DLEVEL=0
CPPFLAGS_LEVEL2 = -DLEVEL=2

#------------------------------------------------------------------------------
# Production rules
.SILENT:
.PHONY: FORCE

all    : MESS $(LIBSO_WP) CPINC

OBJ    : ARRV MARRV BASE MBASE MPI MMPI PETSC MPETSC FV MFV LA MLA
	@echo -e 'Compilation is complete\n'

LINK : $(LIBSO0_WP) $(LIBSO2_WP)

$(LIBSO0_WP) : $(OBJS_0) | $(LIBDIR)
	$(CXX) -shared -Wl,-soname,$(LIBSO0) -o $(LIBSO0_WP) $(OBJS_0) $(LDFLAGS) $(LDLIBS)
	@echo -e '--------------------------------\n''Library libMAC0.so generated !\n''--------------------------------'	

$(LIBSO2_WP) : $(OBJS_2) | $(LIBDIR)
	$(CXX) -shared -Wl,-soname,$(LIBSO2) -o $(LIBSO2_WP) $(OBJS_2) $(LDFLAGS) $(LDLIBS)
	@echo -e '--------------------------------\n''Library libMAC2.so generated !\n''--------------------------------'	


ifeq ($(MAC_LEVEL),0)
$(LIBSO_WP) : $(LIBSO0_WP) FORCE | $(LIBDIR)
	@cd $(LIBDIR) && ln -sf $(LIBSO0) $(LIBSO)
else ifeq ($(MAC_LEVEL),2)
$(LIBSO_WP) : $(LIBSO2_WP) FORCE | $(LIBDIR)
	@cd $(LIBDIR) && ln -sf $(LIBSO2) $(LIBSO)
else
$(LIBSO_WP) : $(LIBSO0_WP) FORCE | $(LIBDIR)
	@cd $(LIBDIR) && ln -sf $(LIBSO0) $(LIBSO)
endif

ARRV : $(OBJS_ARRV_0) $(OBJS_ARRV_2)

BASE  : $(OBJS_BASE_0) $(OBJS_BASE_2)

MPI  : $(OBJS_MPI_0) $(OBJS_MPI_2)

PETSC   : $(OBJS_PETSC_0) $(OBJS_PETSC_2)

FV     : $(OBJS_FV_0) $(OBJS_FV_2)

LA   : $(OBJS_LA_0) $(OBJS_LA_2)
	
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
	@echo -e '************************************************************\n''MAC compiled with' $(CXX) $(CPPFLAGS) $(CXXFLAGS)'\n''LEVEL=0 ->' $(LIBSO0) '\n''LEVEL=2 ->' $(LIBSO2) '\n************************************************************'	

$(OBJDIR)/%.0.o : %.cc
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(CXX) $(CPPFLAGS) $(CPPFLAGS_LEVEL0) $(CXXFLAGS) -c $< -o $@

$(OBJDIR)/%.2.o : %.cc
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(CXX) $(CPPFLAGS) $(CPPFLAGS_LEVEL2) $(CXXFLAGS) -c $< -o $@

$(LIBDIR):
	@mkdir -p $@
	
$(INCDIR):
	@mkdir -p $@
	
CPINC  : $(INCDIR)
	cp $(ALLINC) $(INCDIR)/	

-include $(DEPS)

clean:
	rm -rf $(OBJDIR)/*.0.o $(OBJDIR)/*.0.d
	rm -rf $(OBJDIR)/*.2.o $(OBJDIR)/*.2.d
	rm -f $(LIBSO0_WP) $(LIBSO2_WP)
	rm -rf $(INCDIR)/*.hh
	rm -rf $(INCDIR)/*.icc
 
