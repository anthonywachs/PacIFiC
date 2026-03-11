CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)
CPPFLAGS += $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS)
LDLIBS += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS)

OBJDIR := $(FLUID_OBJDIR)
BINDIR := $(FLUID_BINDIR)
TARGET := fluid
TARGET_WP := $(BINDIR)/$(TARGET)

ALLINC =

# External dependencies
MAC_INCDIRS := \
  $(MAC_INCDIR)
CPPFLAGS += $(addprefix -I,$(MAC_INCDIRS))
LDLIBS += -L$(MAC_LIBDIR) -lMAC
LDFLAGS += -L$(MAC_LIBDIR) -Wl,-rpath,$(MAC_LIBDIR)

GRAINS_INCDIRS := \
  $(GRAINS_INCDIR)
CPPFLAGS += $(addprefix -I,$(GRAINS_INCDIRS))
LDLIBS += -L$(GRAINS_LIBDIR) -lGrains
LDFLAGS += -L$(GRAINS_LIBDIR) -Wl,-rpath,$(GRAINS_LIBDIR)

# Source and include directories
SRCDIRFSI := FluidSolidInterface/src
INCDIRFSI := FluidSolidInterface/include
INCDIRS += $(INCDIRFSI)
ALLINC += $(wildcard $(INCDIRFSI)/*.hh)

SRCDIRMISC := Misc/src
INCDIRMISC := Misc/include
INCDIRS += $(INCDIRMISC)
ALLINC += $(wildcard $(INCDIRMISC)/*.hh)

SRCDIRPP := PostProcessing/src
INCDIRPP := PostProcessing/include
INCDIRS += $(INCDIRPP)
ALLINC += $(wildcard $(INCDIRPP)/*.hh)

SRCDIRDS := Solvers/DirectionSplitting/src
INCDIRDS := Solvers/DirectionSplitting/include
INCDIRS += $(INCDIRDS)
ALLINC += $(wildcard $(INCDIRDS)/*.hh)

SRCDIRDLMFD := Solvers/DLMFD/src
INCDIRDLMFD := Solvers/DLMFD/include
INCDIRS += $(INCDIRDLMFD)
ALLINC += $(wildcard $(INCDIRDLMFD)/*.hh)

SRCDIRVPH := Solvers/Viscoplastic_HeatTransfer/src
INCDIRVPH := Solvers/Viscoplastic_HeatTransfer/include
INCDIRS += $(INCDIRVPH)
ALLINC += $(wildcard $(INCDIRVPH)/*.hh)

# Source files per package
SRCFILES_WP_FSI := $(wildcard $(SRCDIRFSI)/*.cc)
SRCFILESFSI := $(subst $(SRCDIRFSI)/,,$(SRCFILES_WP_FSI))

SRCFILES_WP_MISC := $(wildcard $(SRCDIRMISC)/*.cc)
SRCFILESMISC := $(subst $(SRCDIRMISC)/,,$(SRCFILES_WP_MISC))

SRCFILES_WP_PP := $(wildcard $(SRCDIRPP)/*.cc)
SRCFILESPP := $(subst $(SRCDIRPP)/,,$(SRCFILES_WP_PP))

SRCFILES_WP_DS := $(wildcard $(SRCDIRDS)/*.cc)
SRCFILESDS := $(subst $(SRCDIRDS)/,,$(SRCFILES_WP_DS))

SRCFILES_WP_DLMFD := $(wildcard $(SRCDIRDLMFD)/*.cc)
SRCFILESDLMFD := $(subst $(SRCDIRDLMFD)/,,$(SRCFILES_WP_DLMFD))

SRCFILES_WP_VPH := $(wildcard $(SRCDIRVPH)/*.cc)
SRCFILESVPH := $(subst $(SRCDIRVPH)/,,$(SRCFILES_WP_VPH))

# Objects
OBJS_FSI := $(SRCFILESFSI:%.cc=$(OBJDIR)/%.o)
OBJS_MISC := $(SRCFILESMISC:%.cc=$(OBJDIR)/%.o)
OBJS_PP := $(SRCFILESPP:%.cc=$(OBJDIR)/%.o)
OBJS_DS := $(SRCFILESDS:%.cc=$(OBJDIR)/%.o)
OBJS_DLMFD := $(SRCFILESDLMFD:%.cc=$(OBJDIR)/%.o)
OBJS_VPH := $(SRCFILESVPH:%.cc=$(OBJDIR)/%.o)
OBJS := $(OBJS_FSI) $(OBJS_MISC) $(OBJS_PP) $(OBJS_DS) $(OBJS_DLMFD) $(OBJS_VPH)
DEPS := $(OBJS:.o=.d)

VPATH := $(SRCDIRFSI):$(SRCDIRMISC):$(SRCDIRPP):$(SRCDIRDS):$(SRCDIRDLMFD):$(SRCDIRVPH)

CPPFLAGS += $(addprefix -I,$(INCDIRS))
CXXFLAGS += -O3 -g -MMD -MP

.SILENT:

all: MESS $(TARGET_WP)

OBJ: FSI MFSI MISC MMISC PP MPP DS MDS DLMFD MDLMFD VPH MVPH
	@echo -e 'Compilation is complete\n'

$(TARGET_WP): $(OBJS) | $(BINDIR)
	$(CXX) -o $(TARGET_WP) $(OBJS) $(LDFLAGS) $(LDLIBS)
	@echo -e '-------------------------------\n''Executable fluid generated !\n''-------------------------------'

FSI: $(OBJS_FSI)

MISC: $(OBJS_MISC)

PP: $(OBJS_PP)

DS: $(OBJS_DS)

DLMFD: $(OBJS_DLMFD)

VPH: $(OBJS_VPH)

MFSI:
	@echo -e '*** Package FluidSolidInterface compiled\n'

MMISC:
	@echo -e '*** Package Misc compiled\n'

MPP:
	@echo -e '*** Package PostProcessing compiled\n'

MDS:
	@echo -e '*** Package Solvers/DirectionSplitting compiled\n'

MDLMFD:
	@echo -e '*** Package Solvers/DLMFD compiled\n'

MVPH:
	@echo -e '*** Package Solvers/Viscoplastic_HeatTransfer compiled\n'

MESS:
	@echo -e '************************************************************\n''FLUID compiled with' $(CXX) $(CPPFLAGS) $(CXXFLAGS)'\n************************************************************'

$(OBJDIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	@echo 'Compiling '$<
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BINDIR):
	@mkdir -p $@

clean:
	rm -rf $(OBJDIR)
	rm -f $(TARGET_WP)

-include $(DEPS)
