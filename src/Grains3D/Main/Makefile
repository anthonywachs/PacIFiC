# Compilation & linking flags
CXX = $(PACIFIC_MPICXX)
CC = $(PACIFIC_MPICC)

CPPFLAGS += -I$(LIBGRAINS_INCDIR)
CPPFLAGS += $(PACIFIC_MPI_SKIP_CPPFLAGS) $(PACIFIC_XERCESC_CFLAGS) $(PACIFIC_ZLIB_CFLAGS)

CXXFLAGS += -O3 -g -MMD -MP

LDFLAGS += -L$(LIBGRAINS_LIBDIR) -Wl,-rpath,$(LIBGRAINS_LIBDIR)
LDLIBS += -lGrains3D
LDLIBS += $(PACIFIC_XERCESC_LFLAGS) $(PACIFIC_ZLIB_LFLAGS)


# Target directories
LIBGRAINS_SO := $(LIBGRAINS_LIBDIR)/libGrains3D.so

TARGET := $(MAIN_BINDIR)/grains3d
SRC := src/main.cpp
OBJ := $(LIBGRAINS_OBJDIR)/main.o
DEP := $(OBJ:.o=.d)


# -----------------------------------------------------------------------------
# Production rules
.PHONY: all clean

all: $(TARGET)

$(MAIN_BINDIR):
	@mkdir -p $@

$(OBJ): $(SRC) | $(MAIN_BINDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJ) $(LIBGRAINS_SO)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS) $(LDLIBS)
	@echo -e '-----------------------\n''Software grains3d built !\n''-----------------------'

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)

-include $(DEP)
