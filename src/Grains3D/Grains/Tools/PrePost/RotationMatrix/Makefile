# Compilation & linking flags
CXX = $(PACIFIC_CXX)
CC = $(PACIFIC_MPICC)

CXXFLAGS += -O3 -g -MMD -MP
CPPFLAGS = -I$(LIBGRAINS_INCDIR)

LDFLAGS += -L$(LIBGRAINS_LIBDIR) -Wl,-rpath,$(LIBGRAINS_LIBDIR)
LDLIBS += -lGrains

TARGET := $(MAIN_BINDIR)/rotationmatrix
SRC := src/PP_rotationmatrix.cpp
OBJ := $(LIBGRAINS_OBJDIR)/PP_rotationmatrix.o
DEP := $(OBJ:.o=.d)


# -----------------------------------------------------------------------------
# Production rules
.PHONY: all clean

all: $(TARGET)

$(MAIN_BINDIR):
	@mkdir -p $@

$(OBJ): $(SRC) | $(MAIN_BINDIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	$(CXX) -o $@ $(OBJ) $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)

