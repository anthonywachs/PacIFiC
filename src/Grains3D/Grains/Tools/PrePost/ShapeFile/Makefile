# Compilation & linking flags
CXX = $(PACIFIC_CXX)
CC = $(PACIFIC_MPICC)

CXXFLAGS += -O3 -g -MMD -MP

TARGET := $(MAIN_BINDIR)/shapefile
SRC := src/PP_shapefile.cpp
OBJ := $(LIBGRAINS_OBJDIR)/PP_shapefile.o
DEP := $(OBJ:.o=.d)


# -----------------------------------------------------------------------------
# Production rules
.PHONY: all clean

all: $(TARGET)

$(MAIN_BINDIR):
	@mkdir -p $@

$(OBJ): $(SRC) | $(MAIN_BINDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(TARGET): $(OBJ)
	$(CXX) -o $@ $(OBJ)

clean:
	rm -f $(OBJ) $(DEP) $(TARGET)

