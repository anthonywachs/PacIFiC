
#standard targets:
# all   : builds the target (TARGET variable) in BINDIR
# obj   : builds the objects only
# clean : clean the dependencies, the objects files and the target
# cleandepend : clean the dependencies
# cleanobj    : clean the objects files
# cleantarget : clean the target file

.PHONY : default_target all dirs
.PHONY : clean cleandepend cleanobj cleantarget
.SUFFIXES :
.SUFFIXES : .c .cc .d .f .F
default_target: all

###############################################################################
#default link mode (cc,c,f,...)
###############################################################################
LINK_MODE=cc

###############################################################################
#declare here the final target and the start sources (which may be in
#different dirs).
###############################################################################
TARGET = exe

MKDIR    = mkdir -p
RMDIR    = rm -fr
DEPFILTER = perl -e 'while(<>) {/^\s*\/usr\// && next;s{^($(*F)).o\s*:}{$(@D)/\1.o $(@D)/\1.d :}; print;}print "\n";'

###############################################################################
#object files definition
NOTDIRSRC := $(notdir $(SRC))
VPATH:=$(subst ' ',':',$(sort $(dir $(SRC))))
OBJ := $(patsubst %.c,   %.o, $(filter %.c,   $(NOTDIRSRC)))\
       $(patsubst %.cc,  %.o, $(filter %.cc,  $(NOTDIRSRC)))\
       $(patsubst %.cpp, %.o, $(filter %.cpp, $(NOTDIRSRC)))\
       $(patsubst %.F,  %.o, $(filter %.F,  $(NOTDIRSRC)))\
       $(patsubst %.f,  %.o, $(filter %.f,  $(NOTDIRSRC)))

OBJ_X = $(addprefix $(BINDIR), $(OBJ))

#target file definition
TARGET_X += $(addprefix $(BINDIR)/, $(TARGET))

#(sub)directories creation in BINDIR
ifdef OBJ_X
DIRS := $(shell $(MKDIR) $(sort $(dir $(OBJ_X))))
endif

#dependencies inclusion
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),cleandepend)
ifneq ($(MAKECMDGOALS),cleanobj)
ifneq ($(MAKECMDGOALS),cleantarget)
ifdef OBJ_X
include $(OBJ_X:.o=.d)
endif
endif
endif
endif
endif

#standard targets
all:$(TARGET_X)
obj: $(OBJ_X)
cleandepend: ; $(RM) $(OBJ_X:.o=.d)
cleanobj:    ; $(RM) $(OBJ_X)
cleantarget: ; $(RM) $(TARGET_X)
clean       : cleandepend cleanobj cleantarget

#the target cmd depends on the suffix (.a, .so or other (executable))
$(TARGET_X):$(OBJ_X)
ifeq ($(suffix $(TARGET_X)),.a)
	$(AR) $(ARFLAGS) $@ $?
else
	$(RM) $(TARGET_X)
ifeq ($(suffix $(TARGET_X)),$(DYNAMIC_LIB_EXT))
	 cd $(BINDIR) ; $(LD.so) $(LDFLAGS.so) $(OBJ) $(PRECOMP_OBJ) $(LDLIBSSO) $(LDLIBS)
else 
	$(LINK.$(LINK_MODE)) -o $@ $(LDFLAGS) $(OBJ_X) $(PRECOMP_OBJ) $(LDLIBS)
endif
endif
###############################################################################
# Production rules ... Don't look further !
###############################################################################

$(BINDIR)%.d:%.F   ; touch $@
$(BINDIR)%.d:%.f   ; touch $@
$(BINDIR)%.d:%.c 
	@echo "Building dependencies for $<"
	@$(MKDEP.c) $< | $(DEPFILTER) > $@
$(BINDIR)%.d:%.cc 
	@echo "Building dependencies for $<"
	@$(MKDEP.c) $< | $(DEPFILTER) > $@
$(BINDIR)%.d:%.cpp 
	@echo "Building dependencies for $<"
	@$(MKDEP.c) $< | $(DEPFILTER) > $@

$(BINDIR)%.o:%.f  ;	$(COMPILE.f)   $< -o $@
$(BINDIR)%.o:%.F  ;	$(COMPILE.f)   $< -o $@
$(BINDIR)%.o:%.c  ;	$(COMPILE.c)   $< -o $@
$(BINDIR)%.o:%.cc ;	$(COMPILE.cc)  $< -o $@
$(BINDIR)%.o:%.cpp ;	$(COMPILE.cc)  $< -o $@
