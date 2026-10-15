empty :=
space := $(empty) $(empty)
comma := ,

# Convert GCC/Clang-style linker forwarding flags to NVCC style.
# Example:
#   -Wl,-rpath,/path -> -Xlinker -rpath -Xlinker /path
to_nvcc_ldflags = \
  $(filter-out -Wl$(comma)%,$(1)) \
  $(foreach grp,$(patsubst -Wl$(comma)%,%,$(filter -Wl$(comma)%,$(1))), \
    $(foreach opt,$(subst $(comma),$(space),$(grp)),-Xlinker $(opt)))
