# ----------------
# Targets
# ----------------
MAKEFLAGS += --no-print-directory

.SILENT:

update:
	make -f GrainsRecipes MESS
	make -f GrainsRecipes CPINC
	make -f GrainsRecipes START
	make -j4 --output-sync=line -f GrainsRecipes OBJ
	make -f GrainsRecipes LINK
	make -f GrainsRecipes BASILISK	
	
clean:
	make -f GrainsRecipes clean
	
cpinc:
	make -f GrainsRecipes CPINC
	
depend:
	make -f GrainsRecipes depend
