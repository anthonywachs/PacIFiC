# ----------------
# Standard targets
# ----------------
all: exe0 exe2

exe0: 
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo '~~~ Compiling exe0 linked to libmac0 ~~~'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'	
	make -j4 -f Makefile.mac CCC=$(MAC_FULL_EXT) exe0

exe2: 
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	@echo '~~~ Compiling exe2 linked to libmac2 ~~~'
	@echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
	make -j4 -f Makefile.mac CCC=$(MAC_FULL_EXT) exe2
	
clean:
	$(RM) lib/Linux-$(MAC_FULL_EXT)
	@echo 'Executables for architecture $(MAC_FULL_EXT) deleted'
	
cleanall: 
	$(RM) lib
	@echo 'Executables for all architectures deleted'	


# ----	
# Help
# ----		
help:
	@echo 'Below are the various targets:'
	@echo '   all (default) $(BANG) compile exe0 and exe2'
	@echo '   exe0          $(BANG) compile exe0 linked to libmac0'
	@echo '   exe2          $(BANG) compile exe2 linked to libmac2'	
	@echo '   clean         $(BANG) delete lib/Linux-$(MAC_FULL_EXT)'		
	@echo '   cleanall      $(BANG) delete lib'			

	
##################################################################
# internal commands                                              #
##################################################################
RM := rm -rf
BANG := \#
