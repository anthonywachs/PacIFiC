# ----------------
# Standard targets
# ----------------
update:
	@echo 'Updating Grains3D Doxygen doc'	
	cd ../Grains/ ; make CPINC; cd ../DoxyGRAINS/
	doxygen Doxyfile

cleanhtml:
	@echo 'Deleting html Doxygen files'	
	$(RM) html/* 

cleanlatex:
	@echo 'Deleting latex Doxygen files'	
	$(RM) latex/*

cleanall: cleanhtml cleanlatex

install: cleanall update
	@echo 'Full Grains platform built!'
	

# ----	
# Help
# ----		
help:
	@echo 'Below are the various targets:'
	@echo '   STANDARD TARGETS:'
	@echo '      update (default) $(BANG) update Grains3D Doxygen doc'
	@echo '      cleanhtml        $(BANG) delete html Doxygen files'
	@echo '      cleanlatex       $(BANG) delete latex Doxygen files'	
	@echo '      cleanall         $(BANG) delete both html and latex Doxygen filelatex Doxygen file'		
	@echo '      install          $(BANG) clean and generate Doxygen files'			

	
##################################################################
# internal commands                                              #
##################################################################
TOUCH := touch
RM := rm -rf
INSTALL_XERCES := ./install.sh
XERCES_DIR := XERCES-2.8.0
XERCES_SOURCE := XERCES-2.8.0/src/xercesc
INSTALL_DTD := ./installdtd.sh
BANG := \#
