default:
	@echo 'No default for this Makefile, run make help for targets'
	
help:
	@echo "------------------------------------------------------------------------------------"
	@echo "Targets:"
	@echo "   clean : cleans basilisk"
	@echo "   download : downloads basilisk from basilisk website using darcs"
	@echo "   update : updates basilisk from basilisk website using darcs"
	@echo "   compile : compiles basilisk with the selected architecture"
#	@echo "   downloadmesaglu : downloads the recommended versions of GLU and MESA from the web"
#	@echo "   mesaglu : compiles the recommended versions of GLU and MESA"
#	@echo "   deletemesaglu : deletes the recommended versions of GLU and MESA"
#	@echo "   cleanall : deletemesaglu clean"
	@echo "------------------------------------------------------------------------------------"

 
clean:
ifeq ($(shell test -d $(BASILISK) && echo -n yes),yes)
	$(RM) $(BASILISK_COMPILED)
	@echo "Compiled basilisk $(BASILISK_COMPILED) deleted"
else
	@echo "No compiled basilisk $(BASILISK_COMPILED) to delete"
endif


download:
ifeq ($(shell test -d $(BASILISK_HEAD) && echo -n yes),yes)
	@echo "Basilisk head directory $(BASILISK_HEAD) already exists, do \"make update\" instead"
else
	darcs get --lazy http://basilisk.fr/basilisk
endif

	
update:
ifeq ($(shell test -d $(BASILISK_HEAD) && echo -n yes),yes)
	cd $(BASILISK_HEAD); darcs pull; cd $(OCTREE_HOME)
else
	@echo "Basilisk head directory $(BASILISK_HEAD) does not exist yet, do \"make download\" first"
endif


compile: clean	
	$(MKDIR) $(BASILISK_COMPILED)
	$(CPDIR) $(BASILISK_HEAD)/src $(BASILISK_COMPILED)/
	cd $(BASILISK); $(CPDIR) config.gcc config; make -k clean; make; cd $(OCTREE_HOME)

		
#downloadmesaglu: deletemesaglu
#	$(MKDIR) $(OCTREE_HOME)/MESAGLU
#	cd $(OCTREE_HOME)/MESAGLU; wget http://basilisk.fr/src/gl/mesa-17.2.4.tar.gz; $(EXPAND) mesa-17.2.4.tar.gz; $(RM) mesa-17.2.4.tar.gz; cd $(OCTREE_HOME)
#	cd $(OCTREE_HOME)/MESAGLU; wget --no-check-certificate https://archive.mesa3d.org/glu/glu-9.0.1.tar.gz; $(EXPAND) glu-9.0.1.tar.gz; $(RM) glu-9.0.1.tar.gz; cd $(OCTREE_HOME)
#
#	
#mesaglu:
#ifeq ($(shell test -d $(OCTREE_HOME)/MESAGLU && echo -n yes),yes)
#	cd $(OCTREE_HOME)/MESAGLU/mesa-17.2.4; ./configure --includedir=$(OCTREE_HOME)/MESAGLU/include --libdir=$(OCTREE_HOME)/MESAGLU/lib$(BASILISK_BITS_DEFAULT)-$(BASILISK_SERCOMPIL_ENV)-$(BASILISK_SERCOMPIL_VERSION) --enable-osmesa --with-gallium-drivers=swrast --disable-driglx-direct --disable-dri --disable-gbm --disable-egl; make; make install; cd $(OCTREE_HOME)
#	cd $(OCTREE_HOME)/MESAGLU/glu-9.0.1; ./configure --includedir=$(OCTREE_HOME)/MESAGLU/include --libdir=$(OCTREE_HOME)/MESAGLU/lib$(BASILISK_BITS_DEFAULT)-$(BASILISK_SERCOMPIL_ENV)-$(BASILISK_SERCOMPIL_VERSION); make; make install; cd $(OCTREE_HOME)
#else
#	@echo "MESAGLU directory does not exist yet, do \"make downloadmesaglu\" first"
#endif
#
#
#deletemesaglu:
#ifeq ($(shell test -d $(OCTREE_HOME)/MESAGLU && echo -n yes),yes)
#	$(RM) MESAGLU
#	@echo "MESAGLU deleted"
#else
#	@echo "No MESAGLU to delete"	
#endif	
#
#	
#cleanall: deletemesaglu clean
	
		
##################################################################
# internal commands                                              #
##################################################################
RM    := rm -fr
MKDIR := mkdir -p
EXPAND := tar xzvf
CPDIR := cp -r

