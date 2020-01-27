update:
	cd Main/src ; make ; cd ../..

depend:
	cd Main/src ; make depend ; cd ../..
	cd Grains/src ; make depend ; cd ../..

xerces:
	cd $(XERCES_DIR) ; $(INSTALL_XERCES) ; cd ..

dtd:
	cd Main/dtd ; $(INSTALL_DTD) ; cd ../..

install:
	make xerces ; ./makeARCH create ; make update ; make dtd
	@echo 'ENVIRONNEMENT GRAINS GENERE !'

installdepend: xerces depend update dtd
	@echo 'ENVIRONNEMENT GRAINS GENERE !'

cleanxerces:
	cd $(XERCES_SOURCE) ; make clean ; cd ../../..
	cd $(XERCES_DIR) ; $(RM) lib$(GRAINS_BITS_EXT)-$(GRAINS_SERCOMPIL_ENV)-$(GRAINS_SERCOMPIL_VERSION) ; cd ..
	@echo 'XERCES supprime'

clean:
	cd Main/src ; make clean ; cd ../..
	@echo 'GRAINS supprime'

cleandtd:
	cd Main/dtd ; $(RM) Grains*.dtd ; cd ../..

cleanall:
	make cleanxerces ; make clean ; ./makeARCH "rm -r"; make cleandtd


##################################################################
# internal commands                                              #
##################################################################
RM := rm -rf
INSTALL_XERCES := ./install.sh
XERCES_DIR := XERCES-2.8.0
XERCES_SOURCE := XERCES-2.8.0/src/xercesc
INSTALL_DTD := ./installdtd.sh
# DO NOT DELETE
