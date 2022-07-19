update:
	( cd Position/src; make; cd ../.. )
	( cd ShapeFile/src; make; cd ../.. )
	( cd RotationMatrix/src; make; cd ../.. )
	@echo 'PrePost tools built!'

depend:
	( cd Position/src; make depend; cd ../.. )
	( cd ShapeFile/src; make depend; cd ../.. )	
	( cd RotationMatrix/src; make depend; cd ../.. )		

clean:
	( cd Position/src; make clean; cd ../.. )
	( cd ShapeFile/src; make clean; cd ../.. )
	( cd RotationMatrix/src; make clean; cd ../.. )
	@echo 'PrePost tools cleaned'
	
cleandepend:
	( cd Position/src; rm -f Makedepend; cd ../.. )
	( cd ShapeFile/src; rm -f Makedepend; cd ../.. )
	( cd RotationMatrix/src; rm -f Makedepend; cd ../.. )
	@echo 'PrePost tools Makedepend files deleted'			

createdepend:
	( cd Position/src; touch Makedepend; cd ../.. )
	( cd ShapeFile/src; touch Makedepend; cd ../.. )
	( cd RotationMatrix/src; touch Makedepend; cd ../.. )
	@echo 'PrePost tools Makedepend created'	
