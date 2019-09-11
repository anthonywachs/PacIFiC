all:../src/MAC_Yacc.cc ../src/MAC_Lex.cc 

../src/MAC_Lex.cc:Gram.l
	flex -PMAC -d Gram.l
	mv lex.MAC.c ../src/MAC_Lex.cc

../src/MAC_Yacc.cc:Gram.y
	bison -p MAC -d Gram.y -o y.tab.c
	mv y.tab.c ../src/MAC_Yacc.cc
	mv y.tab.h ../include/MAC_Yacc.hh

clean:
	/bin/rm -f lex.MAC.c y.tab.*  y.output ../src/MAC_Lex.cc ../src/MAC_Yacc.cc ../include/MAC_Yacc.hh
