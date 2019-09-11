# Compiler package

package Compiler ;
use strict ;
use diagnostics ;
use Defs ;
use File::Spec ;

BEGIN {
  require Exporter;
  
  use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
  $VERSION = 1.0;
  @ISA = qw(Exporter);
  @EXPORT = qw(new run native_compiler add_macrodef add_library add_library_path add_include print do_profiling );
}

use vars @EXPORT_OK ;

sub native_compiler {
  my $sys = shift @_ ;
  my $result = "gcc" ;
  if( $sys =~ "SunOS" ) { $result = "CC" } ;
  if( $sys =~ "OSF1" ) { $result = "cxx" } ;
  $result ;
}

sub new {
  my ( $class, $name, $sys, $opt ) = @_ ;
  my $self = {} ;
  bless $self, $class ;
  $self->{compiler} = $name ;
  $self->{system} = $sys ;
  $self->{optimisation} = $opt ;
  my $sub="init_$name" ;
  $self->{CPP} = "cpp" ;
  $self->{CPPFLAGS} = () ;
  $self->{VARS} = {} ;
  $self->{DYNAMIC_LIB_EXT} = ".so" ;
  $self->$sub ( )
    or Defs::Error( "Unrecognized compiler system $name" ) ; 
  $self ;
}

sub init_cxx {
  my ( $self ) = @_ ;
  $self->{opt} = "-O" ;
  $self->{dbg} = "-g" ;
  $self->{cc} = "cc" ;
  $self->{CXX} = "cxx" ;
  $self->add_macrodef( "-D__USE_STD_IOSTREAM" ) ;
  my $opt = $self->{optimisation} ;
  $self->{CXXFLAGS} = "$self->{$opt} " ;
  $self->{FFLAGS} = $self->{CXXFLAGS} ;
  $self->{LDFLAGS} = $self->{CXXFLAGS} ;
  $self->{FC} = "f77" ;
  $self->{LD} = "cxx" ;
  $self->{LIBS} = [ "Ufor", "Futil", "for", "m" ] ;
  if( -f "/usr/lib/libfl.a" ) 
    { push @{$self->{LIBS}}, "fl" }
  else
    { push @{$self->{LIBS}}, "l" }
  $self->{LIBSDIR} = ["/usr/shlib"] ;
  $self->{AR} = "cxx -shared -o \$@ " ;
  $self->{ARFLAGS} = "$self->{$opt}" ;

  $self->{CPP} = "cxx" ;
  $self->{LDLIBSSO} = "" ;

  1;
}

sub init_gcc {
  my ( $self ) = @_ ;
  my $gcc_version = `gcc --version` or Defs::Error "Unable to find gcc compiler" ;
  $self->{opt} = "-O3" ;
  $self->{dbg} = "-g" ;
  $self->{cc} = "gcc" ;
  $self->{CXX} = "g++" ;
  my $opt = $self->{optimisation} ;
  $self->{CXXFLAGS} = "$self->{$opt} -fPIC" ;
  $self->{FFLAGS} = $self->{CXXFLAGS} ;
  $self->{LDFLAGS} = $self->{CXXFLAGS} ;
  my $warnings=" -Wall -Wno-ctor-dtor-privacy -pedantic -W ".
    "-Wcast-qual -Wwrite-strings -Wconversion -Winline" ;
  if( $gcc_version =~ /2.95/ ) { 
    $warnings .= " -Wno-unused" 
  }
  else {
    $warnings .= " -Wshadow -Wno-unused-parameter" 
  }
  $self->{CXXFLAGS} .= $warnings ;
  $self->{FC} = "g77" ;
  $self->{LD} = "g++" ;
  $self->{LIBS} = [ "g2c", "m" ] ;
  if( -f "/usr/lib/libfl.a" ) 
    { push @{$self->{LIBS}}, "fl" }
  else
    { push @{$self->{LIBS}}, "l" }
  $self->{LIBSDIR} = [] ;
  $self->{AR} = "g++ -shared -o \$@ " ;
  $self->{ARFLAGS} = "$self->{$opt}" ;

  $self->{CPP} = "gcc" ;
  $self->{LDLIBSSO} = "" ;
  1;
}

sub finalize_cxx {
  my ( $self ) = @_ ;
  my $rflags = " -rpath " ;
  map { $self->{LIBRARIES} .= " -L$_ $rflags $_" ; } @{$self->{LIBSDIR}} ;
  if( $self->{profiling} ) {
    my $opt = " -pg" ;
    $self->{CXXFLAGS} .= $opt ;
    $self->{LDFLAGS} .= $opt ;
  }
  $self->{MKDEP} = "cxx -M -noimplicit_include \$(CPPFLAGS)" ;

}

sub finalize_gcc {
  my ( $self ) = @_ ;
  my $rflags = " -Xlinker -rpath -Xlinker " ;
  if( $self->{system} eq "SunOS" ){ 
    $rflags = " -Xlinker -R -Xlinker" ; 
    $self->{AR} = "g++ -Wl,-G -o \$@ " ;
  }   
  map { $self->{LIBRARIES} .= " -L$_ $rflags $_" ; } @{$self->{LIBSDIR}} ;
  if( $self->{profiling} ) {
    my $opt = " -pg" ;
    $self->{CXXFLAGS} .= $opt ;
    $self->{LDFLAGS} .= $opt ;
  }
  if( $self->{system} =~ /CYGWIN/ ) {
  	$self->{DYNAMIC_LIB_EXT} = ".dll" ;
  	$self->{AR} = "g++ -shared -Wl,--out-implib=\$@.a -o \$(\@D)/\$(subst lib,cyg,\$(\@F)) " ;
  	$self->{ARFLAGS} = " -Wl,--export-all-symbols -Wl,--enable-auto-import -Wl,--whole-archive " ;
  	$self->{LDLIBSSO} .= " -Wl,--no-whole-archive -lg2c -lfl" ;
  }
}

sub finalize_CC {
  my ( $self ) = @_ ;
  map { $self->{LIBRARIES} .= " -L$_ -R$_" ; }  @{$self->{LIBSDIR}} ;
  $self->{LDLIBSSO} .= " -lCstd -lCrun" ;
  if( $self->{profiling} ) {
    my $opt = " -xprofile=tcov" ;
    $self->{CXXFLAGS} .= $opt ;
    $self->{LDFLAGS} .= $opt ;
  }
}


sub init_CC {
  my ( $self ) = @_ ;
  $self->{opt} = "-fast -xtarget=ultra2" ;
  $self->{dbg} = "-g" ;

  $self->{cc} = "cc" ;
  $self->{CXX} = "CC" ;
  my $opt = $self->{optimisation} ;
  $self->{CXXFLAGS} = "$self->{$opt} -KPIC" ;
  $self->{FC} = "f77" ;
  $self->{LD} = "CC" ;
  $self->{LDFLAGS} = $self->{CXXFLAGS} ;
  $self->{FFLAGS} = $self->{CXXFLAGS} ;
  $self->{LIBS} = [ "F77", "M77", "sunmath", "l" ] ;
  $self->{LIBSDIR} = [] ;
  $self->{AR} = "CC -G -o \$@ " ;
  $self->{ARFLAGS} = " $self->{$opt} " ;
  1;
}

sub finalize {
  my $self = shift @_ ;
  $self->{MKDEP} = "gcc -M \$(CPPFLAGS)" ;
  my $sys = $self->{compiler} ;
  my $sub="finalize_$sys" ;
  $self->$sub() ; 
  map { $self->{LIBRARIES} .= " -l$_" ; } @{$self->{LIBS}} ;
}

sub run {
  my ( $self, $makefile, $target ) = @_ ;
  $self->finalize() ;
  my @args = ( "-f $makefile", 
	       "CCC=\"$self->{CXX}\"",
	       "AR_CMD=\"$self->{AR}\"",
	       "FC=\"$self->{FC}\"",
	       "LIB=\"$self->{LIBRARIES}\"",
	       "CCFLAGS=\"$self->{CXXFLAGS}\"",
	       "LDFLAGS=\"$self->{LDFLAGS}\"",
	       "CPPFLAGS=\"@{$self->{CPPFLAGS}}\"",
	       "FFLAGS=\"$self->{FFLAGS}\"",
	       "archiveName=$target",
	       "$target" ) ;
  print "arg is @args\n" ;
  system "make @args" ;
}

sub print {
  my ( $self ) = @_ ;
  $self->finalize() ;
  foreach (keys %{$self->{VARS}} ) {
    print "$_ = $self->{VARS}{$_}\n" ;
  }
  print "CC       = $self->{cc}\n" ;
  print "CXX      = $self->{CXX}\n" ;
  print "CXXFLAGS = $self->{CXXFLAGS}\n" ;
  print "CPP      = $self->{CPP}\n" ;
  print "CPPFLAGS = ".shift @{$self->{CPPFLAGS}} ;
  map { print "\\\n  $_" ; } @{$self->{CPPFLAGS}} ; print"\n" ;
  print "FC       = $self->{FC}\n" ;
  print "FFLAGS   = $self->{FFLAGS}\n" ;
  print "LD       = $self->{LD}\n" ;
  print "LDFLAGS  = $self->{LDFLAGS}\n" ;
  print "LDLIBS   = $self->{LIBRARIES}\n" ;
  print "LDLIBSSO = $self->{LDLIBSSO}\n" ;
  print "AR       = $self->{AR}\n" ;
  print "ARFLAGS  = $self->{ARFLAGS}\n" ;
  print "DYNAMIC_LIB_EXT = $self->{DYNAMIC_LIB_EXT}\n" ;
  print "MKDEP.c  = $self->{MKDEP}\n" ;
  print "MKDEP.cc = $self->{MKDEP}\n" ;
}

sub add_library {
  my ( $self, @libs ) = @_ ;
  unshift @{$self->{LIBS}}, @libs ;
}

sub add_include {
  my ( $self, @libs ) = @_ ;
  map { $self->add_macrodef(" -I$_" ) ; } @libs ;
}

sub add_library_path {
  my ( $self, @libs ) = @_ ;
  push @{$self->{LIBSDIR}}, @libs ;
}

sub add_macrodef {
  my ( $self, $macro ) = @_ ;
  push @{$self->{CPPFLAGS}}, $macro ;
}

sub do_profiling {
  my ( $self ) = @_ ;
  $self->{profiling}=1 ;
}

sub add_var {
  my ( $self, $name, $val ) = @_ ;
  $self->{VARS}{$name}=$val ;
}

1;
