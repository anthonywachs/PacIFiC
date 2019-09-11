# MAC Global definitions

package Defs ;

use strict ;
use diagnostics ;
use File::Spec ;
use Text::Wrap ;
use POSIX ;

BEGIN {
  require Exporter;
  
  use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
  $VERSION = 1.0;
  @ISA = qw(Exporter);
  @EXPORT = qw( package UnitTests script_cmd bin_tool Error string_vector );
  @EXPORT_OK = qw( $mac_home $mac_path );

  use vars qw( $mac_home $mac_path ) ;

  $mac_home = $ENV{MAC_HOME} or Defs::Error( "MAC_HOME not set" ) ;
  $mac_path = File::Spec->catfile( $mac_home, "tools", "Perl" ) ;

}

sub search{
  my ($rep,$underRep ) = @_ ;
  my $dirstd = File::Spec->catfile( $mac_home, $rep ) ;
  opendir DIRHAND, $dirstd ;
  my @LIST;
  map { push @LIST, "$dirstd/$_" if -d "$dirstd/$_/$underRep" } File::Spec->no_upwards( readdir DIRHAND ) ;
  @LIST ;
}


sub package{
  my $rootdir = shift @_ ;

  my %package = ( "BASE" => [ "$mac_home/Base" ],
		  "LA"  => [ "$mac_home/LinearAlgebra" ],
		  "FV"  => [ "$mac_home/FiniteVolume" ],
		  "AV"  => [ "$mac_home/ArrayVector" ]		  
		) ;
  %package ;
  #map {print "$_ @{$package{$_}}\n" ; } keys(%package) ;
}

sub externalAPIpackage{
  my $rootdir = shift @_ ;

  my %package = ( "MPI" => [ "$mac_home/ExternalAPI/MPI" ],
		  "PETSC"  => [ "$mac_home/ExternalAPI/PETSc_3.2.0" ]	  
		) ;
  %package ;
}

sub UnitTests{
  my @res = search( "UnitTests", "tests" ) ;
  @res ;
}

## --Search for MAC_HOME/tools/Perl/<command>.pl script
sub script{
  my $tool = shift @_ ; 
  my $mac_perl=File::Spec->catfile( $Defs::mac_home, 'tools', 'Perl' ) ;
  my $abstool=File::Spec->catfile( $mac_perl, "$tool.pl" ) ;
  Defs::Error( "Unknown command : $tool" ) unless -f $abstool ;
  my @res = ( "perl", "-w", "-I", "$mac_perl", "$abstool" ) ;
  @res ;
}

sub bin_tool{
  my $tool = shift @_ ;
  my $machine = POSIX::uname() ; chomp $machine ;
  my $exe = 
    File::Spec->join( $Defs::mac_home, "tools", $tool, "lib", $machine, "exe" ) ;
  Defs::Error( "Bad tool $tool \n ($exe not found)" ) unless -x $exe ;
  $exe ;
}

sub Error{
print STDOUT "\n" ;
print STDOUT "-------------------------------------------------\n" ;
print STDOUT "|        Fatal error in \"mac\" utility          \n" ;
print STDOUT "-------------------------------------------------\n" ;
print STDOUT wrap ( "| ", "| ", "@_\n" ) ;
print STDOUT "-------------------------------------------------\n" ;
exit(1) ;
}

sub string_vector {
  my $key = shift @_ ;
  my $result = ' ' ;
  if( scalar @_ ) {
    $result = " $key=<" ;
    map { $result.= " \"$_\" " } @_ ;
    $result.="> " ;
  }
  $result ;
}

1;
