package Util ;

use strict ;
use diagnostics ;
use File::Spec ;

BEGIN {
  require Exporter;
  
  use vars qw($VERSION @ISA @EXPORT @EXPORT_OK);
  $VERSION = 1.0;
  @ISA = qw(Exporter);
  @EXPORT = qw( get_field );
}

sub get_field{ 
  my ( $file, $pattern ) = @_ ;
  open HAND, $file or Defs::Error( "Unable to find $file" ) ;
  my $result ;
  foreach (<HAND>) { if(/$pattern/) { $result = $1 ; last ; } }
  $result ;
}

sub absolute_pathname {
  my( $dir ) = @_ ;
  my @dirs = File::Spec->splitdir( File::Spec->rel2abs( $dir ) ) ;
  my $i = 1 ;
  do {
    if( $dirs[$i] eq '..' ) {
      splice( @dirs, $i-1, 2 ) ;
      --$i ;
    } 
    elsif( $dirs[$i] eq '.' ) {
      splice( @dirs, $i, 1 ) ;
    }
    else {
      ++$i ;
    }
  } while( $i <= $#dirs ) ;
  my $result = File::Spec->join( @dirs ) ;
  return $result ;
}

sub common_root {
  my( @dir_elms ) = @_ ;
  my( @resu ) = () ;

  #---debug map{ print"dir_elms $_\n" } @dir_elms ;
  my $keep_going = 1 ;
  while( $keep_going ) {

    my $suivant = undef ;
    for( my $i = 0 ; $i <= $#dir_elms ; $i++ ) {
      my( @dir ) = File::Spec->splitdir( $dir_elms[$i] ) ;
      #---debug print "splitted dir : " ;
      #---debug for ( @dir ) { print "$_;" }
      my( $nn ) = shift( @dir ) ;
      #---debug print " ... tentative : $nn\n" ;
      $dir_elms[$i] = File::Spec->join( @dir ) ;
      if( $i == 0 ) {
        $suivant = $nn ;
      } elsif( defined( $suivant ) && ( $nn ne $suivant ) ) {
        $keep_going = 0 ;
        $suivant = undef ;
      }
      if( $#dir == 0 ) {
        $keep_going = 0 ;
      }
    }
    if( defined( $suivant ) ) {
       push @resu, $suivant ;
       #---debug map{ print" $_\n" } @resu ;
    }

  }
  return( File::Spec->join( @resu ) ) ;
}


sub win_relative_path {
  my ($file,$project) = @_ ;
  my $result = "" ;
  if( $file =~ /^\$\(MAC_HOME\)/ ) {
    $result = $file ;
	$result =~ s/\//\\/g ;
  } else {
    my @filetab = () ;
    if( Arch::is_posix() == 1 ) {
      @filetab = split /\//, absolute_pathname( $file ) ;
    } else {
      @filetab = split /\\/, absolute_pathname( $file ) ;
    }  
    my @projecttab = () ;
    if( Arch::is_posix() == 1 ) {
      @projecttab = split /\//, absolute_pathname($project) ;
    } else {
      @projecttab = split /\\/, absolute_pathname($project) ;
    }
    shift @filetab ;
    shift @projecttab ;
    my $ok = 1 ;
    while( $ok ) {
      if( scalar( @filetab ) && scalar( @projecttab ) ) {
        if( "$filetab[0]" eq "$projecttab[0]" ) {
	      shift @filetab ;
	      shift @projecttab ;
        } else { last ; }
      } else { last ; }
    } 
    shift  @projecttab ;
    map { $result .= "..\\" ; } ( @projecttab ) ;
    my $least = pop @filetab ;
    map { $result .= "$_\\" ; } ( @filetab ) ;
    if ( $least ) {
      $result .= $least ;
    } else {
      $result = "." ;
    }
  }
  $result ;
}

1;
