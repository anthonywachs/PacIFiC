package Arch ;

use strict ;
use diagnostics ;
use File::Spec ;
use POSIX ;

my $mak_dir = "./etc";
my $prefix = "";
my $suffix = ".mak";
my $verbose=0;
my $mac_arch_file = 'arch_file.cfg';
my $mac_arch_dir = $ENV{MAC_ARCH_DIR} ;
my $has_mac_arch_dir = ( $mac_arch_dir && -d $mac_arch_dir );
my $unspecified = "undef" ;

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub architecture_name {
  my $compiler = shift ;

  my @xx = ( undef, $unspecified, $unspecified ) ;
  if( $verbose ) { $verbose = 2 ; }
  my @arch_1 = Arch::query_best_arch( $compiler ) ;
  exit(1) unless defined $arch_1[0] ;
  if( $verbose ) { $verbose = 1 ; }
  Arch::set_prefix('extra-');
  my @arch_2 = Arch::query_best_arch( $compiler ) ;
  exit(1) unless defined $arch_2[0] ;
  Defs::Error( "inconsistent query" ) unless( $arch_2[3] eq $arch_1[3] ) ;
  if( ( ! defined $arch_2[0] ) || ( $arch_1[0] < $arch_2[0] ) ) {
    for( my $i=1 ; $i<=2 ; ++$i ) {
      if( $arch_1[$i] ne $unspecified ) {
        $xx[$i] = $arch_1[$i] ;
      } elsif( $arch_2[$i] ne $unspecified ) {
        $xx[$i] = $arch_2[$i] ;
      }
    }
  } else {
    for( my $i=1 ; $i<=2 ; ++$i ) {
      if( $arch_2[$i] ne $unspecified ) {
        $xx[$i] = $arch_2[$i] ;
      } elsif( $arch_1[$i] ne $unspecified ) {
        $xx[$i] = $arch_1[$i] ;
      }
    }
  }
  if( $arch_1[3] ne "" ) {
    $xx[1] = $arch_1[3] ;
  }
  my $arch_fin = $xx[1]."-".$xx[2] ;

  return $arch_fin ;
} ;

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub query_best_arch {
  my $compiler = shift || "gcc";
  my ($sysname, $realhostname, $release) = POSIX::uname();
  $release  =~ s/[-(].*$//; #remove training comments on the release
  $sysname  =~ s/[-(].*$//; #remove training comments on the release
  my $hostname = $realhostname ;
  my $archdeb = resolve( $realhostname ) ;
  if( $archdeb ) {
    $hostname = $archdeb ;
  }
  if( $verbose == 2 ) {
    print STDERR "Hostname: $hostname " ;
    if( $archdeb ) {
      print STDERR "(substitution for $realhostname)" ;
    }
    print STDERR "\n" ;
  }
  my @checks1 = ( $hostname
		, $hostname
		, $sysname."_".$release
		, $sysname
		, $unspecified
		, $sysname."_".$release
		, $sysname
	        );
  my @checks2 = ( $compiler
		, $unspecified
		, $compiler
		, $compiler
		, $compiler
		, $unspecified
		, $unspecified
	        );
  if( $verbose == 2) { print STDERR "Compiler: $compiler\n" ; }
  if( $verbose ) {
    if( ! $prefix ) {
      print STDERR "Architecture-Makefile searched in:\n" ;
    } else {
      print STDERR "Extra-Makefile searched in:\n" ;
    }
    print STDERR "   $mac_arch_dir\n" if ( $has_mac_arch_dir ) ;
    print STDERR "   $mak_dir\n";
  }
  my $i_resu = undef ;
  for( my $i=0 ; $i<=$#checks1 ; $i++ ) {
    my $ff = get_makefile( $checks1[$i], $checks2[$i] ) ;
    if( -f $ff ) {
      if( ! defined $i_resu ) { $i_resu = $i } ;
      if( $verbose ) {
        print STDERR "\t*  $checks1[$i]-$checks2[$i]\t$ff\n" ;
      }
    } else {
      if( $verbose ) {
        print STDERR "\t   $checks1[$i]-$checks2[$i]\n" ;
      }
    }
  }
  if( defined $i_resu ) {
     return ( $i_resu, $checks1[$i_resu], $checks2[$i_resu], $archdeb ) ;
  }
  if ($verbose) {
    print STDERR "FAILED: unable to find an " ;
    if( ! $prefix ) {
      print STDERR "Architecture-Makefile" ;
    } else {
      print STDERR "Extra-Makefile" ;
    }
    print STDERR ".\n\n";
  }

  return;
}

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub set_verbose {$verbose=shift;}
sub get_verbose {$verbose;}

sub set_prefix {$prefix=shift;}
sub get_prefix {$prefix;}

sub set_suffix {$suffix=shift;}
sub get_suffix {$suffix;}

sub set_mak_dir {$mak_dir=shift;}
sub get_mak_dir {$mak_dir;}

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub get_makefile {
  my ($check1,$check2) = @_ ;
  my $xxx ;
  if( ($check1 ne $unspecified) && ($check2 ne $unspecified) ) {
    $xxx .= $check1."-".$check2 ;
  } elsif( $check1 ne $unspecified ) {
    $xxx .= $check1 ;
  } elsif( $check2  ne $unspecified ) {
    $xxx .= $check2 ;
  }
  if ( $has_mac_arch_dir ) {
    my $dd = File::Spec->catfile($mac_arch_dir, $prefix . $xxx . $suffix);
    return $dd if (-f $dd );
  }
  File::Spec->catfile($mak_dir, $prefix . $xxx . $suffix)
}

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub resolve {
  my $host = shift ;
  # on parcourt les répertoires $mac_arch_dir et $mak_dir à la recherche d'un 
  #   fichier $mac_arch_file
  # une fois trouvé on cherche une substitution du host donné en paramètre 
  #   dans ce fichier
  # le fichier $mac_arch_file est un fichier sur 2 colonnes séparées par un ou 
  #   plusieurs blancs, les commentaires débutent par #:
  #	1 : expression régulière mappant des hostnames
  #	2 : nom du host équivalent (alias)
  #	..: le reste est du commentaire
  #	une fois l'alias du host trouvé => on s'en va
  #exemple de fichier
  #sinux1 pinux    # sinux1 vers pinux
  #sinux\d+ sinux  # les autres noeuds sinux vers sinux
  #pinux\d+ pinux  # les noeuds pinux vers pinux

  my @dirs = ();
  push @dirs, $mac_arch_dir if ($has_mac_arch_dir);
  push @dirs, $mak_dir;
  my $first = "yes" ;
  foreach my $dir (@dirs) {
     my $arch_file = File::Spec->catfile($dir, $mac_arch_file);
     next unless (-f $arch_file );
     if( $verbose == 2 ) {
        if( $first eq "yes" ) {
           print STDERR "Possible hostname substitution searched in:\n" ; 
        }
        print STDERR "   $arch_file\n"
     }
     $first = "no" ;
     open(AFILE,"<$arch_file");
     while(<AFILE>) {
       chomp;
       s/\#.*$//;
       my ($regexp, $alias) = split;
       return $alias if ($regexp && $host =~ /^$regexp\z/) ;
     }
     close(AFILE);
  }
  return "" ;
}

#----------------------------------------------------------------------------
#
#----------------------------------------------------------------------------
sub is_posix {
  my $result = 1;
  my ($sysname) = POSIX::uname();
  if( $sysname =~ /^[Ww][Ii][Nn]/ )
  {
    $result = 0;
  }
  return $result;
}

1;
