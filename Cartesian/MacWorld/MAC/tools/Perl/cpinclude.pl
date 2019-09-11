package main ;
use strict ;
use Defs;
use Arch;
use File::Basename ;
use File::Find ;
use File::Spec ;
use Getopt::Long;
use Pod::Usage;
use Util;

sub cpinc{
  my ($ref_hash,$root,$verbosity) = @_;
  my %packageslist = %{$ref_hash} ;
  my $dest = "$root/include";
  my $cmd = '';
  
  for (keys %packageslist) 
  {
    if ( $verbosity ) {print "Copying from @{$packageslist{$_}} to $dest\n";}
  
    # Copy *.hh
    $cmd = "cp @{$packageslist{$_}}/include/*.hh $dest/"; 
    system( $cmd );
  
    # Copy *.icc and *.h if they exist
    my $found_icc = 0;
    my $found_h = 0;
    my $dir = "@{$packageslist{$_}}/include";
    my $filename = '';
    my $dirs = '';
    my $suffix = '';  
    opendir(DIR, $dir) or die $!;
    while (my $file = readdir(DIR)) {
      ($filename, $dirs, $suffix) = fileparse("$file", qr/\.[^.]*/);
      if ( $suffix eq ".icc" ) { $found_icc += 1 ;}
      if ( $suffix eq ".h" ) { $found_h += 1 ;}
    }
    closedir(DIR);
    if ( $found_icc ) {
      $cmd = "cp @{$packageslist{$_}}/include/*.icc $dest/"; 
      system( $cmd );
    }
    if ( $found_h ) {
      $cmd = "cp @{$packageslist{$_}}/include/*.h $dest/"; 
      system( $cmd );
    }    
  }  
}


#-------------------------------------------------------------------------
# MAIN start here
#-------------------------------------------------------------------------
my $man = '' ;
my $H = '' ;
my $verbosity = '' ;
my $src_mac = '' ;
my $src_api = '' ;
my $root=$Defs::mac_home ;
my $cmd = '';
my %macPackages ;
my %extAPIPackages ;

my $result = GetOptions (man => \$man
			 , help => \$H
			 , verbose =>  \$verbosity
			 , mMAC =>  \$src_mac
			 , mAPI => \$src_api ) ;

pod2usage( -exitstatus => 0, -verbose => 2) if $man;
pod2usage( -verbose => 1 ) if $H  ;


# Clean include directory
$cmd = "rm -rf $root/include/*";
if ( $verbosity ) { print "$cmd\n";}
system( $cmd ); 


# Copy include files *.hh and *.icc from all packages to include
if( $src_mac ) { 
  %macPackages=Defs::package($root) ;
} else {
  Defs::Error( "MAC packages not specified" ) ;  
}
cpinc( \%macPackages, $root, $verbosity );

if( $src_api ) { 
  %extAPIPackages=Defs::externalAPIpackage($root) ;  
  cpinc( \%extAPIPackages, $root, $verbosity );
} 

exit;


##
# POD Documentation
#
__END__

=head1 NAME

cpinclude - copy include files from all subdirectories to include

=head1 SYNOPSIS

=over

mac cpinclude [-help|-man]

mac cpinclude [options...] 

=back

=head1 DESCRIPTION

For admin purposes only, copy include files from all subdirectories to include
for later link a MAC-based application to the library

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-mMAC>

Add the directories of the MAC repository to the list
of paths to copy include files. 

=item B<-mAPI>

Add the directories of the ExternalAPI repository to the list
of paths to copy include files.

=back

=head1 EXAMPLES

=over

=item C<mac cpinclude -mMAC -mAPI>

=back

=cut

