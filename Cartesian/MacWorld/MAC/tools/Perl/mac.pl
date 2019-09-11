# MAIN

package main ;

use diagnostics ;
use strict ;
use Pod::Usage;
use File::Spec ;
use Defs ;

## Command recovering

pod2usage( -verbose => 1 ) unless $ARGV[0] ;
if ( $ARGV[0] eq "-h" || $ARGV[0] eq "-help" ) {
  pod2usage( -verbose => 1 ) ;
  shift @ARGV ;
}
if ( $ARGV[0] eq "-m" || $ARGV[0] eq "-man" ) {
  pod2usage( -verbose => 2 ) ;
  shift @ARGV ;
}
my $verbose = 0 ;
if ( $ARGV[0] eq "-v" || $ARGV[0] eq "-verbose" ) {
  $verbose=1;
  shift @ARGV;
}

my $cmd = shift @ARGV ;
pod2usage( -verbose => 0 ) if !$cmd ;

## --Run script

my @script_cmd = Defs::script( "$cmd" ) ;
print "Running " . "@script_cmd" . " ". "@ARGV" ."\n" if $verbose;
#exec ( @script_cmd, @ARGV ) or Defs::Error "Unable to execute $cmd" ;
# Peut-être... à voir avec le boss...
my $exit_statut = system ( @script_cmd, @ARGV ) ;
exit $exit_statut

##
# POD Documentation
#
__END__

=head1 NAME

mac - management of MAC-based applications

=head1 SYNOPSIS

mac [-help|-man]

mac [options...] command [command-options-and-arguments]

command is either :

  depend         generation of makefiles
  build          generation of executables, objects and libraries
  run            execution of an application
  arch           display of the available compiler architecture

[command-options-and-arguments] : depends on the chosen command.

=head1 DESCRIPTION  

C<mac> is a collection of utilities devoted to the management
of MAC-based applications.

=head1 OPTIONS

The following options are those of C<mac> itself. See the
manpages of the various commands for their specific options.

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=back

=cut
