# MAIN

package main ;
use strict ;
use Getopt::Long;
use Pod::Usage;
use Defs;
use POSIX ;
use Arch;
my $verbosity = '' ;

##
# Command line Parsing

my $man = '' ;
my $H = '' ;
my $target = '' ;
my $machine = POSIX::uname() ; chomp $machine ;
my $make = "make" ;
my $link_mode = "cc" ;
my $makefile = '' ;

my @with=();
my @without=();


#prepend the MAC_BUILD to ARGV ...
unshift(@ARGV, split(/ /,$ENV{MAC_BUILD})) if ($ENV{MAC_BUILD});

my $arch_exe = "exe" ;
if( Arch::is_posix() == 0 ) {
  $arch_exe .= ".exe" ;
} 
my $result = GetOptions ('exe' => sub { $target = $arch_exe ; }
			 , 'man' => \$man
			 , 'help' => \$H
			 , 'make=s' => \$make
			 , 'archive=s' => \$target
			 , 'object=s' =>  \$target
                         , 'link_mode=s' => \$link_mode
			 , 'verbose' =>  \$verbosity
			 , 'with=s' => \@with
			 , 'without=s' => \@without
                         , 'makefile=s' => \$makefile
			) ;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
pod2usage(-verbose => 1 ) if $H  ;
if( scalar @ARGV != 1 || !$result || !$target ) {
  pod2usage( -verbose => 0 ) ;  
}
@with = split(/,/,join(',',@with));
@without = split(/,/,join(',',@without));

my $lib_path = shift @ARGV ;

if( $makefile eq "" ) {
  $makefile = File::Spec->catfile( $lib_path, "Makefile" ) ;
}
if( ! -f $makefile ) {
  my $mesg = "Unknown file : $makefile\n\n" ;
  $mesg .= "Use \"mac --depend\" to build the makefile \"Makefile\" \n" ;
  $mesg .= "in the directory \"$lib_path\" prior to the\n" ;
  $mesg .= "execution of \"mac --depend\"" ;
  Defs::Error( $mesg ) ;
}

if( $verbosity ) {
  print "Compile $target on $machine \n" ;
  print "\tin $lib_path \n" ;
  print "\twith $makefile\n" ;
  print "\tenabling extra : " . join(',',@with) . "\n" if (@with);
  print "\tdisabling extra : " . join(',',@without) . "\n" if (@without);
}

Defs::Error( "Unable to find Makefile option : $makefile" ) unless -f $makefile ;
##
my @makecmd = ( $make, "-f", $makefile, "TARGET=$target", "LINK_MODE=$link_mode" ) ;
map {push @makecmd, "WITH_".uc($_)."=1"} @with;
map {push @makecmd, "WITH_".uc($_)."=0"} @without;
if( $verbosity ) {
   print "Executing command : \n   " ;
   for ( @makecmd ) { print " $_" ; }
   print "\n"
}
exec(@makecmd);
##
# POD Documentation
#
__END__

=head1 NAME

build - generation of executables, objects and libraries related to MAC-based applications

=head1 SYNOPSIS

mac build [-help|-man]

mac build [options...] -exe F<bindir>

mac build [options...] -object F<filename.o> F<bindir>

mac build [options...] -archive F<archname.so> F<bindir>

=head1 DESCRIPTION

The generation of executables, objects and libraries associated to
a particular MAC-based application can be simply performed by using
the GNU C<make> utility. The first step consists in building
a suitable makefile with the C<mac depend> utility. The second
step consists in running GNU C<make>  with suitable arguments,
a task whose responsiblity is assigned to C<mac build>.

Essentially, GNU C<make> is invoked with the
following instructions:

=over

=item 1.

read the makefile called F<Makefile>, located in the directory F<bindir>
(unless otherwise specified with the C<-makefile> options);

=item 2.

update the target determined by the mutually exclusive instructions
C<-exe,-object,-archive>;

=item 3.

use a specific set of make options that are determined by
the calling options of C<mac build>.

=back

=head1 OPTIONS

=over

=item B<-h, -help>

Print a brief help message and exit.

=item B<-m, -man>

Print the manual page and exit.

=item B<-v, -verbose>

Execute verbosely.

=item B<-e, -E, -exe>

Ask C<make> to
update the target defined by 
the executable of the MAC-based application.

=item B<-o, -O, -object> F<filename.o>

Ask C<make> to update the target
defined by the module object F<filename.o>
(typically, C<filename> is the name of one of
the classes making up
the considered application).

=item B<-a, -A, -archive> F<archname.ext>

Ask C<make> to update the target
defined by the archive F<archname.ext>
(the associated command depends on
the extension F<.ext>).

=item B<-link_mode> link_mode

Set the linker:
    link_mode=cc : use the c++ linker (default)
    link_mode=c  : use the c linker
    link_mode=f  : use the fortran linker

=item B<-make> makename

Use the command of name C<makename> as 
the make command (default: C<make>).
This option is typically used for systems
on which the GNU C<make> utility is
called C<gmake>.

=item B<-makefile> makefile

Set the makefile used for all tasks of compilation,
assembly and linking. Default is: F<bindir>/Makefile.

=item B<-with> EXTlist

Notify that the current application DOES require the
packages of C<EXTlist> so that the archives
associated to those packages
should be added to the list of files to link.
C<EXTlist> is a comma separated list denoting external APIs
that might be used by some components of MAC.
This option is meaningful only for the targets
of the options C<-exe> and C<-archive>.
Note that C<mac depend> defines some defaults in the
generated makefile for all possible external libraries.

=item B<-without> EXTlist

Notify that the current application DOES NOT require the
packages of C<EXTlist> so that the archives
associated to those packages
should NOT be added to the list of files to link.
C<EXTlist> is a comma separated list denoting external APIs
that might be used by some components of MAC.
This option is meaningful only for the targets
of the options C<-exe> and C<-archive>.
Note that C<mac depend> defines some defaults in the
generated makefile for all possible external libraries.

=back

=head1 ARGUMENT

=over

=item B< F<bindir> >

Directory containing the makefile generated by C<mac depend>,
called F<Makefile>. Any file produced by the execution
of C<mac build> (that is: when updating a target of F<Makefile>)
will be created in that directory (objects, libraries, executable,
dependency files ...). The executable will be called F<exe>.

=back

=head1 EXAMPLES

=over

=item C<mac depend -l mac1 dbg bin/dbg .>

=item C<mac build -exe bin/dbg>

The current application is made of all the header and source files
located in any subdirectory of the working directory. The
generated F<Makefile>, the executable F<exe> and all the files
created during the compiling process will be located in the
subdirectory F<bin/dbg> of the working directory.

=item C<mac build -without petsc,opengl -exe lib>

All the files generated for and during the compilation process
are located in the subdirectory lib of the working directory
(and possibly in the working directory itself).
Linking will be performed without the archives associated
to the PETSc and OpenGL libraries.

=item C<mac build -object myclass.o /usr/smith/appli/bin>

Update the target F<myclass.o> of the makefile F<Makefile>
located in the directory F</usr/smith/appli/bin>, that
is: build object file F<myclass.o> in that directory.

=back

=head1 ENVIRONMENT

It is possible to store arguments and options, overwritable by the command
line arguments, in the environment variable MAC_BUILD.

=cut
