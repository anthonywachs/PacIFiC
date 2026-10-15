#!/usr/bin/perl

#-------------------------------------------------------------------
# Appel : MakeAdd.pl rep_obj fic_in fic_out
# 
# Ce script est utilise pour ajouter le chemein menant aux objets au
# fichier cree par "makedepend".
# Il insere le nom du repertoire au debut de toutes les lignes 
# commencant par les caracteres :
#  - a -> z,
#  - A -> Z,
#  - et _ 
#
# - rep_obj :repertoire ou se trouvent les objets,
# - fic_in : fichier en lecture,
# - fic_out : fichier en ecriture.
#-------------------------------------------------------------------
use File::Basename;

$Add = $ARGV[0];
$FicIn = $ARGV[1];
$FicOut = $ARGV[2];

open(IN,"$FicIn");
open(OUT,">$FicOut");
while (<IN>)
{
    chop;
    ($Cible,$Target) = split(/:/,$_);
     $Cible =~ s/ //g;
		#s/^(\w+.*)/$Add$1/g;
    if ( $Cible =~ /.*\.o$/) { 
      my $baseCible = '';
      $baseCible = basename($Cible);
      print OUT "$Add$baseCible :$Target\n" ;
    }
}
close (IN);
close (OUT);


