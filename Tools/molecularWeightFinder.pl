#!usr/bin/perl

=head1 NAME

multipleRegressionPrediction.pl - This program is used to quickly find the molecular weight of a list 
of given compounds.

=head1 SYNOPSIS
 
This program searched pubchem for each compound in the given list and displays the top 5 results for each.
The user can then choose to select one as correct, search google and manually enter the molecular weight, or
skip that compound.

=head1 USAGE

Usage: molecularWeightFinder.pl [OPTIONS]

=head1 Optional Arguments:

Displays the quick summary of the program options.

=head2 --outfile FILE

This option takes the file FILE and uses that
as the output file. This file will be a list
of chemicals and their molecular weight, 1 per
line, separated by tabs. 

=head2 --infile FILE
This option takes the file FILE and uses that
as the input file. This file should be a list
of chemicals, 1 per line, and nothing else.

=head2 --version

Displays the version information.

=head2 --help

Displays the quick summary of program options.

=head1 OUTPUT

Outputs a file listing the inputted chemicals and their molecular weights

=head1 SYSTEM REQUIREMENTS

=over

=item * Perl (version 5.8.5 or better) - http://www.perl.org

=back

=head1 CONTACT US

    If you have trouble installing and excecuting multipleRegressionPrediction.pl, 
    please contact us at
    
    btmcinnes at vcu dot edu.

=head1 Author

 Jack N. McDowell

=head1 COPYRIGHT

Copyright (c) 2017

 Jack N. McDowell
 Bridget T. McInnes
                     
This program is free software; you can redistribute it and/or modify it 
under the terms of the GNU General Public License as published by the Free 
Software Foundation; either version 2 of the License, or (at your option) 
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT 
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with 
this program; if not, write to:

 The Free Software Foundation, Inc.,
 59 Temple Place - Suite 330,               
 Boston, MA  02111-1307, USA.          
  
=cut

##########################################################################

                        #   CODE STARTS HERE

##########################################################################
#  reference the getOption cpan page
use Getopt::Long;
use Data::Dumper;
my $args = join(" ", @ARGV);
eval(GetOptions("version", "help", "infile=s", "outfile=s")) or die ("Please check the above mentioned option(s).\n");

#   if help is defined, print out help
if( defined $opt_help ) {
    $opt_help = 1;
    print "This program takes as input a value for the independent \n";
    print "variable (reaction temperature by default), and based on\n"; 
    print "several regressions calculated from the provided data, estimates\n";
    print "the resulting value of the dependent variable (mean particle size by default).\n\n";

    print "Usage: molecularWeightFinder.pl [OPTIONS] \n\n";

    print "OPTIONS:\n\n";

    print "--outfile FILE              This option takes the file FILE and uses that\n";
    print "                            as the output file. This file will be a list\n";
    print "                            of chemicals and their molecular weight, 1 per \n";
    print "                            line, separated by tabs. This defaults to names.txt\n\n";

    print "--infile FILE               This option takes the file FILE and uses that\n";
    print "                            as the input file. This file should be a list\n";
    print "                            of chemicals, 1 per line, and nothing else.\n";
    print "                            This defaults to masses.txt\n\n";

    print "--version                   Prints the version number\n\n";

    print "--help                      Prints this help message.\n\n";
    exit;
}

# if version is requested, show version
if( defined $opt_version ) {
    $opt_version = 1;
    print '$Id: molecularWeightFinder.pl,v 1.0 2017/06/28 10:41 AM jack Exp $';
    print "\n";
    exit;
}

use v5.14;
system("java -jar ./MolarMassFinder.jar " . $args);