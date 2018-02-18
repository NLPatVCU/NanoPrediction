#!usr/bin/perl

=head1 NAME

methodCharting.pl - This program draws a graph where each reaction is a point colored according to the
method used in the reaction. Both the independent, or x, variable and the dependent, or y, variable can be
specified for the gra[h. 

=head1 SYNOPSIS
 
This program draws a graph where each reaction is a point colored according to the
method used in the reaction. Both the independent, or x, variable and the dependent, or y, variable can be
specified for the graph. Once the program is run, a GUI appears with the graph. By mousing over any item
on the legend, all other colors are set to 30% opacity. Additionally, the graph can be saved by pressing the 
"Save to File" button and the colors can be changed using the "Reset Colors" button.

=head1 USAGE

Usage: methodCharting.pl [OPTIONS]

=head1 Optional Arguments:

Displays the quick summary of the program options.

=head2 --dependent STRING

This option is used to set the dependent or y variable on the resulting
graph. If no dependent is specified, mean particle diameter is used.

=head2 --color STRING

This option is used to set the column used to color each point on the resulting
graph. If no color is specified, method is used.

=head2 --independent STRING

This option is used to set the independent or x variable on the resulting
graph. If no independent is specified, reaction temperature is used.

=head2 --file STRING

This option is used to set the file from which data will be read.

=head2 --outliers

If this option is used, then outliers will not be excluded from the graph.

=head2 --version

Displays the version information.

=head2 --help

Displays the quick summary of program options.

=head1 SYSTEM REQUIREMENTS

=over

=item * Perl (version 5.8.5 or better) - http://www.perl.org

=back

=head1 CONTACT US

    If you have trouble installing and excecuting methodCharting.pl, 
    please contact us at
    
    btmcinnes at vcu dot edu.

=head1 Author

 Jack N. McDowell

=head1 COPYRIGHT

Copyright (c) 2015

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

eval(GetOptions("version", "help", "dependent=s", "independent=s", "outliers", "file=s", "color=s")) or die ("Please check the above mentioned option(s).\n");

#   if help is defined, print out help
if( defined $opt_help ) {
    $opt_help = 1;
    &showHelp;
    exit;
}

# if version is requested, show version
if( defined $opt_version ) {
    $opt_version = 1;
    &showVersion();
    exit;
}

my $showEquations = $opt_showEquations;
my $dependent = $opt_dependent || "Mean_Particle_Diameter_nm_TEM";
my $independent = $opt_independent || "Reaction_Temperature_C";
my $outlier_str = defined $opt_outliers ? " " . $opt_outliers : "";
my $file_name = $opt_file || "Cu Syntheses_2017_NAL.csv";
my $color = $opt_color || "Method";

use v5.14;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);
use warnings;

#reads the text file
open(my $file, $file_name);
my $line_num = 0;
my @keys;       #represents a list of all columns in the table
my @table;      #represents the table of all values
my %methods = ();
while(<$file>){
    $_ =~ s/[^\t,\. a-zA-Z0-9_#]//g;
    my @line = split("\t", $_);
    if($line_num++ == 0){
        @keys = @line;
        my ($in, $dep, $col) = ("", "", "");
        for(my $i = 0; $i < (scalar @keys); $i++){
            #remove special characters and spaces
            $keys[$i] =~ s/[^ a-zA-Z_#]//g;
            chomp($keys[$i]);
            $keys[$i] =~ s/ +/_/g;
            $keys[$i] =~ s/(_+$)//g;
            if($keys[$i] eq $independent){
                $in = "1";
            } if($keys[$i] eq $dependent) {
                $dep = "1";
            } if($keys[$i] eq $color) {
                $col = "1";
            }
        }
        for(my $i = 0; $i < (scalar @keys); $i++){
            splice @keys, $i--, 1 if length $keys[$i] < 1;
        }
        if(!$in){
            die("Independent variable $independent does not exist!\n Options are " . Dumper(\@keys) . ".\n");
        } elsif (!$dep) {
            die("Dependent variable $dependent does not exist!\n Options are " . Dumper(\@keys) . ".\n");
        } elsif (!$col) {
            die("Color variable $color does not exist!\n Options are " . Dumper(\@keys) . ".\n");
        }
    } else {
        my $row = {};
        for(my $i = 0; $i < (scalar @keys); $i++){
            $$row{$keys[$i]} = $line[$i];
            $$row{$keys[$i]} =~ s/[\(|\)|']//;
            if($keys[$i] eq $color and length(($line[$i])) > 2){
                $methods{$line[$i]}++;
            }
        }
        push(@table, $row);
    }
}
#Find the independent variable values and dependent variable values.
my @xpoints;
my @ypoints;
my @methods;
for(my $i = 0; $i < (scalar @table); $i++){
    if(looks_like_number($table[$i]->{$independent}) and looks_like_number($table[$i]->{$dependent}) and $table[$i]->{$color}){
        push @xpoints, $table[$i]->{$independent} + 0;
        push @ypoints, $table[$i]->{$dependent} + 0;
        push @methods, $table[$i]->{$color};
    }
}
my $xs = join(",", @xpoints);
my $ys = join(",", @ypoints);
my $methodsList = join(",", @methods);
$methodsList =~ s/ /_/g;
system("java -jar ./MethodCharting.jar " . $methodsList . " " . $xs . " " . $ys . " " . $independent . " " . $dependent . " " . $outlier_str);

##############################################################################            

#  function to output the version number                                                   

##############################################################################              

sub showVersion {
    say '$Id: methodCharting.pl,v 1.0 2017/06/21 12:23 jack Exp $';
}

##############################################################################              

#  function to output "ask for help" message when user's goofed                              

##############################################################################               

sub askHelp {
    print STDERR "Type methodCharting.pl --help for help.\n";
}

##############################################################################                 

#  function to output help messages for this program                                         

##############################################################################                

sub showHelp() {
    print "methodCharting.pl - This program draws a graph where each reaction is a\n";
    print "point colored according to the method used in the reaction. Both the independent,\n";
    print "or x, variable and the dependent, or y, variable can be specified for the graph.\n\n";

    print "Usage: methodCharting.pl [OPTIONS] \n\n";

    print "OPTIONS:\n\n";

    print "--dependent DEPENDENT       This option takes the string DEPENDENT\n";
    print "                            and uses the column that this string denotes\n";
    print "                            as the dependent variable in the generated graph.\n";
    print "                            If none is specified, mean particle diameter is used\n\n";

    print "--file FILENAME             This option takes the string FILENAME and looks in\n";
    print "                            the folder methodCharting.pl is in for a file with\n";
    print "                            that name. It then uses the data from that file for the graph.\n";
    print "                            If none is specified, Cu Syntheses_2017_NAL.csv is used\n\n";

    print "--independent INDEPENDENT   This option takes the string INDEPENDENT\n";
    print "                            and uses the column that this string denotes\n";
    print "                            as the independent variable in the generated graph.\n";
    print "                            If none is specified, reaction temperature is used\n\n";

    print "--color COLOR               This option takes the string COLOR\n";
    print "                            and uses the column that this string denotes\n";
    print "                            as the coloring variable in the generated graph.\n";
    print "                            If none is specified, method is used\n\n";

    print "---outliers                 When this option is used, outliers will not be excluded\n";
    print "                            from the graph.\n\n";

    print "--version                   Prints the version number\n\n";

    print "--help                      Prints this help message.\n\n";
}
