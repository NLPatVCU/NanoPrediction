#!usr/bin/perl

=head1 NAME

multipleRegressionPrediction.pl - This program predicts the value of the given dependent variable 
(mean particle size by default) for the given value of the given independent variable (reaction temperature by default).

=head1 SYNOPSIS
 
This program takes as input a value for the independent variable (reaction temperature by default), 
and based on several regressions calculated from the provided data, estimates the resulting dependent 
variable (mean particle size by default).

=head1 USAGE

Usage: multipleRegressionPrediction.pl [OPTIONS]

=head1 Optional Arguments:

Displays the quick summary of the program options.

=head2 --dependent STRING

This option is used to set the dependent variable for predictions. 
The value should be a column in the spreadsheet with everything but underscores, letters,
and whitespace removed and whitespace replaced wtih underscores.

=head2 --independent STRING

This option is used to set the independent variable for predictions. 
The value should be a column in the spreadsheet with everything but underscores, letters,
and whitespace removed and whitespace replaced wtih underscores. To add multiple
independent variables, seperate them with commas and no spaces.

=head2 --file FILE

This option is used to set the location of the file where all of the data 
is located. By default, this is set to Cu Syntheses_2017_NAL.csv.

=head2 --showEquations

Displays the equations of best fit.

=head2 --version

Displays the version information.

=head2 --help

Displays the quick summary of program options.

=head1 OUTPUT

Outputs the expected value of the dependent variable for each given value of the independent variable

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
use Statistics::Regression;

eval(GetOptions("version", "help", "dependent=s", "independent=s", "file=s", "showEquations")) or die ("Please check the above mentioned option(s).\n");

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
my $fileLocation = $opt_file || "Cu Syntheses_2017_NAL.csv";
my $independentStr = $opt_independent || "Reaction_Temperature_C";
my @independent = split(",", $independentStr);

use v5.14;
use Data::Dumper;
use Scalar::Util qw(looks_like_number);

say "Dependent Variable: " . $dependent;
say "Independent Variables: " . join(", ", @independent);

#reads the text file
open(my $file, $fileLocation);
my $line_num = 0;
my @keys;       #represents a list of all columns in the table
my @table;      #represents the table of all values
my %independentHash;
my %methods = ();
for(my $i = 0; $i < (scalar @independent); $i++){
    $independentHash{$independent[$i]} = 1;
}
while(<$file>){
    my @line = split("\t", $_);
    if($line_num++ == 0){
        @keys = @line;
        if(-e "masses.txt"){
            push @keys, ("Solvent_Mass", "Precursor_Mass", "Reducing_Agent_Mass", "Stabilizer_Mass", "Precursor_to_Reducing_Agent_Ratio");
        }
        my ($in, $dep) = (0, "");
        for(my $i = 0; $i < (scalar @keys); $i++){
            #remove special characters and spaces
            $keys[$i] =~ s/[^ a-zA-Z_0-9#]//g;
            chomp($keys[$i]);
            $keys[$i] =~ s/ +/_/g;
            $keys[$i] =~ s/(_+$)//g;
            if($independentHash{$keys[$i]}){
                $in++;
            } if($keys[$i] eq $dependent) {
                $dep = "1";
            }
        }
        for(my $i = 0; $i < (scalar @keys); $i++){
            splice @keys, $i--, 1 if length $keys[$i] < 1;
        }
        if($in != scalar @independent){
            die("One or more independent variables do not exist!\n Options are " . Dumper(\@keys) . ".\n");
        } elsif (!$dep) {
            die("Dependent variable $dependent does not exist!\n Options are " . Dumper(\@keys) . ".\n");
        }
        
    } else {
        my $row = {};
        for(my $i = 0; $i < (scalar @keys); $i++){
            $$row{$keys[$i]} = $line[$i];
            if($keys[$i] eq "Method" and length(($line[$i])) > 0){
                $methods{$line[$i]}++;
            }
        }
        push(@table, $row);
    }
}

my %molecularWeights = ();
if(-e "masses.txt"){
    open(my $weights, "masses.txt");
    while(<$weights>){
        chomp;
        my @row = split("\t", $_);
        $molecularWeights{$row[0]} = $row[1];
    }
    foreach my $row (@table){
        if(defined $row->{"Solvent"} and defined $molecularWeights{$row->{"Solvent"}}){
            $row->{"Solvent_Mass"} = $molecularWeights{$row->{"Solvent"}};
        }
        if(defined $row->{"Precursor"} and defined $molecularWeights{$row->{"Precursor"}}){
            $row->{"Precursor_Mass"} = $molecularWeights{$row->{"Precursor"}};
        }
        if(defined $row->{"Reducing_agent"} and defined $molecularWeights{$row->{"Reducing_agent"}}){
            $row->{"Reducing_Agent_Mass"} = $molecularWeights{$row->{"Reducing_agent"}};
        }
        if(defined $row->{"Stabilizer"} and defined $molecularWeights{$row->{"Stabilizer"}}){
            $row->{"Stabilizer_Mass"} = $molecularWeights{$row->{"Stabilizer"}};
        }
        if(defined $row->{"Precursor_Conc_mM"} and defined $row->{"Precursor_Volume_mL"} and defined $row->{"Precursor_Mass"} and (not defined $row->{"Precursor_Mass_mg"} or not looks_like_number($row->{"Precursor_Mass_mg"}))){
            $row->{"Precursor_Mass_mg"} = calcMass(calcMoles($row->{"Precursor_Conc_mM"}, $row->{"Precursor_Volume_mL"}), $row->{"Precursor_Mass"})
                if looks_like_number($row->{"Precursor_Conc_mM"}) and looks_like_number($row->{"Precursor_Volume_mL"}) and looks_like_number($row->{"Precursor_Mass"});
        }
        if(defined $row->{"Reducing_Agent_Conc_mM"} and defined $row->{"Reducing_Agent_Volume_mL"} and defined $row->{"Reducing_Agent_Mass"} and (not defined $row->{"Reducing_Agent_Mass_mg"} or not looks_like_number($row->{"Reducing_Agent_Mass_mg"}))){
            $row->{"Reducing_Agent_Mass_mg"} = calcMass(calcMoles($row->{"Reducing_Agent_Conc_mM"}, $row->{"Reducing_Agent_Volume_mL"}), $row->{"Reducing_Agent_Mass"})
                if looks_like_number($row->{"Reducing_Agent_Conc_mM"}) and looks_like_number($row->{"Reducing_Agent_Volume_mL"}) and looks_like_number($row->{"Reducing_Agent_Mass"});
        }
        if(defined $row->{"Precursor_Conc_mM"} and defined $row->{"Precursor_Mass_mg"} and defined $row->{"Precursor_Mass"} and (not defined $row->{"Precursor_Volume_mL"} or not looks_like_number($row->{"Precursor_Volume_mL"}))){
            $row->{"Precursor_Volume_mL"} = calcMass((1000 * $row->{"Precursor_Mass_mg"} / $row->{"Precursor_Mass"}), $row->{"Precursor_Conc_mM"})
                if looks_like_number($row->{"Precursor_Conc_mM"}) and looks_like_number($row->{"Precursor_Mass_mg"}) and looks_like_number($row->{"Precursor_Mass"});
        }
        if(defined $row->{"Reducing_Agent_Conc_mM"} and defined $row->{"Reducing_Agent_Mass_mg"} and defined $row->{"Reducing_Agent_Mass"} and (not defined $row->{"Reducing_Agent_Volume_mL"} or not looks_like_number($row->{"Reducing_Agent_Volume_mL"}))){
            $row->{"Reducing_Agent_Volume_mL"} = calcMass((1000 * $row->{"Reducing_Agent_Mass_mg"} / $row->{"Reducing_Agent_Mass"}), $row->{"Reducing_Agent_Conc_mM"})
                if looks_like_number($row->{"Reducing_Agent_Conc_mM"}) and looks_like_number($row->{"Reducing_Agent_Mass_mg"}) and looks_like_number($row->{"Reducing_Agent_Mass"});
        }
        if(defined $row->{"Reducing_Agent_Conc_mM"} and defined $row->{"Reducing_Agent_Volume_mL"} and defined $row->{"Reducing_Agent_Mass"} and (not defined $row->{"Reducing_Agent_Mass_mg"} or not looks_like_number($row->{"Reducing_Agent_Mass_mg"}))){
            $row->{"Reducing_Agent_Mass_mg"} = calcMass(calcMoles($row->{"Reducing_Agent_Conc_mM"}, $row->{"Reducing_Agent_Volume_mL"}), $row->{"Reducing_Agent_Mass"})
                if looks_like_number($row->{"Reducing_Agent_Conc_mM"}) and looks_like_number($row->{"Reducing_Agent_Volume_mL"}) and looks_like_number($row->{"Reducing_Agent_Mass"});
        }
        if(defined $row->{"Precursor_Mass_mg"} and defined $row->{"Reducing_Agent_Mass_mg"}){
            $row->{"Precursor_to_Reducing_Agent_Ratio"} = $row->{"Precursor_Mass_mg"} / $row->{"Reducing_Agent_Mass_mg"} if looks_like_number($row->{"Precursor_Mass_mg"}) and looks_like_number($row->{"Reducing_Agent_Mass_mg"})
        }
    }
}
print "What method is the reaction? ";
my $method;
while(!($methods{chomp($method = <STDIN>)})){
    last if $methods{$method} || $method eq "All";
    say "Options are \"All\" to select all methods or " . Dumper(\{%methods});
    print "What method is the reaction? ";
}
my $points = 0;
#Find the independent variable values and dependent variable values.
for(my $i = 0; $i < (scalar @table); $i++){
    my $valid = 1;
    foreach my $variable (@independent) {
        if(!$table[$i]->{$variable} || !looks_like_number($table[$i]->{$variable})){
            $valid = 0;
        }
    }
    my $args = [1];
    foreach my $variable (@independent) {
        push($args, $table[$i]->{$variable});
    }
    if($valid and $table[$i]->{$dependent} and looks_like_number($table[$i]->{$dependent}) and ($table[$i]->{"Method"} eq $method || $method eq "All")){
        addPoint($args, $table[$i]->{$dependent});
    }
    if($table[$i]->{$dependent} and looks_like_number($table[$i]->{$dependent}) and ($table[$i]->{"Method"} eq $method || $method eq "All")){
        addPossiblePoint($args, $table[$i]->{$dependent});
    }
}
$dependent =~ s/_/ /g;
for(my $i = 0; $i < scalar @independent; $i++){
    $independent[$i] =~ s/_/ /g;
}

say "$points valid points available.";
unshift @independent, "Constant";
if(($points) <= 1){
    die("Not enough numerical points available!\n");
}

#Prediction calculation
my %results;
$results{"Logarithmic"} = logarithmicRegression();
$results{"Exponential"} = exponentialRegression();
$results{"Linear"} = linearRegression();
$results{"Quadratic"} = quadraticRegression();

#evaluates different types of equations
my %evaluators;
$evaluators{"Logarithmic"} = \&logarithmic;
$evaluators{"Exponential"} = \&exponential;
$evaluators{"Quadratic"} = \&quadratic;
$evaluators{"Linear"} = \&linear;
my $combinedEquation = genEquation();
#takes user input and dependent variable

LOOP: while(1){
    my @x = ();
    foreach my $xval (@independent) {
        if($xval ne "Constant"){
            print "Enter " . $xval. " or type stop to stop: ";
            my $input = <STDIN>;
            chomp($input);
            if(!looks_like_number($input)){
                last LOOP;
            }
            push(@x, $input);
        }
    }
    foreach my $regression (sort { -abs($results{$a}[0]) <=> -abs($results{$b}[0]) } keys %results) {
        say "Using a " . $regression . " regression (correlation coefficient " . $results{$regression}[0] . "), estimated $dependent = " . $evaluators{$regression}->(\@x, $results{$regression});
    }
    say "Using the combined equation, estimated $dependent = " . evaluate($combinedEquation, \@x);
    say "\n";
}

#Regression code
my @xpoints;
my @ypoints;

my @possibleXs;
my @possibleYs;
sub addPoint {
    push(@xpoints, $_[0]);
    push(@ypoints, $_[1]);
    $points++;
}
sub addPossiblePoint {
    push(@possibleXs, $_[0]);
    push(@possibleYs, $_[1]);
}

#Provide the neccesary information to linearize data for regression
sub logarithmicRegression {
    my $variables = shift || \@independent;
    my @xs = @{shift || \@xpoints};
    my @ys = @{shift || \@ypoints};
    my $logarithmicRegression = Statistics::Regression->new("Logarithmic", $variables);
    for(my $i = 0; $i < $points; $i++){
        my $xvalues = [1];
        for(my $j = 1; $j < scalar @{$xs[$i]}; $j++){
            $$xvalues[$j] = $xs[$i][$j] > 0 ? log($xs[$i][$j]) : "NaN";
        }
        $logarithmicRegression->include($ys[$i], $xvalues);
    }
    my @coefficients = $logarithmicRegression->theta();
    unshift @coefficients, $logarithmicRegression->rsq();
    my $equationText = "Logarithmic:  \ty =";
    for(my $i = 2; $i < scalar @coefficients; $i++){
        $equationText .= " " . $coefficients[$i] . " * ln(x" . ($i - 2) . ") +";
    }
    $equationText .= " " . $coefficients[1] . "";
    say $equationText . (" " x (110 - length $equationText, 5)[110 - length $equationText < 5]) . "r^2 = $coefficients[0]" if ($showEquations and $variables eq \@independent);
    return \@coefficients;
}
sub exponentialRegression {
    #linearizes data by taking the log of the y value
    my $variables = shift || \@independent;
    my @xs = @{shift || \@xpoints};
    my @ys = @{shift || \@ypoints};
    my $exponentialRegression = Statistics::Regression->new("Exponential", $variables);
    for(my $i = 0; $i < $points; $i++){
        $exponentialRegression->include(log($ys[$i]), $xs[$i]);
    }
    my @coefficients = $exponentialRegression->theta();
    unshift @coefficients, $exponentialRegression->rsq();
    my $equationText = "Exponential:  \ty = e^(";
    for(my $i = 2; $i < scalar @coefficients; $i++){
        $equationText .= " " . $coefficients[$i] . " * x" . ($i - 2) . " +";
    }
    $equationText .= " " . $coefficients[1] . ")";
    say $equationText . (" " x (110 - length $equationText, 5)[110 - length $equationText < 5]) . "r^2 = $coefficients[0]" if ($showEquations and $variables eq \@independent);
    return \@coefficients;
}
sub quadraticRegression { #Assumes all values, both x and y, are non-negative
    my $variables = shift || \@independent;
    my @xs = @{shift || \@xpoints};
    my @ys = @{shift || \@ypoints};

    my $quadraticRegression = Statistics::Regression->new("Quadratic", $variables);
    for(my $i = 0; $i < $points; $i++){
        $quadraticRegression->include(sqrt($ys[$i]), $xs[$i]);
    }
    my @coefficients = $quadraticRegression->theta();
    unshift @coefficients, $quadraticRegression->rsq();
    my $equationText = "Quadratic:    \ty = (";
    for(my $i = 2; $i < scalar @coefficients; $i++){
        $equationText .= " " . $coefficients[$i] . " * x" . ($i - 2) . " +";
    }
    $equationText .= " " . $coefficients[1] . ")^2";
    say $equationText . (" " x (110 - length $equationText, 5)[110 - length $equationText < 5]) . "r^2 = $coefficients[0]" if ($showEquations and $variables eq \@independent);
    return \@coefficients;
}
sub linearRegression {
    my $variables = shift || \@independent;
    my @xs = @{shift || \@xpoints};
    my @ys = @{shift || \@ypoints};
    my $linearRegression = Statistics::Regression->new("Linear", $variables);
    for(my $i = 0; $i < $points; $i++){
        $linearRegression->include($ys[$i], $xs[$i]);
    }
    my @coefficients = $linearRegression->theta();
    unshift @coefficients, $linearRegression->rsq();
    my $equationText = "Linear:       \ty =";
    for(my $i = 2; $i < scalar @coefficients; $i++){
        $equationText .= " " . $coefficients[$i] . " * x" . ($i - 2) . " +";
    }
    $equationText .= " " . $coefficients[1];
    say $equationText . (" " x (110 - length $equationText, 5)[110 - length $equationText < 5]) . "r^2 = $coefficients[0]" if ($showEquations and $variables eq \@independent);
    return \@coefficients;
}

#Delinearize the regression
sub quadratic {
    my $x = shift;
    my $m = shift;
    my $sum = 0;
    for(my $i = 2; $i < scalar @{$m}; $i++){
        $sum += $$m[$i] * $$x[$i - 2];
    }
    return sqrt($sum + $$m[1]);
}
sub logarithmic {
    my $x = shift;
    my $m = shift;
    my $sum = 0;
    for(my $i = 2; $i < scalar @{$m}; $i++){
        $sum += $$m[$i] * log($$x[$i - 2]);
    }
    return $sum + $$m[1];
}
sub exponential {
    my $x = shift;
    my $m = shift;
    my $sum = 0;
    for(my $i = 2; $i < scalar @{$m}; $i++){
        $sum += $$m[$i] * $$x[$i - 2];
    }
    return exp($sum + $$m[1]);
}

#Find the value of a point on a line
sub linear {
    my $x = shift;
    my $m = shift;
    my $sum = 0;
    for(my $i = 2; $i < scalar @{$m}; $i++){
        $sum += $$m[$i] * $$x[$i - 2];
    }
    return $sum + $$m[1];
}

#Finds the best regression type
sub regressionType {
    my $independentVar = shift;
    my @independents = ();
    for(my $i = 0; $i < scalar @independent; $i++){
        if($independent[$i] eq $independentVar){
            for(my $j = 0; $j < scalar @possibleXs; $j++){
                $independents[$j] = $possibleXs[$j][$i];
            }
        }
    }
    my @xs = ();
    my @ys = ();
    for(my $i = 0; $i < (scalar @independents); $i++){
        if($independents[$i] and looks_like_number($independents[$i])){
            push(@xs, [1, $independents[$i]]);
            push(@ys, $ypoints[$i]);
        }
    }
    my %options;
    $options{"Logarithmic"} = logarithmicRegression(["Constant", $independentVar], \@xs, \@ypoints);
    $options{"Exponential"} = exponentialRegression(["Constant", $independentVar], \@xs, \@ypoints);
    $options{"Linear"} = linearRegression(["Constant", $independentVar], \@xs, \@ypoints);
    $options{"Quadratic"} = quadraticRegression(["Constant", $independentVar], \@xs, \@ypoints);
    my $regression = (sort { -abs($options{$a}[0]) <=> -abs($options{$b}[0]) } keys %options)[0];
    return [$regression, $options{$regression}];
}

#Given a list of values, returns a list indicating what percentage of the sum of all values in the list each element was.
sub percentage {
    my @values = @{$_[0]};
    my $sum = 0;
    foreach my $value (@values){
        $sum += $value;
    }
    for(my $i = 0; $i < scalar @values; $i++){
        $values[$i] /= $sum;
    }
    return \@values;
}

sub genEquation {
    my %variables = ();
    my $displayEquation = "y = ";
    foreach my $xvar (@independent) {
        next if($xvar eq "Constant");
        $variables{$xvar} = regressionType($xvar);
    }
    my @weights = ();
    foreach my $xvar (@independent) {
        next if($xvar eq "Constant");
        push @weights, $variables{$xvar}[1][0]
    }
    @weights = @{percentage(\@weights)};
    my $equation = "";
    my $i = 0;
    foreach my $xvar (@independent) {
        next if $xvar eq "Constant";
        $equation .= "(" . $variables{$xvar}[0] . " " . ($variables{$xvar}[1][2] * $weights[$i]) . " " . ($variables{$xvar}[1][1] * $weights[$i]) . ") ";
        if($displayEquation ne "y = "){
            $displayEquation .= " + ";
        }
        if($variables{$xvar}[0] eq "Exponential"){
            $displayEquation .= "e^(" . ($variables{$xvar}[1][2] * $weights[$i]) . " * x$i + " . ($variables{$xvar}[1][1] * $weights[$i]) . ")";
        } elsif($variables{$xvar}[0] eq "Logarithmic"){
            $displayEquation .= ($variables{$xvar}[1][2] * $weights[$i]) . " * ln(x$i) + " . ($variables{$xvar}[1][1] * $weights[$i]);
        } elsif($variables{$xvar}[0] eq "Quadratic"){
            $displayEquation .= "(" . ($variables{$xvar}[1][2] * $weights[$i]) . " * x$i + " . ($variables{$xvar}[1][1] * $weights[$i]) . ")^2";
        } elsif($variables{$xvar}[0] eq "Linear"){
            $displayEquation .= ($variables{$xvar}[1][2] * $weights[$i]) . " * x$i + " . ($variables{$xvar}[1][1] * $weights[$i++]);
        }
    }
    my $error = 0;
    my $meanError = 0;
    my $yMean = 0;
    foreach my $y (@ypoints){ $yMean += $y; }
    $yMean /= (scalar @ypoints);
    for(my $i = 0; $i < $points; $i++){
        shift $xpoints[$i];
        $error += ($ypoints[$i] - evaluate($equation, $xpoints[$i])) ** 2;
        unshift $xpoints[$i], 1;
        $meanError += ($ypoints[$i] - $yMean) ** 2;
    }
    say "Combined equation r^2 = " . (1 - $error / $meanError);
    if($showEquations){
        say $displayEquation;
    }
    return $equation;
}

sub evaluate {
    my @equation = split(/\) \(/, substr($_[0], 1, length($_[0]) - 3));
    my $xdata = $_[1];
    my $index = 0;
    my $sum = 0;
    foreach my $piece (@equation){
        my @pieceBreakdown = split(/ /, $piece);
        my $coef = $pieceBreakdown[1];
        my $intercept = $pieceBreakdown[2];
        if($pieceBreakdown[0] eq "Logarithmic"){
            $sum += $coef * log($$xdata[$index++] > 0 ? $$xdata[$index - 1] : 1) + $intercept;
        } elsif($pieceBreakdown[0] eq "Exponential"){
            $sum += exp($coef * $$xdata[$index++] + $intercept);
        } elsif($pieceBreakdown[0] eq "Linear"){
            $sum += $coef * $$xdata[$index++] + $intercept;
        } elsif($pieceBreakdown[0] eq "Quadratic"){
            $sum += sqrt($coef * $$xdata[$index++] + $intercept);
        }
    }
    return $sum;
}

sub calcMoles { #in moles
    my $molarity = $_[0] / 1000; #in M
    my $volume = $_[1] / 1000; #in L
    return ($molarity) * ($volume);
}

sub calcMass {
    my $moles = $_[0]; #in moles
    my $molecularWeight = $_[1] * 1000; #in g
    return $moles * $molecularWeight;
}

sub calcVolume { #in L
    my $moles = $_[0]; #in moles
    my $molarity = $_[1] / 1000; #in M
    return ($moles / $molarity) * 1000;
}

sub calcMolecularWeight {
    my $moles = $_[0];
    my $mass = $_[1] / 1000;
    return $mass / $moles;
}

##############################################################################            

#  function to output the version number                                                   

##############################################################################              

sub showVersion {
    say '$Id: multipleRegressionPrediction.pl,v 1.0 2017/06/21 12:23 jack Exp $';
}


##############################################################################              

#  function to output "ask for help" message when user's goofed                              

##############################################################################               

sub askHelp {
    print STDERR "Type multipleRegressionPrediction.pl --help for help.\n";
}


##############################################################################                 

#  function to output help messages for this program                                         

##############################################################################                

sub showHelp() {
    print "This program takes as input a value for each independent \n";
    print "variable (reaction temperature by default), and based on\n"; 
    print "several regressions calculated from the provided data, estimates\n";
    print "the resulting value of the dependent variable (mean particle size by default).\n\n";

    print "Usage: multipleRegressionPrediction.pl [OPTIONS] \n\n";

    print "OPTIONS:\n\n";

    print "--dependent DEPENDENT       This option takes the string DEPENDENT\n";
    print "                            and uses the column that this string denotes\n";
    print "                            as the dependent variable. Type this in with \n";
    print "                            special characters removed and spaces as underscores\n\n";

    print "--independent INDEPENDENT   This option takes the string INDEPENDENT\n";
    print "                            and uses the column(s) that this string denotes\n";
    print "                            as the independent variable(s). Type this in with \n";
    print "                            special characters removed and spaces as underscores\n";
    print "                            To add multiple independent variables, separate them\n";
    print "                            with commas and no spaces\n\n";

    print "--file FILE                 This option takes the file FILE and tries to open\n";
    print "                            a file at that location. The data in the specified\n";
    print "                            file is used by the program for its data. If none is\n";
    print "                            specified, Cu Syntheses_2017_NAL.csv is used\n\n";

    print "--showEquations             Shows the equations of best fit\n\n";

    print "--version                   Prints the version number\n\n";

    print "--help                      Prints this help message.\n\n";
}