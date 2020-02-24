#!/usr/bin/perl -w
use strict;
use IO::File;

my $fname="summary_hbm.txt";
#Columns of summary file are
#Event	Gene	Sig. RBPs	Sig. Gene Expression	Sig. Sex	Tissue	Dimorphic


#This script extracts the data that we will need for
#figure 4c (expression vs RBP coefficients in mammary tissue).


my $out=new IO::File(">mt.txt") or die "$!";
my $lv= new IO::File(">lv.txt") or die "$!";


my $fh=new IO::File($fname) or die "$!";
my $header=<$fh>; # skip the header line
while (<$fh>) {
    chomp;
    my @F=split(m/\t/);
    my $gene=$F[1];
    my $rbpList=$F[2];
    my $expr=$F[3];
    my $sex=$F[4];
    my $tissue=$F[5];
   
   
    my $dimorphic=$F[6];
    print "$_\n";
    my $rbp=parseRbpCoefficients($rbpList);
    next unless (defined $rbp);
    next unless  ($dimorphic eq "Yes");
    if ($tissue =~ m/Mammary/) {
        print $out "$rbp\t$expr\n";
    } elsif ($tissue =~ m/Ventricle/) {
	print $lv "$rbp\t$expr\n";
    }
}

$out->close();
$lv->close();



sub parseRbpCoefficients {
    my $str=shift;
    #print "STR $str\n";
    my @B=split(",",$str);
    my $sum=0;
    return undef unless (defined $str && length($str)>0);
    foreach my $b(@B) {
	if ($b =~ m/\w+\(([\d\.\-]+)\)/) {
	    $sum += abs($1);
	} else {
	    die "CRAP WHAT IS THIS $b\n";
	}
    }
    return $sum;
}
