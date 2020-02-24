#!/usr/bin/perl -w
use strict;
use IO::File;

my $fname="summary_hbm.txt";
#Columns of summary file are
#Event	Gene	Sig. RBPs	Sig. Gene Expression	Sig. Sex	Tissue	Dimorphic

my $out=new IO::File(">rbpCoeff.txt") or die "$!";
my $out2=new IO::File(">exprCoeff.txt") or die "$!";
my $out3 = new IO::File(">sexCoeff.txt") or die "$!";

my $fh=new IO::File($fname) or die "$!";
my $header=<$fh>; # skip the header line
while (<$fh>) {
    chomp;
    my @F=split(m/\t/);
    my $gene=$F[1];
    my $rbpList=$F[2];
    my $expr=$F[3];
    my $sex=$F[4];
    my $dimorphic=$F[6];
    my @rbp=parseRbpCoefficients($rbpList);
    next unless (@rbp);
    if ($dimorphic eq "Yes") {
	foreach my $r (@rbp) {
	    print $out "$r\tbiased\n";
	}
	print $out2 "$expr\tbiased\n";
	print $out3 "$sex\tbiased\n";
    } elsif ($dimorphic eq "No"){
	foreach my $r (@rbp) {
	    print $out "$r\tnonbiased\n";
	}
	print $out2 "$expr\tnonbiased\n";
	print $out3 "$sex\tnonbiased\n";
    } else {    
	die "Could not parse dimorphic status \"$dimorphic\"";
    }
}

$out->close();
$out2->close();
$out3->close();


sub parseRbpCoefficients {
    my $str=shift;
    #print "STR $str\n";
    my @B=split(",",$str);
    my @ret;
    return @ret unless (defined $str && length($str)>0);
    foreach my $b(@B) {
	if ($b =~ m/\w+\(([\d\.\-]+)\)/) {
	    push(@ret,$1);
	} else {
	    die "CRAP WHAT IS THIS $b\n";
	}
    }
    return @ret;
}
