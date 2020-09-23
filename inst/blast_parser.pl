#!/usr/bin/perl


open(A,$ARGV[0]);
@txt = <A>;
close A;

open(B,$ARGV[1]);
@header = <B>;
close B;

foreach $l (@txt)
{
	chomp $l;


	if($l =~ /Query/)
	{
		next;
	}
	if($l =~ /Subject_(\d+)\s+/)
	{
		$tmp = "$header[$1-1]";
		chomp $tmp;
		print $tmp;
	}
	if($l =~ /Subject.{13}(.*)\s+\d+/)
	{
		$aln =$1;
		$aln =~ s/\s/-/g;
		print "\t$aln\n";
	}
}
