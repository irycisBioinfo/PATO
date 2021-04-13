#!/usr/bin/perl

$output = $ARGV[0];
open(A,$ARGV[1]);

open(OUT,">$output.gff");

while($l = <A>)
{
	if($l =~ "#FASTA")
	{
		close OUT;
		open(OUT,">$output.fasta");
	}else{
		
		print OUT $l;
	}
}
