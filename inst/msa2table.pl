#!/usr/bin/perl

open(A,$ARGV[0]);
@txt = <A>;
close A;



for($i =0; $i<scalar(@txt)-1; $i = $i+1)
{
	if(
	{
		$l1 = $txt[$i];
		$l2 = $txt[$i+1];
		chomp $l1;
		chomp $l2;
	
		$l1 =~ s/>//;
		print "$l1\t$l2\n";
	}
		
}
