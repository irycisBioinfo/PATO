#!/usr/bin/perl

open(A,$ARGV[0]);
@txt = <A>;
close A;


for($i=0; $i<scalar(@txt);$i++)
{
	if($txt[$i] =~ /^>/ && $txt[$i+1] =~ /^>/)
	{
		$num++;
		open(OUT,">core_$num.aln");
	}else{
		print OUT $txt[$i];
	}
}
