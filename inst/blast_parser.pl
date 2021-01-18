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


	if($l =~ /(Query_1\s+1\s+)/)
	{
		
		$offset = length($1)
	}
	if($l =~ /Subject_(\d+)\s+/)
	{
		$tmp = "$header[$1-1]";
		chomp $tmp;
		print $tmp;
		
		chomp $l;
		$seq = substr($l,$offset,length($l));
		$seq =~ s/\s+\d+//;
		print "\t$seq\n";
		
	}
	
}
