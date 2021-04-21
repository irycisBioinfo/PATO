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


	if($l =~ /(Query_1\s+1\s+)(\S+)\s+(\d+)/)
	{

		$offset = length($1);

		$seq_length = length($2);
	}
	if($l =~ /Subject_(\d+)\s+/)
	{
		$tmp = "$header[$1-1]";
		chomp $tmp;
		print $tmp;

		chomp $l;
		$seq = substr($l,$offset,$seq_length);
		$seq =~ s/\s/-/g;
		print "\n$seq\n";

	}
}
