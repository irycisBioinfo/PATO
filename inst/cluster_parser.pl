#!/usr/bin/perl

open(A,"dnds.fasta");

$prev ="";
$i=0;
$j=0;
$k=1;

$folder = "folder_0";
mkdir $folder;

if(-e "Correspondences.tsv")
{
	unlink("Correspondences.tsv");
}
open(T,">Correspondences.tsv");

while($l = <A>)
{

	if($j == 2000)
	{
		$folder = "folder_".$k;
		$j=0;
		$k++;
		mkdir $folder
	}

	if($l =~ />/ && $prev =~ />/)
	{
	  {

  		$file = "family_".$i;
  		$i++;
  		$j++;
  		close O;
  		open(O,">$folder/$file");

  		$index = $l;
  		chomp $index;
  		$index =~ s/>//;
  		print T "$index\t$folder/$file\n";
	  }

	}

	$prev = $l;


	if($i ==0)   ##For the first line
	{
		next;
	}
	print O $l;

}
