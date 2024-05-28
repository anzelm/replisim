#!/usr/bin/env perl 
#===============================================================================
#
#         FILE: sort_origins.pl
#
#        USAGE: ./sort_origins.pl  
#        preprocess the origins file 
#

while (<>)
{
	$oo = $_;
	@a = split;
	$w = $a[0] * 1000000000 + $a[1] + $a[2] ;
	$val { $oo } = $w;
}

for ( sort { $val{$a} <=> $val{$b} } keys %val )
{
	print;
}
