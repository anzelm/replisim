#!/usr/bin/env perl 
#===============================================================================
#
#        FILE: post_correlate.pl   
#        USAGE: ./post_correlate.pl   max_time saturation_time
#		Correlate RepliSim result with experimental replicated tracks.
#
################################################################################

$tav_limit = $ARGV[0];
$tav_step = 50;
$tav_step = $ARGV[1] if ($ARGV[1]) ;

while (<STDIN>)
{
	#	uncomment if one chromosome only
	#	next unless ( /^22/ );
	chomp;
	@a = split;
	$tav = $a[2];
	$tst = ( $a[2] > $tav_step ) ? 1 : 0 ;
	$trx = $a[4];

	next if ($tav < 0 );
	next if ($tav > $tav_limit );
	push @tav, $tav;
	push @trx, $trx;
	push @tst, $tst;
	++ $nrws;
}

$cr = corr ( \@tav, \@trx );
$cs = corr ( \@tst, \@trx );

print $cr, "\t", $cs, "\t($tav_limit $nrws)\n";


exit;


sub corr {
    my ($aref, $bref) = @_;
    my $n; my @zna, @znb, @za, @zb;
    if ( $#$aref != $#$bref ) { return "NA" };
    $n = $#$bref + 1;
    @za = @{$aref}; @zb = @{$bref};
    @zna = normalize_sigma (@za);
    @znb = normalize_sigma (@zb);
    my $ss = 0;
    for (0..$n-1) { $ss += $zna[$_] * $znb[$_]; };
    $ss /= $n; return $ss;
};

sub normalize_sigma {
        my @num = @_; my $sum = 0;
        for (@num) { $sum += $_ ; };
        $sum /= (1.0 + $#num);
        for (0..$#num) { $num[$_] -= $sum ; } ;
        my $ssq = 0;
        for (@num) { $ssq += $_ * $_ ; };
        $ssq /= (1.0 + $#num); $ssq = sqrt($ssq);
        for (0..$#num) { $num[$_] /= $ssq ; } ;
        return @num;
};



