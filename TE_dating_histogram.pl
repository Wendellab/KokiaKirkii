#!/usr/bin/perl

local $\ = "\n";
local $, = "\t";

# Defaluts from repeatexplorer
my $PID = 90;
my $OVL = 0;
my $SCOV = 0;
my $LCOV = 55;

my $prev;
my $histo = {};
while(<>){
    # Parse mgblast output
    next if m/^#/;
    my ($id, $long, $short, $pid, $score) = ([], {}, {}, 0, 0);
    ($id->[0], @$long{qw/length start stop/}, $id->[1], @$short{qw/length start stop/}, $pid, $score) = split;

    # Skip self-hits
    next if $id->[0] eq $id->[1];

    # filter minimum overlap
    $_->{overlap} = abs($_->{start} - $_->{stop}) +1 foreach ($long, $short);
    next unless (($long->{overlap} >= $OVL || $short->{overlap} >= $OVL) && $pid >= $PID);

    # filter minimum coverage
    $_->{cov} = $_->{overlap}*100.0/$_->{length} foreach ($long, $short);
    ($short, $long) = ($long, $short) if($short->{length} > $long->{length}); 
    next unless ($short->{cov} >= $SCOV && $long->{cov} >=$LCOV);

    # Skip lower HSPs from same Hit
    $id = join("\t", sort @$id);
    next if $prev && $prev eq $id;
    $prev = $id;

    $histo->{$pid}++;
}    

#Print histogram
foreach my $pid (sort {$a <=> $b} keys %$histo){
    print $pid, $histo->{$pid};
}
