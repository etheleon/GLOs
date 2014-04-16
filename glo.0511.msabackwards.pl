#!/usr/bin/env perl 

use strict;
use warnings;
use v5.10;
use lib "/export2/home/uesu/perl5/lib/perl5";

use Bio::Seq;
use Bio::SeqIO;
use autodie;

die "usage: glo.0511.msabackwards.pl KXXXXX.msa" unless $#ARGV == 0;

#Reads the fasta file
my $in  = Bio::SeqIO->new(-file => $ARGV[0] , -format => 'Fasta');

say "overlap\tcontig\tmsaposition\tcontigposition";

while (my $seq = $in->next_seq ){
my $seqid 		= $seq->display_id;
if($seqid =~ /contig/) { 
my $sequence 	= $seq->seq;
my ($overlap, $contig) = split(":", $seqid);

my @fasta = split("",$sequence);

my $msacounter = 1;
my $contigcounter = 1;

foreach(@fasta) { 
   if(!/-/) {  
    say "$overlap\t$contig\t$msacounter\t$contigcounter";
    $contigcounter++;
   }
    $msacounter++;
	}
	}
}
