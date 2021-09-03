#!/usr/bin/perl -w

# Copyright 2016-2017 Active Motif, Inc.
#
# Purpose:
# Remove duplicates in position-sorted BAM data based on Molecular IDs
#
# Requirements:
#Before mapping FASTQ data, the sequence of the molecular ID (mID) must
#be appended to FASTQ header string. This ensures that
#the resulting BAM file's QNAME field will contain the
#mID sequence. The QNAME must be separated from the mID by an
#underscore symbol.
#Before de-duping the BAM file(s), they must be coordinate sorted using
#a tool such as samtools sort.
#
# History:
#02/2016 prototype written by Kuan-Bei Chen
#06/29/2016 included strand info to remove duplicates
#07/20/2016 rewritten by Steve Stelman for reduced memory usage
#10/03/2017 added paired-end mode [stelman]
#05/13/2019 small fixes [stelman]
#
# Contact: tech_service@activemotif.com
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

use strict;
use warnings;
use Getopt::Long;

my $Usage = &get_usage;

#define default path to samtools binary.
my $samtools_path = "/home/services/Genpathway/Systems/samtools-0.1.19";

##"/ssoft/spack/humagne/v1/opt/spack/linux-rhel7-x86_S6g1_Mellanox/gcc-7.4.0/samtools-1.9-armssgpe75ofxslcuzcaobarwqgtt7o2/bin/samtools";


# Initialize command line options
my $in_dir = ".";
my $out_dir = ".";
my $output_type = "all";
my $paired_end = 0; #set to non-zero for PE mode
my $help = 0;

# Get command line options
GetOptions (
    'mode=s'=> \$output_type,
    'in=s'=> \$in_dir,
    'out=s'=> \$out_dir,
    'paired'=> \$paired_end,
    'samtools=s'=> \$samtools_path,
    'help'=> \$help,
    ) or die "$0: cannot read command line: @ARGV\n\terror: $!\n";

if ($help) {
    die "$Usage\n";
}

unless ($output_type eq "rm" || $output_type eq "mk" || $output_type eq "all") {
    die "ERROR: -mode parameter \"$output_type\" is invalid\n$Usage\n";
}

unless (-d $in_dir || -f $in_dir) {
    die "ERROR: -in parameter must be directory or .bam file\n$Usage\n";
}

unless (-x "$samtools_path/samtools") {
    die "ERROR: samtools not found in $samtools_path\n$Usage\n";
}

#create output dir if necessary.
unless (-d $out_dir) {
    mkdir $out_dir or die "Can't create output dir $out_dir\n";
}

#add samtools path to $PATH if not already there.
$ENV{'PATH'} = "$samtools_path:" . $ENV{'PATH'} unless ($ENV{'PATH'} =~ /$samtools_path/);

########
# Main #
########

if (-d $in_dir) { # if input is a directoy
    my @infiles = glob("$in_dir/*.bam");

    foreach my $in_file(@infiles) {
	my %counter;
	print "\nProcessing $in_file\n";
	rmDup3($in_file, $out_dir, $output_type, \%counter, $paired_end);
#die;
    }
} elsif (-f $in_dir) { # if input is a file
    my %counter;
    my $in_file = $in_dir;
    die "ERROR: input file is not a .bam file\n" unless ($in_file =~ /\.bam$/);
    print "\nProcessing $in_file\n";
    rmDup3($in_file, $out_dir, $output_type, \%counter, $paired_end);
}

###############
# Subroutines #
###############

sub get_usage {
    my $Usage = <<"END_MESSAGE";
    Usage:
    $0 
	   -mode [rm|mk|all]    (default: all)
	      -in [in_dir|in_file] (default: cwd)
	         -out [out_dir]       (default: cwd)
		    [-paired]            (default: off)
		       -samtools [/path/to/samtools]
		          [-help]
			  
 
END_MESSAGE

return $Usage;
}

sub rmDup3 {
    #Remove/Mark PCR duplicates in sorted BAM using Molecular IDs

    my $infile = shift;
    my $outdir = shift;
    my $output_type = shift;
    my $counter = shift;
    my $paired = shift;
    
    #initialize counters
    $counter->{total} = 0;
    $counter->{kept} = 0;
    $counter->{unmapped} = 0;
    $counter->{dedup} = 0;
    $counter->{deduppair} = 0;
    $counter->{best} = 0;
    $counter->{singleton} = 0;
    
    my $path = "./";
    my ($pre, $root) = $infile =~ /(.*\/)?([^\/]+)\.bam/;
    $path = $pre if ($pre);

    #Define output filenames
    my $out_rmdup = "${outdir}/${root}.rmdup.bam";
    my $out_mkdup = "${outdir}/${root}.mkdup.bam";
    my $out_log = "${outdir}/${root}.log";

    #open input and output pipes using samtools
    open(my $bam_in, '-|', "samtools view -h $infile") or die "Can't open pipe to samtools for $infile :: $!\n";
    open(my $bam_rmdup, "| samtools view -bS - > $out_rmdup") or die "Can't open pipe to samtools for $out_rmdup :: $!\n" if ($output_type eq "all" || $output_type eq "rm");
    open(my $bam_mkdup, "| samtools view -bS - > $out_mkdup") or die "Can't open pipe to samtools for $out_mkdup :: $!\n" if ($output_type eq "all" || $output_type eq "mk");

    open(my $log, ">$out_log") or die "Can't open log $out_log :: $!\n";

    my $prev_rname = undef; #store previous RNAME to know if we're in same group
    my $prev_pos = undef; #store previous POS to know if we're in same group
    my @sam_bucket = (); #stores group of SAM records with same RNAME+POS+MID+Strand (+Tlen for PE) to be evaluated for deduping
    my @skip_bucket = (); #stores SAM records of mate pairs to maintain in sorted BAM file immediately after those in @sam_bucket at same POS
    my %mate_is_deduped;
    
    while(my $sam_in = <$bam_in>) {
	chomp $sam_in;
	if ($sam_in =~ /^\@/) {
	    #pass SAM header to output
	    print $bam_rmdup "$sam_in\n" if ($output_type eq "all" || $output_type eq "rm");
	    print $bam_mkdup "$sam_in\n" if ($output_type eq "all" || $output_type eq "mk");
	} else {
	    #process coordinate-sorted SAM data record, paired-end
	    $counter->{total}++;
	    print "Processed $counter->{total} reads\r" if (($counter->{total} % 10000) == 0); #for progress reporting
	    my %sam; #stores one SAM record
	    &parse_sam_record_w_mid(\%sam, $sam_in, $paired);
	    
	    if ($sam{rname} eq "*") { #unmapped
		#process previous bucket
		&dedup_by_mid(\@sam_bucket, \@skip_bucket, $output_type, $counter, $bam_rmdup, $bam_mkdup, \%mate_is_deduped, $paired) if ($sam_bucket[0]);

		#process unmapped reads to output
		$counter->{kept}++;
		$counter->{unmapped}++;
		print $bam_rmdup "$sam_in\n" if ($output_type eq "all" || $output_type eq "rm");
		print $bam_mkdup "$sam_in\n" if ($output_type eq "all" || $output_type eq "mk");


		# Don't bucket these in paired-end mode
	    } elsif ($paired && ( (($sam{bitflag} & 1) == 0 ) || (($sam{bitflag} & 2) == 0 ) || ($sam{bitflag} & 16)  || ($sam{rnext} ne '=') || ($sam{pnext} < $sam{pos}) ) ) {

		if ($sam_bucket[0]) {
		    #add data to skip bucket to interpolate back into BAM later
		    push(@skip_bucket, $sam{sam});
		} else {
		    #if sam bucket is empty, just print now instead of adding to skip bucket for subsequent interpolation
		    $counter->{kept}++;

		    my @fields = split(/\t/, $sam_in);
		    $fields[1] -= 1024 if ($fields[1] & 1024);

		    print $bam_rmdup join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "rm");
		    print $bam_mkdup join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "mk");
		}


	    } elsif (! $prev_pos) { #first in bucket
		#initialize prior POS and RNAME
		$prev_pos = $sam{pos};
		$prev_rname = $sam{rname};
		
		#start new bucket
		push(@sam_bucket, $sam{sam});
	    } else {

		if ($sam{pos} == $prev_pos && $sam{rname} eq $prev_rname) { #next in bucket
		    #if ((($sam{bitflag} & 16) == 0) && ($sam{bitflag} & 32) && $sam{pos} == $prev_pos && $sam{rname} eq $prev_rname) { #next in bucket
		    #add data to bucket
		    push(@sam_bucket, $sam{sam});
		} else { #new bucket
		    #process previous bucket
		    &dedup_by_mid(\@sam_bucket, \@skip_bucket, $output_type, $counter, $bam_rmdup, $bam_mkdup, \%mate_is_deduped, $paired);

		    #reset prior POS and RNAME
		    $prev_pos = $sam{pos};
		    $prev_rname = $sam{rname};

		    #start new bucket
		    push(@sam_bucket, $sam{sam});
		}
	    }
	}
    }
    #process final bucket
    &dedup_by_mid(\@sam_bucket, \@skip_bucket, $output_type, $counter, $bam_rmdup, $bam_mkdup, \%mate_is_deduped, $paired) if ($sam_bucket[0]);
    
    print "Processed $counter->{total} reads\n";# if (($counter->{total} % 10000) == 0); #for progress reporting

    print $log "COUNTERS:\n";
    print $log "$counter->{total}\tTotal reads\n";
    print $log "$counter->{kept}\tTotal reads retained\n";
    print $log "$counter->{unmapped}\tUnmapped reads retained\n";
    print $log "$counter->{best}\tNon-singltons retained\n";
    print $log "$counter->{singleton}\tSingltons retained\n";
    print $log "$counter->{dedup}\tTotal deduped reads\n";
    print $log "$counter->{deduppair}\tTotal deduped mate pairs\n";
    print $log ($counter->{total} - $counter->{kept}) . "\tCalculated deduped reads\n";
}

sub dedup_by_mid {
    #remove PCR duplicates from bucket of SAM data by Molecular ID
    
    my $sam_bucket = shift; # array ref with SAM data
    my $skip_bucket = shift; # array ref with SAM data
    my $output_type = shift;
    my $counter = shift;
    my $rm_fh = shift; #filehandle for rmdup
    my $mk_fh = shift; #filehandle for mkdup
    my $mate_is_deduped = shift;
    my $paired = shift;

    if (! $sam_bucket->[1]) {
	#print singletons to output
	$counter->{kept}++;
	$counter->{singleton}++;
	print $rm_fh "$sam_bucket->[0]\n" if ($output_type eq "all" || $output_type eq "rm");
	print $mk_fh "$sam_bucket->[0]\n" if ($output_type eq "all" || $output_type eq "mk");
	
	#for paired-end mode, interpolate skipped sam records
	print_skipped_sam_records($skip_bucket, $output_type, $counter, $rm_fh, $mk_fh, $mate_is_deduped) if ($skip_bucket->[0]);
	
    } else {
	#store index to read by mID and MAPQ
	my %data;
	for(my $i = 0; $i < @{$sam_bucket}; $i++) {
	    my %sam;
	    &parse_sam_record_w_mid(\%sam, $sam_bucket->[$i], $paired);
	    $data{$sam{key}}->{$sam{mapq}}->{$i} = undef;
	}

	#find best MAPQ per mID and store index in new hash
	my %best_indices;
      MID: foreach my $key (keys %data) {
	  foreach my $mapq (sort {$b <=> $a} keys %{$data{$key}}) {
	      foreach my $index (sort {$a <=> $b} keys %{$data{$key}->{$mapq}}) {
		  $best_indices{$index} = undef;
		  $counter->{best}++;
		  next MID;
	      }
	  }
      }

	#evaluate each read and check if index exists for deduping
	for(my $i = 0; $i < @{$sam_bucket}; $i++) {
	    my @fields = split(/\t/, $sam_bucket->[$i]);
	    if (exists $best_indices{$i}) {
		$counter->{kept}++;
		$fields[1] -= 1024 if ($fields[1] & 1024); #subtract from bitflag
		print $rm_fh join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "rm");
		print $mk_fh join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "mk");
	    } else {
		$counter->{dedup}++;
		$mate_is_deduped->{$fields[0]}++;
		$fields[1] += 1024 unless ($fields[1] & 1024); #add to bitflag
		print $mk_fh join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "mk");
	    }
	}
	
	#for paired-end mode, interpolate skipped sam records after each processed bucket
	&print_skipped_sam_records($skip_bucket, $output_type, $counter, $rm_fh, $mk_fh, $mate_is_deduped) if ($skip_bucket->[0]);

    }

    #empty bucket arrays for next iteration
    @{$sam_bucket} = ();
    @{$skip_bucket} = ();

}

sub print_skipped_sam_records {
    my $skip_bucket = shift;
    my $output_type = shift;
    my $counter = shift;
    my $rm_fh = shift; #filehandle for rmdup
    my $mk_fh = shift; #filehandle for mkdup
    my $mate_is_deduped = shift;
    
    #process each skipped read for paired-end mode
    for(my $i = 0; $i < @{$skip_bucket}; $i++) {
	my @fields = split(/\t/, $skip_bucket->[$i]);

	#dedup if mate was previously deduped also
	if (exists $mate_is_deduped->{$fields[0]}) {
	    $counter->{deduppair}++;
	    $fields[1] += 1024 unless ($fields[1] & 1024); #add to bitflag
	    print $mk_fh join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "mk");
	    delete $mate_is_deduped->{$fields[0]};
	    #otherwise keep it
	} else {
	    $counter->{kept}++;
	    $fields[1] -= 1024 if ($fields[1] & 1024); #subtract from bitflag
	    print $rm_fh join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "rm");
	    print $mk_fh join("\t", @fields) . "\n" if ($output_type eq "all" || $output_type eq "mk");
	}
    }
}

sub parse_sam_record_w_mid {
    my $data = shift;
    my $record = shift;
    my $paired = shift;
    
    #parse SAM record
    my ($qname, $bitflag, $rname, $pos, $mapq, $cigar, $rnext, $pnext, $tlen, $seq) = (split(/\t/, $record))[0..9];

    #extract MID from QNAME suffix
    my ($mID) = $qname =~ /_([^_]*)(?:\s)?/;

    #assign strand
    my $strand = "f";
    $strand = "r" if ($bitflag & 16);

    #create unique key to remove PCR duplicates
    my $key;
    if ($paired) {
	$key = "${rname}_${pos}_${mID}_${tlen}";
    } else {
	$key = "${rname}_${pos}_${mID}_${strand}";
    }

    #store data record
    $data->{sam} = $record;
    $data->{qname} = $qname;
    $data->{bitflag}= $bitflag;
    $data->{rname}= $rname;
    $data->{pos}= $pos;
    $data->{mapq}= $mapq;
#$data->{cigar}= $cigar;
    $data->{rnext}= $rnext;
    $data->{pnext}= $pnext;
    $data->{tlen}= $tlen;
#$data->{seq}= $seq;
    $data->{mid}= $mID;
    $data->{strand}= $strand;
    $data->{key}= $key;
}

