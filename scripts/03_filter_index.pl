#!/usr/bin/env perl
# ==========================================================
# SoMoSeq preprocessing - step 3: filter index fastq to trimmed reads
# ==========================================================
# creates a filtered index fastq from the zUMIs index file based on valid read IDs from the trimmed R1 fastq file (output of step 2).
# usage: perl 03_filter_index.pl <RAW_DATA_DIR>
# submit via workflow manager with exact RAM and CPU specs; DO NOT run directly on the command line.
# <RAW_DATA_DIR> must contain reads_for_zUMIs_trimmed.re.R1.fastq.gz and
# reads_for_zUMIs.index.fastq.gz. Writes reads_for_zUMIs_filtered.index.fastq.gz, the index file used by step 4 (zUMIs run).

use strict;
use warnings;
use IO::Zlib;

my $dir = $ARGV[0] or die "Usage: $0 <RAW_DATA_DIR>\n";
my $trimmed_r1_file = "$dir/reads_for_zUMIs_trimmed.re.R1.fastq.gz";
my $index_file = "$dir/reads_for_zUMIs.index.fastq.gz";
my $output_file = "$dir/reads_for_zUMIs_filtered.index.fastq.gz";
my $valid_reads_file = "$dir/valid_reads.tmp";


print "Extracting valid read IDs...\n";
open(my $valid_handle, ">", $valid_reads_file) or die "Cannot write to $valid_reads_file: $!";
my $trimmed_r1_handle = IO::Zlib->new($trimmed_r1_file, "rb")
    or die "Cannot open $trimmed_r1_file: $!";

while (my $line = <$trimmed_r1_handle>) {
    if ($. % 4 == 1) {  
        chomp($line);
        $line =~ s/^@//; 
        print $valid_handle "$line\n";  
    }
    if ($. % 1000000 == 0) {
        print STDERR "Processed $. lines from trimmed R1 file...\n";
    }
}
close($trimmed_r1_handle);
close($valid_handle);
print "Valid read IDs written to $valid_reads_file.\n";

# Filterusing valid read IDs
print "Filtering index file...\n";
open(my $valid_reads_fh, "<", $valid_reads_file) or die "Cannot open $valid_reads_file: $!";
my %valid_reads;
while (my $line = <$valid_reads_fh>) {
    chomp($line);
    $valid_reads{$line} = 1;
}
close($valid_reads_fh);

my $index_handle = IO::Zlib->new($index_file, "rb")
    or die "Cannot open $index_file: $!";
my $output_handle = IO::Zlib->new($output_file, "wb")
    or die "Cannot create $output_file: $!";

while (my $line1 = <$index_handle>) {
    my $line2 = <$index_handle>;
    my $line3 = <$index_handle>;
    my $line4 = <$index_handle>;

    chomp($line1);
    my $header = $line1;
    $header =~ s/^@//; 

    if (exists $valid_reads{$header}) {
        print $output_handle "$line1\n$line2$line3$line4";
    }
}

close($index_handle);
close($output_handle);

print "Filtered index file written to $output_file.\n";
