#!/usr/bin/env perl
#
#  FASTQ Splitter  -  a script for partitioning a FASTQ file into pieces
#
#  Version 0.1.2 (April 24, 2014)
#
#  Copyright (c) 2014 Kirill Kryukov
#
#  This software is provided 'as-is', without any express or implied
#  warranty. In no event will the authors be held liable for any damages
#  arising from the use of this software.
#
#  Permission is granted to anyone to use this software for any purpose,
#  including commercial applications, and to alter it and redistribute it
#  freely, subject to the following restrictions:
#
#  1. The origin of this software must not be misrepresented; you must not
#     claim that you wrote the original software. If you use this software
#     in a product, an acknowledgment in the product documentation would be
#     appreciated but is not required.
#  2. Altered source versions must be plainly marked as such, and must not be
#     misrepresented as being the original software.
#  3. This notice may not be removed or altered from any source distribution.
#

use Getopt::Long;
use strict;

$| = 1;

my $start_time = time;

my ($opt_n_parts,$opt_part_size,$opt_measure,$opt_eol,$opt_check,$opt_version,$opt_help);
GetOptions("n-parts=i"   => \$opt_n_parts,
           "part-size=i" => \$opt_part_size,
           "measure=s"   => \$opt_measure,
           "eol=s"       => \$opt_eol,
           "check"       => \$opt_check,
           "version"     => \$opt_version,
           "help"        => \$opt_help)
or die "Can't parse command line arguments\n";

sub show_version
{
    print q{FASTQ Splitter 0.1.2
Copyright (c) 2014 Kirill Kryukov
};
}

sub show_help()
{
    print q{Usage: fastq-splitter.pl [options] <file>...
Options:
    --n-parts <N>        - Divide into <N> parts
    --part-size <N>      - Divide into parts of size <N>
    --measure (all|seq|count) - Specify whether all data, sequence length, or
                           number of sequences is used for determining part
                           sizes ('all' by default).
    --eol (dos|mac|unix) - Choose end-of-line character ('unix' by default).
    --check              - Check FASTQ correctness.
    --version            - Show version.
    --help               - Show help.
};
}

if ($opt_version) { show_version(); }
if ($opt_help) { show_help(); }
if ($opt_help or $opt_version) { exit; }

if (!defined($opt_n_parts) and !defined($opt_part_size))
{
    if (!$opt_help and !$opt_version) { show_version(); show_help(); }
    print STDERR "Error: Splitting method is not specified ('--n-parts' or '--part-size' option)\n";
    exit;
}

if (!@ARGV) { die "File for splitting is not specified\n"; }

if (defined($opt_n_parts) and $opt_n_parts <= 0) { die "Non-positive number of parts\n"; }
if (defined($opt_part_size) and $opt_part_size <= 0) { die "Non-positive part size\n"; }
if (defined($opt_measure) and $opt_measure ne 'all' and $opt_measure ne 'seq' and $opt_measure ne 'count') { die "Unknown value of --measure option\n"; }
if (defined($opt_eol) and $opt_eol ne 'dos' and $opt_eol ne 'mac' and $opt_eol ne 'unix') { die "Unknown value of --eol option\n"; }

my $n_parts = defined($opt_n_parts) ? $opt_n_parts : 0;
my $part_size = defined($opt_part_size) ? $opt_part_size : 0;
my $eol = defined($opt_eol) ? (($opt_eol eq 'dos') ? "\x0D\x0A" : ($opt_eol eq 'mac') ? "\x0D" : "\x0A") : "\x0A";
my $eol_len = length($eol);
my $measure = defined($opt_measure) ? (($opt_measure eq 'count') ? 0 : ($opt_measure eq 'seq') ? 1 : 2) : 2;
my ($total_nseq,$total_seq_length,$line_count);

foreach my $infile (@ARGV) { split_file($infile); }
 
my $end_time = time;
my $elapsed_time = $end_time - $start_time;
print "All done, $elapsed_time second", (($elapsed_time==1)?'':'s'), " elapsed\n";

sub split_file
{
    my ($infile) = @_;
    if (!-e $infile or !-f $infile) { print "Can't find file \"$infile\"\n"; return; }

    my $base = $infile;
    my $ext = '';
    if ($base =~ /^(.+?)(\.(fastq))$/i) { ($base,$ext) = ($1,$2); }

    print "$infile";
    my ($total_measured_size,$n_parts_measured) = (0,0);
    if ($n_parts and $part_size) { $n_parts_measured = $n_parts; }
    else
    {
        ($total_nseq,$total_seq_length,$line_count) = (0,0,0);
        ($total_measured_size,$n_parts_measured) = get_file_size($infile);
        print ": $total_nseq sequences, $total_seq_length bp";
    }

    ($total_nseq,$total_seq_length,$line_count) = (0,0,0);

    if ($part_size)
    {
        print ' => ', ($n_parts ? 'extracting' : 'dividing into'), ' ', $n_parts_measured, ' part', ($n_parts_measured > 1 ? 's' : ''),
              " of <= $part_size ", ($measure ? (($measure > 1) ? 'bytes' : 'bp') : 'sequences'), "\n";
        open(my $IN,'<',$infile) or die "Error: Can't open file \"$infile\"\n";
        my $part = 1;
        my $this_part_size = 0;
        my $num_len = length($n_parts_measured);
        my $out_file = sprintf("%s.part-%0*d%s",$base,$num_len,$part,$ext);
        open(my $OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
        binmode $OUT;
        while (!eof $IN)
        {
            my ($entry_size,$entry) = get_fastq_entry($IN,$opt_check);
            if ($this_part_size and ($this_part_size + $entry_size > $part_size))
            {
                close $OUT;
                $part++;
                if (defined($opt_n_parts) and ($part > $n_parts)) { last; }
                $out_file = sprintf("%s.part-%0*d%s",$base,$num_len,$part,$ext);
                open ($OUT,'>',$out_file) or die "Can't create output file \"$out_file\"\n";
                binmode $OUT;
                $this_part_size = 0;
            }
            print $OUT $entry;
            $this_part_size += $entry_size;
        }
        close $IN;
        close $OUT;
    }
    else
    {
        print " => dividing into $n_parts part", ($n_parts > 1 ? 's' : ''), "\n";
        open(my $IN,'<',$infile) or die "Error: Can't open file \"$infile\"\n";
        my $written = 0;
        my $num_len = length($n_parts);
        for (my $part = 1; $part <= $n_parts; $part++)
        {
            my $should_write = sprintf("%.0f",$total_measured_size * $part / $n_parts);
            my $part_file = $base;
            if ($part_file !~ /\.part-\d+$/) { $part_file .= '.part'; }
            $part_file .= sprintf("-%0*d%s",$num_len,$part,$ext);
            open(my $OUT,'>',$part_file) or die "Error: Can't create file \"$part_file\"\n";
            binmode $OUT;
            while ($written < $should_write)
            {
                my ($entry_size,$entry) = get_fastq_entry($IN,$opt_check);
                print $OUT $entry;
                $written += $entry_size;
            }
            close $OUT;
        }
        close $IN;
    }
}

sub get_file_size
{
    my ($file) = @_;
    open(my $IN,'<',$file) or die "Error: Can't open file \"$file\"\n";
    my ($measured_size,$n_parts_measured,$this_part_size) = (0,1,0);
    while (!eof $IN)
    {
        my ($entry_size) = get_fastq_entry($IN,0);
        $measured_size += $entry_size;
        if ($this_part_size and ($this_part_size + $entry_size > $part_size)) { $this_part_size = $entry_size; $n_parts_measured++; }
        else { $this_part_size += $entry_size; }
    }
    close $IN;
    return ($measured_size,$n_parts_measured);
}

sub get_fastq_entry
{
    my ($IN,$check) = @_;
    my $truncated = "Error: Incomplete FASTQ entry at the end -> looks like truncated input!\n";

    if (eof $IN) { return (0,''); }
    my $name_line = <$IN>;
    $line_count++;
    if (substr($name_line,0,1) ne '@') { die "Error parsing line $line_count: FASTQ entry does not start with '\@':\n$name_line"; }
    my $name = substr($name_line,1);
    $name =~ s/[\x0D\x0A]+$//;
    my $nlen = length($name);

    if (eof $IN) { print STDERR $truncated; return (0,''); }
    my $seq = <$IN>;
    $line_count++;
    $seq =~ s/[\x0D\x0A]+$//;
    my $slen = length($seq);

    if (eof $IN) { print STDERR $truncated; return (0,''); }
    my $plus_line = <$IN>;
    $line_count++;
    if (substr($plus_line,0,1) ne '+') { die "Error parsing line $line_count: Expecting '+', found '$plus_line'"; }
    my $name2 = substr($plus_line,1);
    $name2 =~ s/[\x0D\x0A]+$//;
    my $n2len = length($name2);
    if ($check and $name2 ne '' and $name2 ne $name) { print STDERR "Warning: Misformatted FASTQ entry in input line $line_count: Sequence identifiers don't match:\n\@$name\n+$name2\n"; }

    if (eof $IN) { print STDERR $truncated; return (0,''); }
    my $qual = <$IN>;
    $line_count++;
    $qual =~ s/[\x0D\x0A]+$//;
    my $qlen = length($qual);
    if ($qlen != $slen)
    {
        if (eof($IN)) { print STDERR $truncated; return (0,''); }
        elsif ($check) { print STDERR "Warning: Misformatted FASTQ entry in input line $line_count: quality length ($qlen) differs from sequence length ($slen):\n$seq\n$qual\n"; }
    }

    my $entry_size = $measure ? (($measure > 1) ? ($slen + $qlen + $nlen + $n2len + $eol_len*4 + 2) : $slen) : 1;
    my $entry = '@' . $name . $eol . $seq . $eol . '+' . $name2 . $eol . $qual . $eol;

    $total_nseq++;
    $total_seq_length += $slen;
    return ($entry_size,$entry);
}
