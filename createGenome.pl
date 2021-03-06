#!/usr/bin/perl

=head1 NAME

createGenome.pl

=head1 DESCRIPTION

Create a new genome modifying a reference and adding variations from
Complete Genomics results, we apply SNVs (SNPs and short Indels) and
SV (deletions, inversions, translocations). The result can be a 
diploid (phasing can be applied) or haploid genome.

=head1 USAGE

createGenome.pl [PARAMS]

    Parameter       Description
    -r --reference  Fasta file with reference genome
    -s --snv        Table with SNVs
    -v --sv         Table with SVs
    -o --out        Output file
    -t --type       Use this types of variations, currently supported:
                    In masterVar: [snp, ins, del, sub]
                    In hcSV: [deletion, distal-duplication, inversion,
                    interchromosomal, probable-inversion, tandem-duplication]
                    *complex is always skipped, 'cause it's complex
    -x --sex        define sex [M, F]. Default is M
    -d --diploid    Create a diploid genome
    -h --help       Print this screen and exit
    -v --verbose    Verbose mode
    --version       Print version and exit

=head1 EXAMPLES

    1. Generate a haploid genome
    perl createGenome.pl -r hg19.fa -s masterVar -v highConfindenceSV -o new.fa
    
    2. Use SNPs and create a diplod genome (the second genome is 
    named new.fa_2)
    perl  createGenome.pl -r hg19.fa -s masterVar -t snp -d -o new.fa

=head1 AUTHOR

Juan Caballero, Institute for Systems Biology @ 2012

=head1 CONTACT

jcaballero@systemsbiology.org

=head1 LICENSE

This is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with code.  If not, see <http://www.gnu.org/licenses/>.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Default parameters
my $help      = undef;         # Print help
my $verbose   = undef;         # Verbose mode
my $version   = undef;         # Version call flag
my $reference = undef;
my $out       = undef;
my $snv       = undef;
my $sv        = undef;
my $type      = 'snp,ins,del,sub,deletion,distal-duplication,inversion,interchromosomal,probable-inversion,tandem-duplication';
my $diploid   = undef;
my $sex       =   'M';

# Main variables
my $our_version = 0.1;        # Script version number
my %hap1;
my %hap2;
my %orig;
my %event1;
my %event2;
my %types;

# Calling options
GetOptions(
    'h|help'           => \$help,
    'v|verbose'        => \$verbose,
    'version'          => \$version,
    'r|reference=s'    => \$reference,
    's|snv:s'          => \$snv,
    'v|sv:s'           => \$sv,
    'o|out=s'          => \$out,
    't|type:s'         => \$type,
    'x|sex:s'          => \$sex,
    'd|diploid'        => \$diploid
) or pod2usage(-verbose => 2);

printVersion() if (defined $version);
pod2usage(-verbose => 2) if  (defined $help);
pod2usage(-verbose => 2) if !(defined $reference);
pod2usage(-verbose => 2) if !(defined $out);
pod2usage(-verbose => 2) if !(defined $snv or defined $sv);

# Load sequences to memory
readFasta($reference, \%hap1);
if (defined $diploid) {
    foreach my $chr (keys %hap1) {
        next if ($chr eq 'chrY' and $sex eq 'F');
        next if ($chr eq 'chrX' and $sex eq 'M');
        next if ($chr eq 'chrM');
        $hap2{$chr} = $hap1{$chr};
    }
}

# Defining variation categories
my @types = split(/,/, $type);
foreach my $t (@types) {
    $types{$t} = 1;
}

# Doing changes for SNV and small variations that don't affect the coordinates
if (defined $snv) {
    warn "processing small variations in $snv\n" if (defined $verbose);
    my $snv_fh = defineFH($snv);
    open F, "$snv_fh" or die "cannot open file $snv\n";
    while (<F>) {
        next if  (m/^#/);
        next if  (m/^>/);
        next if  (m/^\n/);
        my ($nul1, $nul2, $chr, $ini, $end, $call, $type, $ref, $all1, $all2, $sco1, $sco2) = split (/\s+/, $_);
        next if  ($call eq 'no-call');
        next if  ($type eq 'ref');
        next if  ($type eq 'complex');
        next if !(defined $types{$type});
        next if !(defined $hap1{$chr});
        my $lenr = length $ref;
        my $len1 = length $all1;
        my $len2 = length $all2;
        if ($lenr == $len1) {
            substr($hap1{$chr}, $ini, $len1) = $all1 if ($all1 ne '?' and $ref ne $all1);
        }
        elsif ($type eq 'del') {
            substr($hap1{$chr}, $ini, $lenr) = 'X' x $lenr if ($len1 >= 1);
        }
        else {
            $event1{$chr}{$ini} = "$type:$ini:$end:$ref:$all1" if($ref ne $all1 and $all1 ne '?');
        }
        
        if (defined $diploid) {
            if ($lenr == $len2) {
                substr($hap2{$chr}, $ini, $len2) = $all2 if ($all2 ne '?' and $ref ne $all2);
            }
            elsif ($type eq 'del') {
                substr($hap2{$chr}, $ini, $lenr) = 'X' x $lenr if ($len2 >= 1);
            }
            else {
                $event2{$chr}{$ini} = "$type:$ini:$end:$ref:$all2" if($ref ne $all2 and $all2 ne '?');
            }
        }
        
    }
    close F;
}

# Doing changes for SV and other variations that don't affect the coordinates
if (defined $sv) {
    my ($s1, $s2, $f1, $f2);
    warn "processing structural variations in $sv\n" if (defined $verbose);
    my $sv_fh = defineFH($sv);
    open F, "$sv_fh" or die "cannot open file $sv\n";
    while (<F>) {
        next if  (m/^#/);
        next if  (m/^>/);
        next if  (m/^\n/);
        my ($nul1, $type, $nul2, $nul3, $frq, $ochr, $oini, $oend, $olen, $odir, $dchr, $dini, $dend, $dlen, $ddir) = split (/\s+/, $_);
        next if !(defined $types{$type});
        next if !(defined $hap1{$ochr});
        if ($type eq 'deletion') {
            if (defined $diploid) {
                if ($frq =~ m/;/) {
                    ($f1, $f2) = split (/;/, $frq);
                    substr($hap1{$ochr}, $oini, $olen) = 'X' x $olen if ($f1 > 0.01);
                    substr($hap2{$ochr}, $oini, $olen) = 'X' x $olen if ($f2 > 0.01);
                }
                else {
                    substr($hap1{$ochr}, $oini, $olen) = 'X' x $olen;
                    substr($hap2{$ochr}, $oini, $olen) = 'X' x $olen;
                }
            }
            else {
                substr($hap1{$ochr}, $oini, $olen) = 'X' x $olen;
            }
        }
        elsif ($type eq 'probable-inversion') {
            if (defined $diploid) {
                if ($frq =~ m/;/) {
                    ($f1, $f2) = split (/;/, $frq);
                    $s1 = substr($hap1{$ochr}, $oini, $olen);
                    $s2 = substr($hap2{$ochr}, $oini, $olen);
                    substr($hap1{$ochr}, $oini, $olen) = reverse $s1 if ($f1 > 0.01);
                    substr($hap2{$ochr}, $oini, $olen) = reverse $s2 if ($f2 > 0.01);
                }
                else {
                    $s1 = substr($hap1{$ochr}, $oini, $olen);
                    $s2 = substr($hap2{$ochr}, $oini, $olen);
                    substr($hap1{$ochr}, $oini, $olen) = reverse $s1;
                    substr($hap2{$ochr}, $oini, $olen) = reverse $s2;
                }
            }
            else {
                $s1 = substr($hap1{$ochr}, $oini, $olen);
                substr($hap1{$ochr}, $oini, $olen) = reverse $s1;
            }
        }
        else {
            if (defined $diploid) {
                if ($frq =~ m/;/) {
                    my ($f1, $f2) = split (/;/, $frq);
                    $event1{$dchr}{$dini} = "$type:$ochr:$oini:$olen:$odir:$dchr:$dini:$dlen:$ddir" if ($f1 > 0.01);
                    $event2{$dchr}{$dini} = "$type:$ochr:$oini:$olen:$odir:$dchr:$dini:$dlen:$ddir" if ($f2 > 0.01);
                }
                else {
                    $event1{$dchr}{$dini} = "$type:$ochr:$oini:$olen:$odir:$dchr:$dini:$dlen:$ddir";
                    $event2{$dchr}{$dini} = "$type:$ochr:$oini:$olen:$odir:$dchr:$dini:$dlen:$ddir";
                }
            }
            else {
                $event1{$dchr}{$dini} = "$type:$ochr:$oini:$olen:$odir:$dchr:$dini:$dlen:$ddir";
            }
        }
    }
    close F;
}

# Now we start altering the sequence with variations (ins & sub in masterVar;
# distal-duplication, interchromosomal, inversion & tandem-duplication in hcSV)
# that will modify the coordinate system.
%orig = %hap1;
foreach my $chr (keys %event1) {
    my @pos = sort {$b<=>$a} (keys %{ $event1{$chr} }); # yes, we do changes in reverse order
    foreach my $pos (@pos) {
        my $event = $event1{$chr}{$pos};
        if    ($event =~ m/^ins:/) { doInsSub   ($event, \$hap1{$chr}); }
        elsif ($event =~ m/^sub:/) { doInsSub   ($event, \$hap1{$chr}); }
        elsif ($event =~ m/^dist/) { doTransLoca($event, \$hap1{$chr}); }
        elsif ($event =~ m/^inte/) { doTransLoca($event, \$hap1{$chr}); }
        elsif ($event =~ m/^inve/) { doTransLoca($event, \$hap1{$chr}); }
        elsif ($event =~ m/^tand/) { doTandemDup($event, \$hap1{$chr}); }
        else {
            # do nothig, we don't know what is it
        }
    }    
}

if (defined $diploid) {
    %orig = %hap2;
    foreach my $chr (keys %event2) {
        my @pos = sort {$b<=>$a} (keys %{ $event2{$chr} }); # yes, we do changes in reverse order
        foreach my $pos (@pos) {
            my $event = $event1{$chr}{$pos};
            if    ($event =~ m/^ins:/) { doInsSub   ($event, \$hap2{$chr}); }
            elsif ($event =~ m/^sub:/) { doInsSub   ($event, \$hap2{$chr}); }
            elsif ($event =~ m/^dist/) { doTransLoca($event, \$hap2{$chr}); }
            elsif ($event =~ m/^inte/) { doTransLoca($event, \$hap2{$chr}); }
            elsif ($event =~ m/^inve/) { doTransLoca($event, \$hap2{$chr}); }
            elsif ($event =~ m/^tand/) { doTandemDup($event, \$hap2{$chr}); }
            else {
                # do nothig, we don't know what is it
            }
        }
    }
}

# Writing final sequence(s)
writeFasta($out,      \%hap1);
writeFasta("$out\_2", \%hap2) if (defined $diploid);

###################################
####   S U B R O U T I N E S   ####
###################################

sub printVersion {
    print "$0 $our_version\n";
    exit 1;
}

sub readFasta {
    my ($fi, $seq_ref) = @_;
    my $chr;
    warn "reading sequences from $fi\n" if (defined $verbose);
    my $fh = defineFH($fi);
    open F, "$fh" or die "cannot open $fi\n";
    while (<F>) {
        chomp;
        if (m/>(.+)/) {
            $chr = $1;
            warn "  .. $chr\n" if (defined $verbose);
        }
        else {
            $$seq_ref{$chr} = $_;
        }
    }
    close F;
}

sub writeFasta {
    my ($fo, $seq_ref) = @_;
    warn "writing new sequences in $fo\n" if (defined $verbose);
    open O, ">$fo" or die "cannot write $fo\n";
    while (my ($chr, $seq) = each %$seq_ref) {
        print O ">$chr\n";
        $seq =~ tr/X//;
        while ($seq) {
            print O substr($seq, 0, 70), "\n";
            substr($seq, 0, 70) = '';
        }
    }
    close O;
}

sub defineFH {
    my ($fi) = @_;
    my $fh   = $fi;
    $fh = "gunzip -c $fi | " if ($fi =~ m/gz$/);
    $fh = "bgrep unzip2 -c $fi | " if ($fi =~ m/bz2$/);
    my $res = undef;
    return $res;
}

sub revcomp {
    my ($s) = @_;
    my $r   = reverse $s;
    $r =~ tr/ACGTacgt/TGCAtgca/;
    return $r;
}

sub doInsSub {
    my ($evt, $seq_ref) = @_;
    my ($type, $ini, $end, $ref, $ins) = split (/:/, $evt);
    if ($ini == $end) {
        my $b = substr($$seq_ref, $ini, 1);
        substr($$seq_ref, $ini, 1) = "$b$ins";
    }
    else {
        my $len = $end - $ini;
        my $b   = substr($$seq_ref, $ini, $len);
        substr($$seq_ref, $ini, $len) = "$b$ins";
    }
}

sub doTransLoca {
    my ($evt, $seq_ref) = @_;
    my ($type,$ochr,$oini,$olen,$odir,$dchr,$dini,$dlen,$ddir) = split (/:/, $evt);
    my $orig = substr($orig{$ochr}, $oini, $olen);
    $orig = revcomp($orig) if ($odir ne $ddir);
    substr($$seq_ref, $dini, $dlen) = $orig;
}

sub doTandemDup {
    my ($evt, $seq_ref) = @_;
    my ($type,$ochr,$oini,$olen,$odir) = split (/:/, $evt);
    my $orig = substr($orig{$ochr}, $oini, $olen);
    substr($$seq_ref, $oini, $olen) = $orig x 2;
}
