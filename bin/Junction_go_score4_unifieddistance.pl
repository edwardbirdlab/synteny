#!/usr/bin/perl
use strict;
use warnings;

# Build Go file mapping
my @gos = `ls Go/*`;
my %go_key;
foreach my $sp (@gos) {
    chomp $sp;
    my @split = split(/\./, $sp);
    my @sp_folder = split("/", $split[0]);
    $go_key{$sp_folder[1]} = $sp;
}

# Check arguments
if (@ARGV != 2) {
    die "Usage: $0 <distance_in_genes> <species_name>\n";
}
my $distance = $ARGV[0];
my $species = $ARGV[1];

# Data structures
my %scores;         # For minimum scores
my %all_scores;     # For all scores (to average later)

# Read all gene score files
my @files = glob("*.gene_scores.txt");
foreach my $file (@files) {
    open my $fh, '<', $file or die "Could not open file '$file': $!";
    print "Combining $file\n";
    while (<$fh>) {
        chomp;
        my ($col1, $col2, $score) = split /\t/;
        my $key = "$col1\t$col2";
        next if $score eq "NA";

        # Track minimum
        if (exists $scores{$key}) {
            $scores{$key} = $score if $scores{$key} eq "NA" || $scores{$key} >= $score;
        } else {
            $scores{$key} = $score;
        }

        # Track all scores for average
        push @{ $all_scores{$key} }, $score;
    }
    close $fh;
}

# Write Summary.tsv
open my $out_fh, '>', "Summary.tsv" or die "Could not open file 'Summary.tsv': $!";
foreach my $key (keys %scores) {
    print $out_fh "$key\t$scores{$key}\n";
}
close $out_fh;

# Prepare background and closest output files
my $back = "$species.$distance.bg.txt";
open my $backsave, '>', $back or die "Could not open file '$back': $!";

my $closest_output_file = "$species.${distance}.closest.txt";
open my $closest_fh, '>', $closest_output_file or die "Could not open file '$closest_output_file': $!";

# Write closest scores
my $distance_count = 0;
my @closest_cutoff;
foreach my $key (keys %scores) {
    next if $scores{$key} eq "NA";
    if ($scores{$key} <= $distance) {
        $distance_count++;
        push @closest_cutoff, $key;
        my @splh = split("\t", $key);
        my $gene = $splh[1];
        $gene =~ s/^.*://;       # Remove scaffold: prefix if present
        $gene =~ s/^rna-//;      # Remove rna- prefix
        $gene =~ s/-/_/g;        # Replace hyphens with underscores
        print $closest_fh "$gene\n";
    }
    # Save to background
    my @splh = split("\t", $key);
    my $gene = $splh[1];
    $gene =~ s/^.*://;
    $gene =~ s/^rna-//;
    $gene =~ s/-/_/g;
    print $backsave "$gene\n";
}
close $closest_fh;
close $backsave;

print "Lowest scores written to 'Summary.tsv'.\n$distance_count genes were below the threshold of $distance\n";

# Build farthest list based on average score
my @avg_data;
foreach my $key (keys %all_scores) {
    my @values = @{ $all_scores{$key} };
    next unless @values;
    my $sum = 0;
    $sum += $_ for @values;
    my $avg = $sum / @values;

    my @fields = split(/\t/, $key);
    my $gene = $fields[1];
    $gene =~ s/^.*://;
    $gene =~ s/^rna-//;
    $gene =~ s/-/_/g;

    push @avg_data, [$gene, $avg];
}

# Sort by descending average score and get top N
@avg_data = sort { $b->[1] <=> $a->[1] } @avg_data;
my @last_values = @avg_data[0 .. ($distance_count - 1)];

# Write farthest list
my $farthest_output_file = "$species.${distance}.farthest.txt";
open my $farthest_fh, '>', $farthest_output_file or die "Could not open file '$farthest_output_file': $!";
foreach my $entry (@last_values) {
    my ($gene, $value) = @$entry;
    print $farthest_fh "$gene\n";
}
close $farthest_fh;

print "Results written to $farthest_output_file and $closest_output_file\n";

# Run GO enrichment analysis
print "ChopGO_VTS2_v12.pl -i $farthest_output_file --GO_file $go_key{$species} -bg $back\n";
print "ChopGO_VTS2_v12.pl -i $closest_output_file --GO_file $go_key{$species} -bg $back\n";
`ChopGO_VTS2_v12.pl -i $farthest_output_file --GO_file $go_key{$species} -bg $back`; 
`ChopGO_VTS2_v12.pl -i $closest_output_file --GO_file $go_key{$species} -bg $back`;

