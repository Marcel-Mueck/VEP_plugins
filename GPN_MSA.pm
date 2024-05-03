=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

 Ensembl <http://www.ensembl.org/info/about/contact/index.html>
    
=cut

=head1 NAME

 GPN_MSA

=head1 SYNOPSIS

 mv GPN_MSA.pm ~/.vep/Plugins
 ./vep -i variations.vcf --plugin GPN_MSA,snv=/FULL_PATH_TO_SCORE_FILE/whole_genome_SNVs.tsv.gz
=head1 DESCRIPTION

 A VEP plugin that retrieves MSA scores for variants from one or more
 tabix-indexed MSA data files.

=cut

package GPN_MSA;

use strict;
use warnings;
use Data::Dumper;
use Bio::EnsEMBL::Variation::Utils::Sequence qw(get_matched_variant_alleles);

use Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepTabixPlugin);

# List here all columns in headers that should be included
my %SV_TERMS = (
    insertion => "INS",
    deletion => "DEL",
    duplication => "DUP",
  );

my %INCLUDE_COLUMNS = (
    "GRPN_MSA_score" => {
      "name" => "GRPN_MSA_score",
      "description" => 'GRPN_MSA score.'
    }
);



  my $ALT_NUM = 0;

sub new {
  my $class = shift;
  
  my $self = $class->SUPER::new(@_);

  # Test if tabix exists
  die "\nERROR: tabix does not seem to be in your path\n" unless `which tabix 2>&1` =~ /tabix$/;

  $self->expand_left(0);
  $self->expand_right(0);
  $self->get_user_params();

  my $params = $self->params_to_hash();
  my @files;
  $self->{non_sv_ann_file} = 0;
  # Check files in arguments
  if (!keys %$params) {
    @files = map { $_ ne "1" ? ($_) : () } @{$self->params};
    $self->{force_annotate} = $self->params->[-1] eq "1" ? 1 : 0;
  } else {
    my @param_keys = keys %{$params};

    for my $key ( @param_keys ){
      next if $key eq "force_annotate";
      push @files, $params->{$key};

      $self->{non_sv_ann_file} = 1 if ($key eq "snv" || $key eq "indels");
    }

    $self->{force_annotate} = $params->{force_annotate} ? 1 : 0;
  }

  die "\nERROR: No score files specified\nTip: Add a file after command, example:\nvep ... --plugin GPN_MSA,/FULL_PATH_TO_GPN_MSA_FILE/whole_genome_SNVs.tsv.gz\n" unless @files > 0;
  $self->add_file($_) for @files;

  warn "WARNING: Using snv and/or indels GPN_MSA annotation file with structural variant can increase run time exponentially. ".
        "Consider creating separate input files for SNV/indels and SV and use appropriate GPN_MSA annotation file.\n"
  if $self->{force_annotate};

  my $assembly = $self->{config}->{assembly};

  $self->{header} = ();


  foreach my $file (@files) {
    open IN, "tabix -f -h ".$file." 1:1-1 |";

    my @lines = <IN>;




    # while (my $line = shift @lines) {

    #   next if (rindex $line, "#Chrom", 0);
    #   chomp $line;
    #   $self->{$file} = $line;
    # }

    # # Make sure it has a known prefix in header
    # die "'#Chrom' was not found on header" unless $self->{$file};

    # my $file_check = 0;
    # # Conditional header
    # for (split /\t/, $self->{$file}){
    #   next unless (exists($INCLUDE_COLUMNS{$_}));
    #   $file_check = 1;
    #   $self->{header}{$INCLUDE_COLUMNS{$_}{"name"}} = $INCLUDE_COLUMNS{$_}{"description"} . " " ;
    # }

    # die "\nERROR: $file does not have a known column to be included" unless $file_check;

  }

  close IN;

  return $self;
}

sub variant_feature_types {
    return ['VariationFeature', 'StructuralVariationFeature'];
}

sub feature_types {
  return ['Feature','Intergenic'];
}

sub get_header_info {
  my $self = shift;
  return {
    GRPN_MSA_score => "Score from MSA model",
  }
}

sub run {
  my ($self, $tva) = @_;
  
  my $bvf = $tva->base_variation_feature;

  my ($start, $end, $allele, $ref, $so_term);

  # get allele
  if ($bvf->isa("Bio::EnsEMBL::Variation::VariationFeature")){
    $start = $bvf->{start};
    $end = $bvf->{end};
    $allele = $bvf->alt_alleles->[$ALT_NUM];
    $ref = $bvf->ref_allele_string;

    if (($ALT_NUM + 1) == scalar(@{$bvf->alt_alleles})) {
      $ALT_NUM = 0;
    } else {
      $ALT_NUM += 1;
    };

    return {} unless $allele =~ /^[ACGT-]+$/;

  } else {
    # Do not annotate sv if there is snv/indels annotation file
    return {} if ($self->{non_sv_ann_file} && !$self->{force_annotate});

    $start = $bvf->{start} - 1;
    $end = $bvf->{end};
    $so_term = $bvf->class_SO_term();
    $allele = $SV_TERMS{$so_term};
    $ref = "-";
  };

  my @data =  @{$self->get_data($bvf->{chr}, $start - 2, $end)};

  # Do not annotate if matched lines from annotation file is over threshold
  if(scalar @data > 100000 && !$self->{force_annotate}) {
    my $location = $bvf->{chr} . "_" . $start . "_" . $end;
    warn "WARNING: too many match found (", scalar @data, ") for GPN_MSA variant with location $location. No GPN_MSA annotation will be made. " .
         "Make sure you are not using SNVs/Indels GPN_MSA annotation file with structural variant as input. If you still want to annotate please use force_annotate=1."; 
    return {};
  }
  
  foreach (@data) {

    my $matches = get_matched_variant_alleles(
      {
        ref    => $ref,
        alts   => [$allele],
        pos    => $start,
        strand => $bvf->strand
      },
      {
       ref  => $_->{ref},
       alts => [$_->{alt}],
       pos  => $_->{start},
      }
    );
    #warn Dumper($_->{result})if (@$matches);
    return $_->{result} if (@$matches);
  }
  return {};
}

sub parse_data {
  my ($self, $line, $file) = @_;

  my @values = split /\t/, $line;

  

  my $c = $values[0];
  my $s = $values[1];
  my $ref = $values[2] || "-";
  my $alt = $values[3];
  my $result = $values[4];
  my $end = $values[1];
  # Conditional result
  my %result = ();
  # foreach (keys %INCLUDE_COLUMNS){
  #   next unless (exists($data{$_}));
  #   #$result{$INCLUDE_COLUMNS{$_}{"name"}} = $data{$_};
  #   $result{$INCLUDE_COLUMNS{$_}{"name"}} = $values[4];
  # }

  # do VCF-like coord adjustment for mismatched subs
  # my $end = ($s + length($ref)) - 1;
  # if(length($alt) != length($ref)) {
  #   my $first_ref = substr($ref, 0, 1);
  #   my $first_alt = substr($alt, 0, 1);
  #   if ($first_ref eq $first_alt) {
  #     $s++;
  #     $ref = substr($ref, 1);
  #     $alt = substr($alt, 1);
  #     $ref ||= '-';
  #     $alt ||= '-';
  #   }
  # }

  return {
    ref => $ref,
    alt => $alt,
    start => $s,
    end => $end,
    result => {GRPN_MSA_score => $result}
  };
}

sub get_start {
  return $_[1]->{start};
}

sub get_end {
  return $_[1]->{end};
}

1;
