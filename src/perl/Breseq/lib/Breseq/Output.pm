###
# Pod Documentation
###

=head1 NAME

Breseq::Output.pm

=head1 SYNOPSIS

Various functions for output.

=head1 DESCRIPTION

=head1 AUTHOR

Jeffrey Barrick
<jbarrick@msu.edu>

=head1 COPYRIGHT

Copyright 2009.  All rights reserved.

=cut

###
# End Pod Documentation
###

package Breseq::Output;
use strict;
use CGI qw/:standard *table *Tr *th *td *div/;
use Data::Dumper;
use File::Copy;

require Exporter;
our @ISA = qw( Exporter );
our @EXPORT = qw();

use Breseq::AlignmentOutput;
use Breseq::GenomeDiff;

## these style definitions are included between the
## HEAD tags of every generated page
our $header_style_string = <<ENDOFSTYLE;
body {
	font-family: sans-serif;
	font-size: 11pt;
}
th 	{
	background-color: rgb(0,0,0); 
	color: rgb(255,255,255);
}
table	{
	background-color: rgb(0,0,0); 
	color: rgb(0,0,0);
}

tr	{
	background-color: rgb(255,255,255); 
}

.mutation_in_codon	{
	color:red;
	text-decoration : underline;
}

.mutation_header_row {
	background-color: rgb(0,130,0);
}

.read_alignment_header_row {
	background-color: rgb(255,0,0);
}

.missing_coverage_header_row {
	background-color: rgb(0,100,100);
}

.new_junction_header_row {
	background-color: rgb(0,0,155);
}

.mutation_table_row_0	{
	background-color: rgb(255,255,255);
}
.mutation_table_row_1	{
	background-color: rgb(245,245,245);
}	

.polymorphism_table_row	{
	background-color: rgb(160,255,160);
}

.highlight_table_row	{
	background-color: rgb(192,255,255);
}

.reject_table_row	{
	background-color: rgb(255,200,165);
}

.information_table_row	{
	background-color: rgb(200,255,255);
}

.junction_repeat	{
	background-color: rgb(255,165,0); /* orange */
}
.junction_gene	{
}

ENDOFSTYLE


sub html_index
{
	my ($file_name, $settings, $summary, $ref_seq_info, $gd) = @_;

	open HTML, ">$file_name" or die "Could not open file: $file_name";    

    print HTML start_html(
			-title => "BRESEQ :: Index" . ($settings->{run_name} ne 'unnamed' ? " :: $settings->{run_name}" : ''), 
			-head  => style({type => 'text/css'}, $header_style_string),
	);
	
	## Write fastq read file information
	print HTML breseq_header_string($settings) . p;
    print HTML start_table({-border => 0, -cellspacing => 1, -cellpadding => 5});
	print HTML Tr(th(), th("fastq read file"), th("reads"), th("bases"), th("longest"));
	my $total_bases = 0;
	my $total_reads = 0;
	foreach my $read_file ($settings->read_files)
	{
		my $c = $summary->{sequence_conversion}->{reads}->{$read_file};
		print HTML Tr(
			td(a({-href=>$settings->html_path('error_rates_plot_file_name', {'#'=>$read_file})}, "errors")),
			td($read_file), 
			td({-align=>"right"}, commify($c->{num_reads})), 
			td({-align=>"right"},commify($c->{total_bases})), 
			td($c->{max_read_length} . "&nbsp;bases"), 
		);
		$total_bases += $c->{total_bases};
		$total_reads += $c->{num_reads};
	}
	print HTML Tr({-class=>"highlight_table_row"}, 
		td(),
		td(b("total")), 
		td({-align=>"right"},b(commify($total_reads))), 
		td({-align=>"right"},b(commify($total_bases))), 
		td(b($summary->{sequence_conversion}->{max_read_length} . "&nbsp;bases")), 
	);
	print HTML end_table();

	## Write reference sequence information
	print HTML p;
	print HTML start_table({-border => 0, -cellspacing => 1, -cellpadding => 5});
	print HTML Tr(
		th(),
		th(),
		th("reference sequence"), 
		th("length"), 
		th("description")
	);
	my $total_length = 0;
	foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
	{
		my $c = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id};
		print HTML Tr(
			td(a({-href=>$settings->html_path('coverage_plot_file_name', {'@'=>$seq_id})}, "coverage")), 
			td(a({-href=>$settings->html_path('unique_only_coverage_plot_file_name', {'@'=>$seq_id})}, "distribution")), 
			td($seq_id), 
			td({-align=>"right"},commify($c->{length})), 
			td($c->{definition})
		);
		$total_length+= $c->{length};
	}	
	
	print HTML Tr({-class=>"highlight_table_row"},
		td(),
		td(),
		td(b("total")), 
		td(b({-align=>"right"},commify($total_length))), 
		td()
	);
	
	## junction only reference sequences
	foreach my $seq_id (@{$ref_seq_info->{junction_only_seq_ids}})
	{
		my $c = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id};
		print HTML Tr(
			td({-colspan=>"2", -align=>"center"}, "junction&nbsp;only"), 
			td($seq_id), 
			td({-align=>"right"},commify($c->{length})), 
			td($c->{definition})
		);
		$total_length+= $c->{length};
	}
	
	print HTML end_table();		

	###
	## Mutation predictions
	###
	
	my @muts = $gd->list('SNP', 'INS', 'DEL', 'SUB', 'MOB');
	my $relative_path = $settings->file_name('local_evidence_path');
	$relative_path .= "/" if ($relative_path);
	my $one_ref_seq = scalar(keys %{$ref_seq_info->{ref_strings}}) == 1;
	print HTML p . html_mutation_table_string($gd, \@muts, $relative_path, undef, $one_ref_seq );

	###
	## Unassigned evidence
	###
	
	my @mc = $gd->filter_used_as_evidence($gd->list('MC'));
	if (scalar @mc > 0)
	{
		print HTML p . html_missing_coverage_table_string(\@mc, $relative_path, "Unassigned missing coverage evidence...");
	}
	
	my @jc = $gd->filter_used_as_evidence($gd->list('JC'));
	@jc = grep { !$_->{no_show} } @jc;	

	my @jcu = grep { !$_->{reject} } @jc;	
	if (scalar @jcu > 0)
	{
		print HTML p . html_new_junction_table_string(\@jcu, $relative_path, "Unassigned new junction evidence...");	
	}
	
	###
	## Rejected evidence
	###
	
	my @ra = $gd->filter_used_as_evidence($gd->list('RA'));	
	## don't print ones that overlap predicted deletions or were marked to not show
	@ra = grep { !$_->{deleted} && !$_->{no_show} } @ra;
	
	if (scalar @ra > 0)
	{
		print HTML p . html_read_alignment_table_string(\@ra, $relative_path, "Rejected read alignment evidence...");
	}
	
	my @jcr = grep { $_->{reject} } @jc;
	if (scalar @jcr > 0)
	{	
		print HTML p . html_new_junction_table_string(\@jcr, $relative_path, "Rejected new junction evidence...");	
	}
	
	print HTML end_html;
	close HTML;
}

sub html_header
{
	my ($title) = @_;
	return start_html(
			-title => $title, 
			-head  => style({type => 'text/css'}, $header_style_string),
	);
}

sub html_footer
{
	return end_html;
}

sub html_compare
{
	my ($file_name, $title, $gd, $one_ref_seq, $gd_name_list_ref) = @_;

	open HTML, ">$file_name" or die "Could not open file: $file_name";    

    print HTML start_html(
			-title => $title, 
			-head  => style({type => 'text/css'}, $header_style_string),
	);
	
	my @muts = $gd->mutation_list;
	
	print HTML html_mutation_table_string($gd, \@muts, undef, undef, $one_ref_seq, $gd_name_list_ref);

	print HTML end_html;
	close HTML;
}

sub html_compare_polymorphisms
{
	my ($file_name, $title, $list_ref) = @_;

	open HTML, ">$file_name" or die "Could not open file: $file_name";    

    print HTML start_html(
			-title => $title, 
			-head  => style({type => 'text/css'}, $header_style_string),
	);
		
	print HTML html_read_alignment_table_string($list_ref, undef, undef, 1);

	print HTML end_html;
	close HTML;
}

## Note that we do not overwrite past summary tables
## Instead, we move them. This is in case the script was
## Run multiple times.

sub html_summary_table
{
	my ($file_name, $settings, $summary, $times) = @_;
	
	## Create the current file...
	open HTML, ">$file_name" or die "Could not open file: $file_name";    
    
    print HTML start_html(
			-title => "BRESEQ :: Summary" . ($settings->{run_name} ? " :: $settings->{run_name}" : ''), 
			-head  => style({type => 'text/css'}, $header_style_string),
	);
    
#   print HTML h1("Program Options");
#	print HTML start_table({-width => "100%", -border => 1, -cellspacing => 0, -cellpadding => 3});
#	print HTML Tr(th("Option"), th("Setting"));
#	print HTML Tr(td("Fastq quality score style"), td($settings->{quality_type}));
#	print HTML end_table();
		
	## Write times
	print HTML p . h1("Execution Times");
	print HTML start_table({-width => "100%", -border => 1, -cellspacing => 0, -cellpadding => 3});
	print HTML Tr(th("Step"), th("Time"), th("Elapsed"));
	foreach my $time (@$times)
	{
		print HTML Tr(td($time->{_name}), td($time->{_formatted_time}), td($time->{_formatted_time_elapsed}));
	}
	my $formatted_total_time = time2string($times->[-1]->{_time} - $times->[0]->{_time});
	print HTML Tr({-class=>"highlight_table_row"}, td({-colspan=>2}, b("Total Execution Time")), td(b($formatted_total_time)));
	print HTML end_table();

	print HTML h1("Summary Statistics");

    print HTML h2("Read Sequence Files");
    print HTML start_table({-width => "100%", -border => 1, -cellspacing => 0, -cellpadding => 3});
	print HTML Tr(th("File"), th("Reads"), th("Bases"), th("Max Read Length"));
	my $total_bases = 0;
	my $total_reads = 0;
	foreach my $key (sort keys %{$summary->{sequence_conversion}->{reads}})
	{
		my $c = $summary->{sequence_conversion}->{reads}->{$key};
		print HTML Tr(td($key), td($c->{num_reads}), td($c->{total_bases}), td($c->{max_read_length}));
		$total_bases += $c->{total_bases};
		$total_reads += $c->{num_reads};
	}
	print HTML Tr({-class=>"highlight_table_row"}, td(b("Total")), td(b($total_reads)), td(b($total_bases)), td(b($summary->{sequence_conversion}->{max_read_length})));
	print HTML end_table();

	## Write reference sequence information
    print HTML h2("Reference Sequence Files");
	print HTML start_table({-width => "100%", -border => 1, -cellspacing => 0, -cellpadding => 3});
	print HTML Tr(th("Accession"), th("Length"), th("Description"));
	my $total_length = 0;
	foreach my $seq_id (sort keys %{$summary->{sequence_conversion}->{reference_sequences}})
	{
		my $c = $summary->{sequence_conversion}->{reference_sequences}->{$seq_id};
		print HTML Tr(td($seq_id), td($c->{length}), td($c->{definition}));
		$total_length+= $c->{length};
	}	
	print HTML Tr({-class=>"highlight_table_row"}, td(b("Total")), td(b($total_length)), td());
	print HTML end_table();

	close HTML;
}

sub breseq_header_string
{
	my ($settings) = @_;
	my $output_string = '';
	
	#copy over the breseq_graphic
	my $breseq_graphic_from_file_name = $settings->file_name('breseq_small_graphic_from_file_name');
	my $breseq_graphic_to_file_name = $settings->file_name('breseq_small_graphic_to_file_name');
	
	if (!-e $breseq_graphic_to_file_name)
	{
		copy($breseq_graphic_from_file_name, $breseq_graphic_to_file_name);
	}
	
	$output_string .= start_table({-width => "100%", -border => 0, -cellspacing => 0, -cellpadding => 3});
	$output_string .= start_Tr;
	$output_string .=  td(img({-src=>$settings->html_path('breseq_small_graphic_to_file_name')}));
	$output_string .= start_td({-width => "100%"});
	$output_string .= $settings->{byline};	
	$output_string .= " | " . a({-href=>$settings->html_path('log_file_name')}, 'command line log');
	$output_string .= end_td . end_Tr . end_table;
	
	return $output_string;
}


sub html_genome_diff_item_table_string
{
	my ($gd, $list_ref) = @_;
	
	return '' if (!defined $list_ref) || (scalar @$list_ref) == 0;
	my $first_item = $list_ref->[0];	
				
	##mutation
	if (length($first_item->{type}) == 3)
	{	
		return html_mutation_table_string($gd, $list_ref);
	}
	
	##evidence
	else
	{
		if ($first_item->{type} eq 'MC')
		{
			return html_missing_coverage_table_string($list_ref, undef, undef, 1);
		}
		elsif ($first_item->{type} eq 'RA')
		{
			return html_read_alignment_table_string($list_ref, undef, undef, 1);
		}
		elsif ($first_item->{type} eq 'JC')
		{
			return html_new_junction_table_string($list_ref, undef, undef, 1);
		}
	}
}


sub html_mutation_table_string
{	
		my ($gd, $list_ref, $relative_link, $legend_row, $one_ref_seq, $gd_name_list_ref) = @_;
		$relative_link = '' if (!defined $relative_link);
		my $output_str = '';

		my $q = new CGI;
		$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3});

		#####################
		#### HEADER LINE ####
		#####################

		my $header_text = "Predicted mutation";
		$header_text .= "s" if (scalar @$list_ref > 1);

		
		my @freq_header_list = ("freq");
		@freq_header_list = @$gd_name_list_ref if (defined $gd_name_list_ref);
		
		my $total_cols = 5 + scalar @freq_header_list;
		$total_cols += 1 if (!$one_ref_seq);
		$total_cols += 1 if (!defined $gd_name_list_ref); ## evidence column 
		$output_str.= Tr(th({-colspan => $total_cols, -align => "left", -class=>"mutation_header_row"}, $header_text));
		
		$output_str.= start_Tr();
		$output_str.= th("evidence") if (!defined $gd_name_list_ref); 
		$output_str.= th(nonbreaking("seq id")) if (!$one_ref_seq); 
		
		$output_str.= th(
			[
				"position",
				"mutation",
				@freq_header_list, 
				"annotation", 
				"gene", 
			]
		);	
		$output_str.= th({-width => "100%"}, "description"); 
		$output_str.= end_Tr;

		####################
		#### ITEM LINES ####
		####################

		foreach my $mut (@$list_ref)
		{	
			
			my $evidence_string = '';
			if (!defined $gd_name_list_ref) 
			{
				my $already_added_RA;
				EVIDENCE: foreach my $evidence_item ($gd->mutation_evidence_list($mut))
				{					
					if ($evidence_item->{type} eq 'RA')
					{
						next EVIDENCE if ($already_added_RA);
						$already_added_RA = 1;
					}
					$evidence_string .= "&nbsp;" if ($evidence_string);
					$evidence_string .= a({href => "$relative_link$evidence_item->{_evidence_file_name}" }, $evidence_item->{type});
				}
			}
			
			sub freq_to_string
			{
				my ($mut, $key) = @_;
				$key = "frequency" if (!defined $key);
				my $freq = $mut->{$key};
				return 'H' if ($freq eq 'H');
				return 1 if (!defined $freq);
				return '?' if ($freq eq '?');
				return '' if ($freq == 0);
				my $frequency_string = sprintf("%4.1f%%", $freq*100);
				$frequency_string = "100%" if ($freq == 1); # No "100.0%"
				return $frequency_string;
			}
			
			sub freq_cols
			{
				my @freq_list = @_;
				my $output_str = '';
				
				foreach my $freq (@freq_list)
				{
					$output_str .= td({align=>"right"}, $freq);
				}
				return $output_str;
			}
			
			my $row_class = "normal_table_row";
			my @freq_list;
			if (!defined $gd_name_list_ref)
			{
				if ((defined $mut->{frequency}) && ($mut->{frequency} != 1))
				{
					$row_class = "polymorphism_table_row";	
				}			
				push @freq_list, freq_to_string($mut);
			}
			else
			{
				#"_freq_[name]" keys were made				
				@freq_list = map { freq_to_string($mut, "frequency_$_")  } @$gd_name_list_ref;			
			}
			
			
			if ($mut->{type} eq 'SNP')
			{
				## additional formatting for some variables
				my $aa_codon_change = '';
				if (($mut->{snp_type} ne 'intergenic') && ($mut->{snp_type} ne 'noncoding') && ($mut->{snp_type} ne 'pseudogene'))
				{
					$aa_codon_change .= "$mut->{aa_ref_seq}$mut->{aa_position}$mut->{aa_new_seq}";
					## add color and underlining  

					my $codon_ref_seq = to_underline_red_codon($mut, 'codon_ref_seq');
					my $codon_new_seq = to_underline_red_codon($mut, 'codon_new_seq');

					sub to_underline_red_codon
					{
						my ($mut, $codon_key) = @_;			
						return '' if (!defined $mut->{$codon_key} || !defined $mut->{codon_position} || $mut->{codon_position} eq '');

						my $codon_string;
						my @codon_ref_seq_list = split //, $mut->{$codon_key};
						for (my $i=0; $i<scalar @codon_ref_seq_list; $i++)
						{
							if ($i == $mut->{codon_position})
							{
								$codon_string.= font({-class=>"mutation_in_codon"}, $codon_ref_seq_list[$i]);
							}
							else
							{
								$codon_string.= $codon_ref_seq_list[$i];
							}	
						}	
						return $codon_string;
					}

					$aa_codon_change .= "&nbsp;($codon_ref_seq&rarr;$codon_new_seq)";
				}
				else # ($mut->{snp_type} eq 'NC')
				{
					$aa_codon_change .= $mut->{gene_position};
				}

				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, "$mut->{_ref_seq}&rarr;$mut->{new_seq}");
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, nonbreaking($aa_codon_change));
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;	
			}
			elsif ($mut->{type} eq 'INS')
			{
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, "+$mut->{new_seq}");
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{gene_position}));
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;
			}
			elsif ($mut->{type} eq 'DEL')
			{
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, nonbreaking("&Delta;" . commify($mut->{size}) . " nt"));
				$output_str.= freq_cols(@freq_list);
				my $annotation_str = '';
				$annotation_str = "between $mut->{between}" if ($mut->{between});
				$annotation_str = "$mut->{mediated}-mediated" if ($mut->{mediated});
				$annotation_str = $mut->{gene_position} if (!$annotation_str); 
				$output_str.= td({align=>"center"}, nonbreaking($annotation_str));
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));				
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;
			}
			elsif ($mut->{type} eq 'SUB')
			{
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, nonbreaking("&Delta;$mut->{size} nt&rarr;$mut->{new_seq}"));
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{gene_position}));
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;
			}
			
			elsif ($mut->{type} eq 'CON')
			{
				
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, nonbreaking("$mut->{size} nt&rarr;$mut->{region}"));
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{gene_position}));
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;				
			}

			elsif ($mut->{type} eq 'MOB')
			{
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				my $s;

				my $s_start = '';
				$s_start .= "+" . $mut->{ins_start} if ($mut->{ins_start});
				$s_start .= "&Delta;" . $mut->{del_start} if ($mut->{del_start});
				$s.= $s_start . " :: " if ($s_start);
				
				$s .= "$mut->{repeat_name} (";
				$s .= (($mut->{strand}==+1) ? '+' : (($mut->{strand}==-1) ? '&minus;' : '?'));
				$s .= ")";

				my $s_end = '';
				$s_end .= "&Delta;" . $mut->{del_end} if ($mut->{del_end});
				$s_end .= "+" . $mut->{ins_end} if ($mut->{ins_end});
				$s.= " :: " . $s_end if ($s_end);

				my $dup_str = ($mut->{duplication_size} >= 0) ? "+$mut->{duplication_size}" : "&Delta;" . abs($mut->{duplication_size});
				$s .= " ($dup_str) bp";			
				$output_str.= td({align=>"center"}, nonbreaking($s));
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{gene_position}));
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;
			}
			
			
			elsif ($mut->{type} eq 'INV')
			{
				
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, nonbreaking(commify($mut->{size}) . " nt inversion"));
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, "");
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name_1})) . "&darr;" . i(nonbreaking($mut->{gene_name_2})) );
				$output_str.= td({align=>"left"}, $mut->{gene_product_1} . "&darr;" . $mut->{gene_product_2});
				$output_str.= end_Tr;				
			}
			
			elsif ($mut->{type} eq 'DUP')
			{
				$output_str.= start_Tr({-class=>$row_class});	
				$output_str.= td({align=>"center"}, $evidence_string) if (!defined $gd_name_list_ref);
				$output_str.= td({align=>"center"}, nonbreaking($mut->{seq_id})) if (!$one_ref_seq); 
				$output_str.= td({align=>"right"}, commify($mut->{position}));
				$output_str.= td({align=>"center"}, nonbreaking(commify($mut->{size}) . " nt duplication"));
				$output_str.= freq_cols(@freq_list);
				$output_str.= td({align=>"center"}, "");
				$output_str.= td({align=>"center"}, i(nonbreaking($mut->{gene_name})));
				$output_str.= td({align=>"left"}, $mut->{gene_product});
				$output_str.= end_Tr;				
			}
		}
		
		if ($legend_row)
		{
			$output_str.= start_Tr();	
			$output_str.= td({-colspan=>$total_cols}, b("Evidence codes: RA = read alignment, MC = missing coverage, JC = new junction"));
			$output_str.= end_Tr;	
		}
		
		$output_str.= end_table;	
}

sub html_read_alignment_table_string
{
	my ($list_ref, $relative_link, $title, $show_reject_reason) = @_;
	$relative_link = '' if (!$relative_link);
	$title = "Read alignment evidence..." if (!$title);
	
	my $output_str = '';
	
	my $q = new CGI;
	$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3});
	
	my $link = (defined $list_ref->[0]) && (defined $list_ref->[0]->{_evidence_file_name});
	my $total_cols = $link ? 11 : 10;
	$output_str.= Tr(th({-colspan => $total_cols, -align => "left", -class=>"read_alignment_header_row"}, $title));

	$output_str.= start_Tr();
	if ($link)
	{
		$output_str.= th("&nbsp;"); 
	}
	$output_str.= th("seq&nbsp;id");
	$output_str.= th({colspan => 2}, "position");
	
	$output_str.= th( [
			"change",
			"freq",
			"score", 
			"cov", 
			"annotation", 
			"genes", 
			
		]
	);	
	$output_str.= th({-width => "100%"}, "product"); 
	$output_str.= end_Tr;
	
	foreach my $c (@$list_ref)
	{			
		my $row_class = "normal_table_row";
		if ((defined $c->{frequency}) && ($c->{frequency} != 1))
		{
			$row_class = "polymorphism_table_row";	
		}
		$output_str.= start_Tr({-class=>$row_class});
		
		if ($link)
		{
			$output_str.= td(a({-href=>"$relative_link$c->{_evidence_file_name}"}, '*')); 
		}
		
		my $fisher_p_value = '';
		$fisher_p_value = nonbreaking(sprintf("&nbsp;(%.1E)", $c->{fisher_strand_p_value})) if (defined $c->{fisher_strand_p_value} && $c->{fisher_strand_p_value} ne 'ND');
		
		$output_str.= td({align => "center"}, nonbreaking($c->{seq_id}) );	
		$output_str.= td({align => "right"}, commify($c->{position}) );	
		$output_str.= td({align => "right"}, $c->{insert_position} );
		$output_str.= td({align => "center"}, "$c->{ref_base}&rarr;$c->{new_base}" );	
		$output_str.= td({align => "right"}, sprintf("%4.1f%%", $c->{frequency}*100) );
		$output_str.= td({align => "right"}, sprintf("%.1f", $c->{quality}) );	# . $fisher_p_value	
		my ($top_cov, $bot_cov) = split /\//, $c->{tot_cov};	
		$output_str.= td({align => "center"}, $top_cov + $bot_cov );
		$output_str.= td({align => "center"}, nonbreaking($c->{gene_position}) );	
		$output_str.= td({align => "center"}, i(nonbreaking($c->{gene_name})) );	
		$output_str.= td({align => "left"}, $c->{gene_product} );	
			
		$output_str.= end_Tr;
		
		if ($show_reject_reason) 
		{
			foreach my $reject (Breseq::GenomeDiff::get_reject_reasons($c))
			{
				$output_str.= Tr({-class=>'reject_table_row'}, td({-colspan => $total_cols}, "Rejected: " . decode_reject_reason($reject)));
			}
			
#			if (defined $c->{bias_p_value})
#			{
#				my $bias_p_value = sprintf("%.2E", $c->{bias_p_value});
#				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Combined strand and quality bias " . i("p") . "-value = $bias_p_value"));
#			}
			
			if (defined $c->{fisher_strand_p_value})
			{
				my $fisher_strand_p_value = sprintf("%.2E", $c->{fisher_strand_p_value});
				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, 
					"Strands of reads supporting (+/-):&nbsp;&nbsp;" 
					. b("new") . " base ($c->{new_cov})&nbsp;&nbsp;" 
					. b("ref") . " base ($c->{ref_cov})&nbsp;&nbsp;" 
					. b("total") . " ($c->{tot_cov})"));
				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Fisher's exact test strand distribution " . i("p") . "-value = $fisher_strand_p_value"));
			}

#			if (defined $c->{ks_quality_p_value})
#			{
#				my $ks_quality_p_value = sprintf("%.2E", $c->{ks_quality_p_value});
#				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Kolmogorov-Smirnov test that lower quality scores support polymorphism " . i("p") . "-value = $ks_quality_p_value"));
#			}

			if (defined $c->{ks_quality_p_value_unusual_poly})
			{
				my $ks_quality_p_value = sprintf("%.2E", $c->{ks_quality_p_value_unusual_poly});
				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Kolmogorov-Smirnov test of lower base quality scores supporting polymorphism  " . i("p") . "-value = $ks_quality_p_value"));
			}


			if (defined $c->{ks_quality_p_value_unusual_new})
			{
				my $ks_quality_p_value = sprintf("%.2E", $c->{ks_quality_p_value_unusual_new});
				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Kolmogorov-Smirnov test of lower quality scores supporting new bases " . i("p") . "-value = $ks_quality_p_value"));
			}

			if (defined $c->{ks_quality_p_value_unusual_ref})
			{
				my $ks_quality_p_value = sprintf("%.2E", $c->{ks_quality_p_value_unusual_ref});
				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Kolmogorov-Smirnov test of lower quality scores supporting ref bases " . i("p") . "-value = $ks_quality_p_value"));
			}
			
			if (defined $c->{ks_quality_p_value_unusual_all})
			{
				my $ks_quality_p_value = sprintf("%.2E", $c->{ks_quality_p_value_unusual_all});
				$output_str.= Tr({-class=>'information_table_row'}, td({-colspan => $total_cols}, "Kolmogorov-Smirnov test of lower quality scores supporting all bases " . i("p") . "-value = $ks_quality_p_value"));
			}			
			

		}
	}
	
	$output_str.= end_table;
}

sub html_missing_coverage_table_string
{
	my ($list_ref, $relative_link, $title, $show_reject_reason) = @_;
	$relative_link = '' if (!$relative_link);
	$title = "Missing coverage evidence..." if (!$title);
	
	my $output_str = '';
	
	my $q = new CGI;
	$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3, -width => "100%"});
	
	my $coverage_plots;
	$coverage_plots = (defined $list_ref->[0]) && (defined $list_ref->[0]->{_evidence_file_name});
	
	my $link = (defined $list_ref->[0]) && (defined $list_ref->[0]->{_side_1_evidence_file_name}) && (defined $list_ref->[0]->{_side_2_evidence_file_name});

	my $total_cols = $link ? 11 : 8;
	$output_str.= Tr(th({-colspan => $total_cols, -align => "left", -class=>"missing_coverage_header_row"}, $title));

	$output_str.= start_Tr();
	if ($link)
	{
		$output_str.= th("&nbsp;"); 
		$output_str.= th("&nbsp;"); 
		if ($coverage_plots)
		{
			$output_str.= th("&nbsp;");
		}
	}
	
	$output_str.= th(
		[
			"seq&nbsp;id",
			"start", 
			"end", 
			"size",
			"&larr;cov",
			"cov&rarr;",
			"gene", 
		]
	);	
	$output_str.= th({-width => "100%"}, "description");
	$output_str.= end_Tr;
	
	foreach my $c (@$list_ref)
	{
		## additional formatting for some variables
		my $gene_string = $c->{genes};
	#	$gene_string =~ s/ /&nbsp;/g;
	#	$gene_string =~ s/-/&#8209;/g; #substitute nonbreaking dash
		
		$output_str.= start_Tr;
		
		if ($link)
		{
			$output_str.= td(a({-href=>"$relative_link$c->{_side_1_evidence_file_name}"}, '*')); 
			$output_str.= td(a({-href=>"$relative_link$c->{_side_2_evidence_file_name}"}, '*')); 
			
			if ($coverage_plots)
			{
				$output_str.= td(a({-href=>"$relative_link$c->{_evidence_file_name}"}, '&divide;')); 
			}
		}

		my $start_str = $c->{start};
		$start_str .= "–" . ($c->{start} + $c->{start_range}) if ($c->{start_range} > 0);
		my $end_str = $c->{end};
		$end_str .= "–" . ($c->{end} - $c->{end_range}) if ($c->{end_range} > 0);

		my $size_str = ($c->{end} - $c->{start} + 1);
		$size_str = ($c->{end} - $c->{start} + 1 - $c->{end_range} - $c->{start_range}) . "–" . $size_str if (($c->{end_range} > 0) || ($c->{start_range} > 0));

				
		$output_str.= td(nonbreaking($c->{seq_id})); 
		$output_str.= td({-align=>"right"}, nonbreaking($start_str)); 
		$output_str.= td({-align=>"right"}, nonbreaking($end_str)); 	
		$output_str.= td({-align=>"right"}, nonbreaking($size_str)); 		
		$output_str.= td({-align=>"center"}, nonbreaking("$c->{left_outside_cov} \[$c->{left_inside_cov}\]")); 
		$output_str.= td({-align=>"center"}, nonbreaking("\[$c->{right_inside_cov}\] $c->{right_outside_cov}")); 
		$output_str.= td({align=>"center"}, i(nonbreaking($c->{gene_name})));				
		$output_str.= td({align=>"left"}, $c->{gene_product});		
		$output_str.= end_Tr;
		
		if ($show_reject_reason)
		{
			foreach my $reject (Breseq::GenomeDiff::get_reject_reasons($c))
			{
				$output_str.= Tr({-class=>'reject_table_row'}, td({-colspan => $total_cols}, "Rejected: " . decode_reject_reason($reject)));
			}
		}
	}
	
	$output_str.= end_table;
	return $output_str;
}

sub html_new_junction_table_string
{
	our ($list_ref, $relative_link, $title, $show_reject_reason) = @_;
	$relative_link = '' if (!$relative_link);
	$title = "New junction evidence..." if (!$title);
	
	my $output_str = '';

	my $test_item = $list_ref->[0];
	my $link =  
	       (defined $test_item) 
		&& (defined $test_item->{_side_1_evidence_file_name}) 
		&& (defined $test_item->{_side_2_evidence_file_name})
		&& (defined $test_item->{_new_junction_evidence_file_name})
	;
	
	my $q = new CGI;
	$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3});

	my $total_cols = $link ? 10 : 8;
	$output_str.= Tr(th({-colspan => $total_cols, -align => "left", -class=>"new_junction_header_row"}, $title));
		
	#####################
	#### HEADER LINE ####
	#####################
	
	$output_str.= start_Tr();
	if ($link)
	{
		$output_str.= th({-colspan=>2}, "&nbsp;"); 
	}
	$output_str.= th(
		[
			"seq&nbsp;id",
			"position",
			"overlap",
			"reads", 
			"score",
			"annotation",
			"gene",
		]
	);		
	$output_str.= th({-width => "100%"}, "product"); 
	$output_str.= end_Tr;
	
	####################
	#### ITEM LINES ####
	####################
	
	## the rows in this table are linked (same background color for every two)
	my $row_bg_color_index = 0;
	foreach my $c (@$list_ref)
	{			
		### Side 1
		my $key = 'side_1';			
		my $annotate_key = "junction_" . $c->{"$key\_annotate_key"};
		
		$output_str.= start_Tr({-class=> "mutation_table_row_$row_bg_color_index"});
		$output_str.= td({-rowspan=>2}, a({-href=>"$relative_link$c->{_new_junction_evidence_file_name}"}, "*")) if ($link); 
		{	
			$output_str.= td({-rowspan=>1}, a({-href=>"$relative_link$c->{_side_1_evidence_file_name}"}, "?")) if ($link); 
			$output_str.= td({-rowspan=>1, -class=>"$annotate_key"}, nonbreaking($c->{"$key\_seq_id"}));			
			$output_str.= td( {-align=>"center", -class=>"$annotate_key"}, ($c->{"$key\_strand"} == +1) ? $c->{"$key\_position"} . "&nbsp;=": "=&nbsp;" . $c->{"$key\_position"} );
			$output_str.= td( {-rowspan=>2, -align=>"center"}, $c->{overlap} );
			$output_str.= td( {-rowspan=>2, -align=>"center"}, $c->{total_reads} );
			$output_str.= td( {-rowspan=>2, -align=>"center"}, $c->{score} );
			$output_str.= td( {-align=>"center", -class=>"$annotate_key"}, nonbreaking($c->{"_$key"}->{gene_position}) );
			$output_str.= td( {-align=>"center", -class=>"$annotate_key"}, i(nonbreaking($c->{"_$key"}->{gene_name})) );
			$output_str.= td( {-class=>"$annotate_key"}, $c->{"_$key"}->{gene_product} );
		}
		$output_str.= end_Tr;

		### Side 2
		$key = 'side_2';
		$annotate_key = "junction_" . $c->{"$key\_annotate_key"};
		
		$output_str.= start_Tr({-class=> "mutation_table_row_$row_bg_color_index"});		
		{
			$output_str.= td({-rowspan=>1}, a({-href=>"$relative_link$c->{_side_2_evidence_file_name}"}, "?")) if ($link); 
			$output_str.= td({-rowspan=>1, -class=>"$annotate_key"}, nonbreaking($c->{"$key\_seq_id"}));		
			$output_str.= td( {-align=>"center", -class=>"$annotate_key"}, ($c->{"$key\_strand"} == +1) ? $c->{"$key\_position"} . "&nbsp;=": "=&nbsp;" . $c->{"$key\_position"} );

			$output_str.= td( {-align=>"center", -class=>"$annotate_key"}, nonbreaking($c->{"_$key"}->{gene_position}) );
			$output_str.= td( {-align=>"center", -class=>"$annotate_key"}, i(nonbreaking($c->{"_$key"}->{gene_name})) );
			$output_str.= td( {-class=>"$annotate_key"}, $c->{"_$key"}->{gene_product} );
		}			
		$output_str.= end_Tr;
		
		if ($show_reject_reason)
		{
			foreach my $reject (Breseq::GenomeDiff::get_reject_reasons($c))
			{
				$output_str.= Tr({-class=>'reject_table_row'}, td({-colspan => $total_cols}, "Rejected: " . decode_reject_reason($reject)));
			}
		}
		
		$row_bg_color_index = ($row_bg_color_index+1)%2;
	}
	
	$output_str.= end_table;

}

sub html_evidence_file_name
{
	my ($interval) = @_;
	
	#set up the file name
	my $html_evidence_file_name = 
		"$interval->{prefix}"
		. "_$interval->{seq_id}"
		. "_$interval->{start}" 
		. ((defined $interval->{insert_start}) ?  ".$interval->{insert_start}" : '')
		. "_$interval->{end}"
		. ((defined $interval->{insert_end}) ?  ".$interval->{insert_end}" : '')
		. "_alignment.html";
		
	return $html_evidence_file_name;
}

sub html_evidence_file
{
	my ($settings, $gd, $interval) = @_;

 	$interval->{output_path} = $settings->file_name('evidence_path') . "/$interval->{file_name}";	
	
	my $title = 'BRESEQ :: Results' . ($settings->{run_name} ne 'unnamed' ? " :: $settings->{run_name}" : ''),
	
	open HTML, ">$interval->{output_path}" or die "Could not open file: $interval->{output_path}";

	my $q = new CGI;

	print HTML
		start_html(
			-title => $title, 
			-head  => style({type => 'text/css'},$header_style_string),
	    );

	## print a table for the main item
	## followed by auxiliary tables for each piece of evidence
	
	my $parent_item = $interval->{parent_item};
				
	print HTML html_genome_diff_item_table_string($gd, [$parent_item]) . p;
	my @evidence_list = $gd->mutation_evidence_list($parent_item);
	
	foreach my $type ( 'RA', 'MC', 'JC' )
	{
		my @this_evidence_list = grep {$_->{type} eq $type} @evidence_list;
		next if (scalar @this_evidence_list == 0);
		print HTML html_genome_diff_item_table_string($gd, \@this_evidence_list) . p;
	}
		
	if (defined $interval->{plot})
	{
		print HTML div({-align=>"center"}, img({-src=>$interval->{plot}}));
	}
	elsif ( (defined $interval->{bam_path}) && ($interval->{fasta_path}) )
	{					
		#construct the interval string		
		my $s = '';
		$s .= "$interval->{seq_id}:$interval->{start}";
		$s .= ".$interval->{insert_start}" if (defined $interval->{insert_start});
		$s .= "-$interval->{end}";
		$s .= ".$interval->{insert_end}" if (defined $interval->{insert_end});
		
		my $ao = Breseq::AlignmentOutput->new;
		my $options = {};
		$interval->{quality_score_cutoff} = $settings->{base_quality_cutoff} if (defined $settings->{base_quality_cutoff});
		
		print STDERR "Creating read alignment for region \'$s\'\n";
		print HTML $ao->html_alignment(
			$interval->{bam_path}, 
			$interval->{fasta_path}, 
			$s, 
			$interval,
		);
	}
	print HTML end_html;
	close HTML;
}


sub decode_reject_reason
{
	my ($reject) = @_;
	
	if ($reject eq 'NJ')
	{
		return "Insufficient overlap of new junction sides by reads on both strands.";
	}
	elsif ($reject eq 'EVALUE')
	{
		return "E-value exceeds prediction threshold.";
	}
	elsif ($reject eq 'STRAND')
	{
		return "Prediction not supported by reads on both strands.";
	}
	elsif ($reject eq 'FREQ')
	{
		return "Prediction has frequency below cutoff threshold.";
	}	
	elsif ($reject eq 'COV')
	{
		return "Prediction has coverage below cutoff threshold.";
	}
	elsif ($reject eq 'BIAS_P_VALUE')
	{
		return "Prediction has biased strand and/or quality scores supporting polymorphism.";
	}
	elsif ($reject eq 'KS_QUALITY_P_VALUE_UNUSUAL_POLY')
	{
		return "Prediction has biased quality score distribution for polymorphism bases.";
	}
	elsif ($reject eq 'KS_QUALITY_P_VALUE_UNUSUAL_REF')
	{
		return "Prediction has biased quality score distribution for new bases.";
	}
	elsif ($reject eq 'KS_QUALITY_P_VALUE_UNUSUAL_NEW')
	{
		return "Prediction has biased quality score distribution for ref bases.";
	}
	elsif ($reject eq 'KS_QUALITY_P_VALUE_UNUSUAL_ALL')
	{
		return "Prediction has biased quality score distribution for all bases.";
	}
	elsif ($reject eq 'FISHER_STRAND_P_VALUE')
	{
		return "Prediction has biased read strand distribution supporting polymorphism.";
	}	
	elsif ($reject eq 'POLYMORPHISM_STRAND')
	{
		return "Polymorphism prediction not supported by minimum number of reads on both strands.";
	}
	
	return '';
}

sub create_evidence_files
{
	my ($settings, $gd) = @_;
	
	# gather everything together and then handle all at once
	our @evidence_list;


	sub add_evidence
	{
		my ($evidence_file_name_key, $evidence_item) = @_;		
		$evidence_item->{file_name} = html_evidence_file_name($evidence_item);
		$evidence_item->{item}->{$evidence_file_name_key} = $evidence_item->{file_name};		
		push @evidence_list, $evidence_item;
	}
	
	# Fasta and BAM files for making alignments.
	my $reference_bam_file_name = $settings->file_name('reference_bam_file_name');
	my $reference_fasta_file_name = $settings->file_name('reference_fasta_file_name');

	## hybrids use different BAM files for making the alignments!!!
	my $junction_bam_file_name = $settings->file_name('junction_bam_file_name');
	my $junction_fasta_file_name = $settings->file_name('candidate_junction_fasta_file_name');

	### We make alignments of two regions for deletions: upstream and downstream edges.
	foreach my $item ( $gd->list('MC') )
	{	
		next if ($item->{no_show});
		
		my $parent_item = $gd->parent($item);
		$parent_item = $item if (!$parent_item);
		
		add_evidence( 
			'_side_1_evidence_file_name',
			{ 
				bam_path 	=> $reference_bam_file_name,
				fasta_path 	=> $reference_fasta_file_name,
				seq_id 		=> $item->{seq_id}, 
				start 		=> $item->{start}-1, 
				end 		=> $item->{start}-1, 
				parent_item => $parent_item,
				item 		=> $item,
				prefix 		=> 'MC_SIDE_1',
			}
		);
		
		add_evidence( 
			'_side_2_evidence_file_name',
			{
				bam_path 	=> $reference_bam_file_name,
				fasta_path 	=> $reference_fasta_file_name,
				seq_id 		=> $item->{seq_id}, 
				start		=> $item->{end}+1, 
				end 		=> $item->{end}+1, 
				parent_item => $parent_item,
				item 		=> $item,
				prefix 		=> 'MC_SIDE_2',			
			}
		);

		add_evidence( 
			'_evidence_file_name',
			{
				seq_id 		=> $item->{seq_id}, 
				start		=> $item->{start}, 
				end 		=> $item->{end}, 
				parent_item => $parent_item,
				item 		=> $item,
				prefix 		=> 'MC_PLOT',	
				plot		=> $item->{_coverage_plot_file_name},
			}
		);
	}	

	MUT: foreach my $item ( $gd->list('SNP', 'INS', 'DEL', 'SUB') )
	{
		next if ($item->{no_show});
		
		#this reconstructs the proper columns to draw
		my $start = $item->{position};
		my $end = $start;
		my $insert_start = undef;
		my $insert_end = undef;

		if ($item->{type} eq 'INS')
		{
			$insert_start = 1;
			$insert_end = length($item->{new_seq});			
		}
		elsif ($item->{type} eq 'DEL')
		{
			my $has_ra_evidence;
			foreach my $evidence_item ($gd->mutation_evidence_list($item))
			{
				$has_ra_evidence = 1 if ($evidence_item->{type} eq 'RA');
			}
			## only do deletions if they have within-read evidence
			next MUT if (!$has_ra_evidence);

			$end = $start + $item->{size} - 1;
		}

		##may be a problem here...
		elsif ($item->{type} eq 'SUB')
		{
			$end = $start + length($item->{new_seq}) - 1;
		}

		add_evidence( 
			'_evidence_file_name',
			{
				bam_path 		=> $reference_bam_file_name,
				fasta_path 		=> $reference_fasta_file_name,
				seq_id 			=> $item->{seq_id}, 
				start 			=> $start, 
				end 			=> $end, 
				insert_start 	=> $insert_start, 
				insert_end 		=> $insert_end,
				parent_item 	=> $item, 
				item 			=> $item, 
				prefix 			=> $item->{type},	
			}
		);
				
		#add evidence to 'RA' items as well
		foreach my $evidence_item ( $gd->mutation_evidence_list($item) )
		{
			next if ($evidence_item->{type} ne 'RA');
			$evidence_item->{_evidence_file_name} = $item->{_evidence_file_name};
		}
	}
	
	
	## Still create files for RA evidence that was not good enough to predict a mutation from

	my @ra_list = $gd->list('RA');	
	@ra_list = $gd->filter_used_as_evidence(@ra_list);
	
	RA: foreach my $item ( @ra_list )
	{
		next if ($item->{no_show});
		
		#this reconstructs the proper columns to draw
		add_evidence( 
			'_evidence_file_name',
			{
				bam_path 		=> $reference_bam_file_name,
				fasta_path 		=> $reference_fasta_file_name,
				seq_id 			=> $item->{seq_id}, 
				start 			=> $item->{position}, 
				end 			=> $item->{position}, 
				insert_start 	=> $item->{insert_position}, 
				insert_end 		=> $item->{insert_position},
				parent_item 	=> $item, 
				item 			=> $item, 
				prefix 			=> $item->{type},	
			}
		);
	}

	## This additional information is used for the complex reference line.
	## Note that it is completely determined by the original candidate junction sequence 
	## positions and overlap: alignment_pos and alignment_overlap.

	foreach my $item ( $gd->list('JC') )
	{	
		next if ($item->{no_show});
		
		my $parent_item = $gd->parent($item);
		$parent_item = $item if (!$parent_item);
		
		## regenerate the alignment overlap from the junction_key
		my ($start, $end);
		if ($item->{alignment_overlap} == 0)
		{
			$start = $item->{flanking_left};
			$end = $item->{flanking_left}+1;			
		}
		elsif ($item->{alignment_overlap} > 0)
		{
			$start = $item->{flanking_left}+1;
			$end = $item->{flanking_left}+$item->{alignment_overlap};
		}
		else ## ($item->{overlap} < 0)
		{
			$start = $item->{flanking_left}+1;
			$end = $item->{flanking_left}-$item->{alignment_overlap};
		}

		add_evidence( 
			'_new_junction_evidence_file_name',
			{
				bam_path 	=> $junction_bam_file_name,
				fasta_path 	=> $junction_fasta_file_name,
				seq_id 		=> $item->{key}, 
				start 		=> $start, 
				end 		=> $end, 
				parent_item => $parent_item,
				item 		=> $item,
				prefix 		=> 'JC',
			#### extra information	
				alignment_empty_change_line => 1,			
				alignment_reference_info_list => [
				 	{	
						truncate_end 	=> $item->{flanking_left} + $item->{side_1_overlap}, 
						ghost_end 		=> $item->{side_1_position}, 
						ghost_strand 	=> $item->{side_1_strand},
						ghost_seq_id	=> $item->{side_1_seq_id}
					},
					{	
						truncate_start 	=> $item->{flanking_left} + 1 + abs($item->{alignment_overlap}) - $item->{side_2_overlap}, 
						ghost_start 	=> $item->{side_2_position} , 
						ghost_strand 	=> $item->{side_2_strand},
						ghost_seq_id 	=> $item->{side_2_seq_id}
					}
				],
			}
		);
		## this is the flagship file that we show first when clicking on evidence from a mutation...
		$item->{_evidence_file_name} = $item->{_new_junction_evidence_file_name};
		
		add_evidence( 
			'_side_1_evidence_file_name',
		 	{ 
				bam_path 	=> $reference_bam_file_name,
				fasta_path 	=> $reference_fasta_file_name,
				seq_id 		=> $item->{side_1_seq_id}, 
				start 		=> $item->{side_1_position}, 
				end 		=> $item->{side_1_position}, 
				parent_item => $parent_item,
				item 		=> $item,
				prefix 		=> 'JC_SIDE_1' . "_$item->{side_2_seq_id}_$item->{side_2_position}_$item->{side_2_position}",		## need to be unique
			}
		);

		add_evidence( 
			'_side_2_evidence_file_name',
			{
				bam_path 	=> $reference_bam_file_name,
				fasta_path 	=> $reference_fasta_file_name,
				seq_id 		=> $item->{side_2_seq_id}, 
				start 		=> $item->{side_2_position}, 
				end 		=> $item->{side_2_position}, 
				parent_item => $parent_item,
				item 		=> $item, 
				prefix 		=> 'JC_SIDE_2' . "_$item->{side_1_seq_id}_$item->{side_1_position}_$item->{side_1_position}",
			}
		);
	}

	### now create evidence files
	$settings->create_path('evidence_path');
	print STDERR "Creating HTML evidence files...\n";
	foreach my $e (@evidence_list)
	{			
		print STDERR "Creating evidence file: $e->{file_name}\n" if ($settings->{verbose});
		Breseq::Output::html_evidence_file($settings, $gd, $e);		
	}
		
}

sub save_text_deletion_file
{
	my ($deletion_file_name, $deletions_ref) = @_;

	open DEL, ">$deletion_file_name" or die "Could not open: $deletion_file_name";
	print DEL join("\t", 'seq_id', 'start', 'end') . "\n";
	foreach my $d (@$deletions_ref)
	{
		print DEL join("\t", $d->{seq_id}, $d->{start}, $d->{end}) . "\n"; 
	}
	close DEL;
}


sub draw_coverage
{
	my ($settings, $ref_seq_info, $gd) = @_;
	my @mc = $gd->list('MC');
	my $drawing_format = 'png';
	
	$settings->create_path('coverage_plot_path');
	my $coverage_plot_path = $settings->file_name('coverage_plot_path');	
	my $deletions_text_file_name = $settings->file_name('deletions_text_file_name');
	Breseq::Output::save_text_deletion_file($deletions_text_file_name, \@mc);

	foreach my $seq_id (@{$ref_seq_info->{seq_ids}})
	{
		my $this_complete_coverage_text_file_name = $settings->file_name('complete_coverage_text_file_name', {'@'=>$seq_id});			
		my $res = Breseq::Shared::system("$FindBin::Bin/plot_coverage --drawing-format $drawing_format -t $coverage_plot_path -p $settings->{coverage_plot_path} -i $deletions_text_file_name -c $this_complete_coverage_text_file_name --seq_id=$seq_id");				
		die if ($res);
		
		#need to assign link names that correspond to what the R script is doing
		my $i=1;
		my @this_deletions = grep {$_->{seq_id} eq $seq_id} @mc if ($seq_id);
		foreach my $del (@this_deletions)
		{
			$del->{_coverage_plot_file_name} = "$seq_id\.$i\.$drawing_format";
			$i++;
		}
	}
	$settings->remove_path('deletions_text_file_name');
}

our @execution_times;
sub record_time
{
	my ($name) = @_;
	
	my $this_time = time;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($this_time);	
	my $formatted_time = sprintf "%04d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec;
	my $new_time = { 
		_name => $name, 
		_formatted_time => $formatted_time, 
		_time => $this_time,
		_time_elapsed => 0,
		_formatted_time_elapsed => '',
	 };
	
	##if we had a previous time
	my $time_since_last;
	my $time_since_last_formatted;
	if (scalar @execution_times > 0)
	{
		$time_since_last = $this_time - $execution_times[-1]->{_time};
		$new_time->{_time_elapsed} = $time_since_last;
		$new_time->{_formatted_time_elapsed} = time2string($time_since_last);
	}
	
	push @execution_times, $new_time;
	return $formatted_time;
}

sub time2string
{
    my ($seconds) = @_;
    # Convert seconds to days, hours, minutes, seconds
    my @parts = gmtime($seconds);
    my $ret = '';
    if(sprintf("%4d",$parts[7])>0)
    {
        $ret .= sprintf("%4d",$parts[7]);
        $ret .= sprintf(" %s",($parts[7]>1)?'days':'day');
    }
    if(sprintf("%4d",$parts[2])>0)
    {
        $ret .= sprintf("%4d",$parts[2]);
        $ret .= sprintf(" %s",($parts[2]>1)?'hours':'hour');
    }
    if(sprintf("%4d",$parts[1])>0)
    {
        $ret .= sprintf("%4d",$parts[1]);
        $ret .= sprintf(" %s",($parts[1]>1)?'minutes':'minute');
    }
    if(sprintf("%4d",$parts[0])>0)
    {
        $ret .= sprintf("%4d",$parts[0]);
        $ret .= sprintf(" %s",($parts[0]>1)?'seconds':'second');
    }
    return $ret;
}


## Rather than dealing with our own formats for saving 
## many different kinds of summary statistics,
## just use freeze and thaw to directly save and load
## the perl data structures (generally hashes).
use Storable;
sub save_statistics
{
	my ($file_name, $data) = @_;
	store $data, "$file_name";
}

sub load_statistics
{
	my ($file_name) = @_;
	return retrieve($file_name);
}

sub nonbreaking
{
	my ($text) = @_;
	$text =~ s/-/&#8209;/g; #substitute nonbreaking hyphen
	$text =~ s/–/&#8211;/g; #substitute nonbreaking en dash
	$text =~ s/ /&nbsp;/g; #substitute nonbreaking space
	return $text;
}

sub commify 
{
	my ($input)  = @_;
	$input =~ s{(?<!\d|\.)(\d{4,})}
	{
		my $n = $1;
     	$n=~s/(?<=.)(?=(?:.{3})+$)/,/g;
     	$n;
    }eg;
   return $input;
}



return 1;

