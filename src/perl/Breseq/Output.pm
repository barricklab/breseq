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

.mutation_table_row_0	{
	background-color: rgb(245,245,245);
}
.mutation_table_row_1	{
	background-color: rgb(215,215,215);
}	

.polymorphism_table_row	{
	background-color: rgb(160,255,160);
}

.highlight_table_row	{
	background-color: rgb(192,255,255);
}

.is_insertion_header	{
	background-color: rgb(0,0,255);
	color: rgb(255,255,255);
}
.is_insertion_body	{
	background-color: rgb(192,192,255);
}
.junction_is	{
	background-color: rgb(255,165,0); /* orange */
}
.junction_gene	{
}

ENDOFSTYLE


sub html_index
{
	my ($file_name, $settings, $summary, $ref_seq_info, $annotated_mutations, $annotated_marginal) = @_;

	#copy over the breseq_graphic
	my $breseq_graphic_from_file_name = $settings->file_name('breseq_graphic_from_file_name');
	my $breseq_graphic_to_file_name = $settings->file_name('breseq_graphic_to_file_name');
	system("cp $breseq_graphic_from_file_name $breseq_graphic_to_file_name");

	open HTML, ">$file_name" or die "Could not open file: $file_name";    

    print HTML start_html(
			-title => "BRESEQ :: Index" . ($settings->{run_name} ? " :: $settings->{run_name}" : ''), 
			-head  => style({type => 'text/css'}, $header_style_string),
	);


	print HTML p;
	print HTML start_table({-width => "100%", -border => 0, -cellspacing => 0, -cellpadding => 3});
	print HTML start_Tr;
	print HTML  td(img({-src=>$settings->html_path('breseq_graphic_to_file_name')}));
	print HTML start_td({-width => "100%"});
	print HTML start_div({-style=>"font-size: 14pt;"});
	print HTML b(a({-href=>$settings->html_path('mutations_html_file_name')}, 'mutation predictions')) . br if ($annotated_mutations);
	print HTML b(a({-href=>$settings->html_path('marginal_html_file_name')}, 'marginal predictions')) . br if ($annotated_marginal);
	print HTML b(a({-href=>$settings->html_path('log_file_name')}, 'command line log'));
	print HTML end_div;

	print HTML end_td . end_Tr . end_table;

	## Write fastq read file information
	print HTML p;
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
	
	print HTML hr . $settings->{byline};	
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

sub html_full_table
{
	my ($file_name, $settings, $ref_seq_info, $snps_list_ref, $deletion_list_ref, $hybrid_list_ref) = @_;
	open HTML, ">$file_name" or die "Could not open file: $file_name";
	
	print HTML 
		start_html(
			-title => "BRESEQ :: Predicted Mutations" . ($settings->{run_name} ?  " :: $settings->{run_name}" : ''), 
			-head  => style({type => 'text/css'},$header_style_string),
	    ),
		h2("Within-Read Mutations"),
		html_snp_table_string($snps_list_ref),
		p,
		h2("Missing-Coverage Deletions"),
		html_deletion_table_string($deletion_list_ref),
		p,
		h2("Mosaic-Read Junctions"),
		html_junction_table_string($settings, $hybrid_list_ref, $ref_seq_info),
		end_html;

	close HTML;
}

sub html_snp_table_string
{
	my ($snps_list_ref, $force_link) = @_;
	my $output_str = '';
	
	my $q = new CGI;
	$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3});
	$output_str.= start_Tr();
	
	my $link = $force_link;
	$link = defined $snps_list_ref->[0] && defined $snps_list_ref->[0]->{link} if (!defined $link);
	if ($link)
	{
		$output_str.= th("&nbsp;"); 
	}
	$output_str.= th(
		[
			"seq",
			"pos", 
			"mut", 
			"score", 
#			"cov", 
#			"tot_cov", 
			"gene position", 
			"codon change", 
			"aa change", 
			"gene", 
		]
	);	
	$output_str.= th({-width => "100%"}, "product"); 
	$output_str.= end_Tr;
	
	foreach my $c (@$snps_list_ref)
	{	
		## additional formatting for some variables
		my $position = '';
		$position .= $c->{gene_position} if ($c->{gene_position});
		$position .= "&nbsp;($c->{aa_position})" if ($c->{aa_position});
		$position =~ s/-/&#8209;/g; #nonbreaking dash

		my $aa_change = ($c->{aa_ref_seq}) ? "$c->{aa_ref_seq}&rarr;$c->{aa_new_seq}" : "-";
		
		## add color and underlining  
		
		my $codon_ref_seq = to_underline_red_codon($c, 'codon_ref_seq');
		my $codon_new_seq = to_underline_red_codon($c, 'codon_new_seq');

		sub to_underline_red_codon
		{
			my ($c, $codon_key) = @_;			
			return '' if (!defined $c->{$codon_key} || !defined $c->{codon_position} || $c->{codon_position} eq '');
			
			my $codon_string;
			my @codon_ref_seq_list = split //, $c->{$codon_key};
			for (my $i=0; $i<scalar @codon_ref_seq_list; $i++)
			{
				if ($i == $c->{codon_position})
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
		
		my $codon_change = ($codon_ref_seq) ? "$codon_ref_seq&rarr;$codon_new_seq" : "-"; 
	
		my $display_quality = '';
		$display_quality = $c->{display_quality} if ($c->{display_quality});
		$display_quality =~ s/,/<br>/g;

		my $best_coverage_string = '';
		$best_coverage_string = $c->{best_coverage_string} if ($c->{best_coverage_string});
		$best_coverage_string =~ s/,/<br>/g;
	
		my $total_coverage_string = '';
		$total_coverage_string = $c->{total_coverage_string} if ($c->{total_coverage_string});
		$total_coverage_string =~ s/,/<br>/g;
	
		if ($c->{polymorphism})
		{
			$output_str.= start_Tr({-class=>"polymorphism_table_row"});	
		}
		else
		{		
			$output_str.= start_Tr;	
		}
		
		if ($link)
		{
			$output_str.= td(a({-href=>"$c->{link}"}, '*')); 
		}
		
		if ($c->{polymorphism})
		{
						
			my $display_frequency = sprintf "%4.2f", (1-$c->{frequency})*100;
			my $display_fisher_p_value = sprintf "%.1E", $c->{fisher_strand_p_value};
			
			$output_str.= td(make_nonbreaking($c->{seq_id}));
			$output_str.= td({-align=>"right"}, $c->{gene_shifted_start});
			$output_str.= td(
				[
					code("$c->{ref_seq}&rarr;$c->{new_seq}") . "&nbsp;(\%$display_frequency)",
#					code("$c->{ref_seq}&rarr;$c->{new_seq}") . "&nbsp;(FR=$display_frequency), SFET=$display_fisher_p_value)", 					 
					$display_quality, 
#					$best_coverage_string,
#					$total_coverage_string,
					$position, 
					code($codon_change), 
					code($aa_change), 
					i($c->{gene}), 
					$c->{gene_product}
		
				]
			);		
		}
		else
		{
			$output_str.= td(make_nonbreaking($c->{seq_id}));
			$output_str.= td({-align=>"right"}, $c->{gene_shifted_start});
			$output_str.= td(
				[
					code("$c->{ref_seq}&rarr;$c->{new_seq}"), 
					$display_quality, 
#					$best_coverage_string,
#					$total_coverage_string,
					$position, 
					code($codon_change), 
					code($aa_change), 
					i($c->{gene}), 
					$c->{gene_product}
		
				]
			);
		}
		$output_str.= end_Tr;
	}
	
	$output_str.= end_table;
}

sub html_deletion_table_string
{
	my ($list_ref, $force_link) = @_;
	my $output_str = '';
	
	my $q = new CGI;
	$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3});
	$output_str.= start_Tr();
	
	my $coverage_graphs;
	$coverage_graphs = defined $list_ref->[0] && defined $list_ref->[0]->{coverage_graph_link};
	
	my $link = $force_link;
	if (!defined $link && defined $list_ref->[0])
	{
		$link = (defined $list_ref->[0]->{upstream_interval}->{link}) && (defined $list_ref->[0]->{downstream_interval}->{link});
	}
	if ($link)
	{
		$output_str.= th("&nbsp;"); 
		$output_str.= th("&nbsp;"); 
		if ($coverage_graphs)
		{
			$output_str.= th("&nbsp;");
		}
	}
	
	$output_str.= th(
		[
			"seq",
			"start", 
			"end", 
			"size", 
		]
	);		
	$output_str.= th({-width => "100%"}, "genes"); 
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
			$output_str.= td(a({-href=>"$c->{upstream_interval}->{link}"}, '*')); 
			$output_str.= td(a({-href=>"$c->{downstream_interval}->{link}"}, '*')); 
			
			if ($coverage_graphs)
			{
				$output_str.= td(a({-href=>"$c->{coverage_graph_link}"}, '&divide;')); 
			}
		}
		
		$output_str.= td(make_nonbreaking($c->{seq_id})); 
		$output_str.= td({-align=>"right"}, $c->{start}); 
		$output_str.= td({-align=>"right"}, $c->{end}); 
		$output_str.= td({-align=>"right"}, $c->{size}); 
		$output_str.= td(i($gene_string)); 
			
		$output_str.= end_Tr;
	}
	
	$output_str.= end_table;
}

sub html_junction_table_string
{
	our ($settings, $list_ref, $ref_seq_info, $force_link) = @_;
	my $output_str = '';
		
	## only annotate IS elements for real reads...
	if (!$settings->{marginal_mode})
	{
		for (my $i=0; $i<scalar @$list_ref - 1; $i++)
		{
			my $j1 = $list_ref->[$i];
			my $j2 = $list_ref->[$i+1];		
								
			#must be same IS
			next if (!defined $j1->{is_interval} || !defined $j2->{is_interval});
			next if ($j1->{is_interval}->{is}->{gene} ne $j2->{is_interval}->{is}->{gene});
				
			#must have a non-negative overlap
#			next if ($j1->{overlap} < 0);
#			next if ($j2->{overlap} < 0);

			## positive overlap should be resolved by now
			die if ($j1->{overlap} > 0);
			die if ($j2->{overlap} > 0);

			#must be close together in real coords
			next if (abs($j1->{unique_interval}->{start} - $j2->{unique_interval}->{start}) > 20);
		
			#the first unique coords are going into the IS element
			my $uc1_strand = $j1->{unique_interval}->{strand};
			my $uc2_strand = $j2->{unique_interval}->{strand};
			next if ($uc1_strand != -$uc2_strand);
		
			my $is1_strand = - $j1->{is_interval}->{strand} * $j1->{is_interval}->{is}->{strand} * $j1->{unique_interval}->{strand};
			my $is2_strand = - $j2->{is_interval}->{strand} * $j2->{is_interval}->{is}->{strand} * $j2->{unique_interval}->{strand};
			my $is_strand = ($is1_strand == $is2_strand) ? $is2_strand : '0';
		
			### add additional information to the first match, which will 
			### cause a new line to be drawn in the new junction table
		
			$j1->{is_insertion}->{start} = ($uc1_strand == -1) ? $j2->{unique_interval}->{start} : $j1->{unique_interval}->{start};
			$j1->{is_insertion}->{end} = ($uc1_strand == -1) ? $j1->{unique_interval}->{start} : $j2->{unique_interval}->{start};
			$j1->{is_insertion}->{seq_id} = (defined $j1->{side_1}->{is}) ? $j1->{side_2}->{seq_id} : $j1->{side_1}->{seq_id};
			
			$j1->{is_insertion}->{family} = $j1->{is_interval}->{is}->{gene};
			$j1->{is_insertion}->{is_strand} = $is_strand;
			$j1->{is_insertion}->{is_size} = abs($j1->{is_interval}->{is}->{end} - $j1->{is_interval}->{is}->{start} + 1);
			$j1->{is_insertion}->{after_pos} = $j1->{is_insertion}->{start} - 1;
			$j1->{is_insertion}->{dup_size} = abs($j1->{is_insertion}->{end} - $j1->{is_insertion}->{start}) + 1;
		
			#sometimes the ends of the IS are not quite flush		
			if ($j1->{is_interval}->{strand} == -1)
			{
				$j1->{is_insertion}->{gap_left} = $j1->{is_interval}->{start} - $j1->{is_interval}->{is}->{end} + abs($j1->{overlap});
			}
			else
			{
				$j1->{is_insertion}->{gap_left} = $j1->{is_interval}->{is}->{start} - $j1->{is_interval}->{start} + abs($j1->{overlap});
			}
		
			if ($j2->{is_interval}->{strand} == -1)
			{
				$j1->{is_insertion}->{gap_right} = $j2->{is_interval}->{start} - $j2->{is_interval}->{is}->{end} + abs($j2->{overlap});
			}
			else
			{
				$j1->{is_insertion}->{gap_right} = $j2->{is_interval}->{is}->{start} - $j2->{is_interval}->{start} + abs($j2->{overlap});
			}
			
			if ($j1->{unique_interval}->{strand} *  $j1->{unique_interval}->{read_side} == +1)
			{
				($j1->{is_insertion}->{gap_right}, $j1->{is_insertion}->{gap_left}) = ($j1->{is_insertion}->{gap_left}, $j1->{is_insertion}->{gap_right});
			}
		}
	}
	
	my $q = new CGI;
	$output_str.= start_table({-border => 0, -cellspacing => 1, -cellpadding => 3});
	$output_str.= start_Tr();
	
	my $link = $force_link;
	if (!defined $link && defined $list_ref->[0])
	{
		$link = ( (defined $list_ref->[0]->{side_1}->{link}) && (defined $list_ref->[0]->{side_2}->{link}));
	}
	
	if ($link)
	{
		$output_str.= th({-colspan=>2}, "&nbsp;"); 
	}
	$output_str.= th({-colspan=>2}, "position"); 		
	$output_str.= th(
		[
			"overlap",
			"reads", 
			"gene",
			"coords",
		]
	);		
	$output_str.= th({-colspan=>2, -width => "100%"}, "product"); 
	$output_str.= end_Tr;
	
	## the rows in this table are linked (same background color for every two)
	my $row_bg_color_index = 0;
	
	my $body_countdown = 0;
	foreach my $c (@$list_ref)
	{
		#print STDERR "Hybrid\n";
		#print STDERR Dumper($c);

		## additional formatting for some variables
	#	my $gene_string = $c->{gene_string};
	#	$gene_string =~ s/ /&nbsp;/g;
		
		#special line for IS insertions
		if ($c->{is_insertion})
		{
			my $isi = $c->{is_insertion};
			
			$output_str.= start_Tr({-class=> "is_insertion_header"});
			$output_str.= td({-colspan=>1}, '') if ($link); 
			$output_str.= td({-colspan=>1}, '');
			$output_str.= td({-colspan=>1, -align=>"center"}, $isi->{seq_id});
			$output_str.= td({-colspan=>1, -align=>"center"}, $isi->{after_pos});
			$output_str.= td({-colspan=>1, -align=>"center"}, "+$isi->{dup_size} bp");
			$output_str.= td({-colspan=>1}, '');
			
			my $s = '';
			$s .=  "+$isi->{gap_left}&nbsp;bp&nbsp;::&nbsp;" if ($isi->{gap_left} > 0);
			$s .=  "&Delta;" . abs($isi->{gap_left}) . "&nbsp;bp&nbsp;::&nbsp;" if ($isi->{gap_left} < 0);
			$s .= "$isi->{family}&nbsp;(";
			$s .= (($isi->{is_strand}==+1) ? '+' : (($isi->{is_strand}==-1) ? '&minus;' : '0'));
			$s .= ")";
			$s .=  "&nbsp;::&nbsp;+$isi->{gap_right}&nbsp;bp" if ($isi->{gap_right} > 0);
			$s .=  "&nbsp;::&nbsp;&Delta;" . abs($isi->{gap_right}) . "&nbsp;bp" if ($isi->{gap_right} < 0);
			$output_str.= td({-colspan=>2, -align=>"center"}, $s);
			$output_str.= td({-colspan=>1, -align=>"center"}, '');			
			$output_str.= end_Tr;
			
			$body_countdown = 2;
		}
		
		my $key = 'side_1';			
		if ($body_countdown > 0)
		{
			$output_str.= start_Tr({-class=> "is_insertion_body"});
		}
		else
		{
			$output_str.= start_Tr({-class=> "mutation_table_row_$row_bg_color_index"});
		}

		my $redundant_annotate_key = 'junction_gene';
		$redundant_annotate_key = 'junction_is' if ($c->{$key}->{redundant} == 1);
		
		$output_str.= td({-rowspan=>2}, a({-href=>"$c->{link}"}, "*")) if ($link); 
		$output_str.= td({-rowspan=>1, -class=>"$redundant_annotate_key"}, a({-href=>"$c->{side_1}->{link}"}, "?")) if ($link); 
		$output_str.= td({-rowspan=>1, -class=>"$redundant_annotate_key"}, make_nonbreaking($c->{$key}->{seq_id}));			
		$output_str.= td( {-align=>"center", -class=>"$redundant_annotate_key"}, ($c->{$key}->{strand} == +1) ? $c->{$key}->{start} . "&nbsp;=": "=&nbsp;" . $c->{$key}->{start} );
		$output_str.= td( {-rowspan=>2, -align=>"center"}, $c->{overlap} );
		$output_str.= td( {-rowspan=>2, -align=>"center"}, "$c->{total_reads}" );
		## gene data
		{
			my $info = $c->{$key}->{$c->{$key}->{annotate_key}};				
			my $gene_string = i($info->{gene});
			$gene_string =~ s/-/&#8209;/g; #substitute nonbreaking dash		
					
			my $coord_string = $info->{interval};
			$coord_string =~ s/-/&#8209;/g; #substitute nonbreaking dash		
					
			$output_str.= td( {-align=>"center", -class=>"junction_$c->{$key}->{annotate_key}"}, $gene_string );
			$output_str.= td( {-align=>"center", -class=>"junction_$c->{$key}->{annotate_key}"}, $coord_string );
			$output_str.= td( {-class=>"junction_$c->{$key}->{annotate_key}"}, $info->{product} );
		}
		## /gene data
		
		$output_str.= end_Tr;

		$key = 'side_2';
		if ($body_countdown > 0)
		{
			$output_str.= start_Tr({-class=> "is_insertion_body"});
		}
		else
		{
			$output_str.= start_Tr({-class=> "mutation_table_row_$row_bg_color_index"});
		}
		
		$redundant_annotate_key = 'junction_gene';
		$redundant_annotate_key = 'junction_is' if ($c->{$key}->{redundant} == 1);
		
		$output_str.= td({-rowspan=>1, -class=>"$redundant_annotate_key"}, a({-href=>"$c->{side_2}->{link}"}, "?")) if ($link); 
		$output_str.= td({-rowspan=>1, -class=>"$redundant_annotate_key"}, make_nonbreaking($c->{$key}->{seq_id}));		
		$output_str.= td( {-align=>"center", -class=>"$redundant_annotate_key"}, ($c->{$key}->{strand} == +1) ? $c->{$key}->{start} . "&nbsp;=": "=&nbsp;" . $c->{$key}->{start} );
		## gene data
		{
			my $info = $c->{$key}->{$c->{$key}->{annotate_key}};				
			my $gene_string = i($info->{gene});
			$gene_string =~ s/-/&#8209;/g; #substitute nonbreaking dash		
					
			my $coord_string = $info->{interval};
			$coord_string =~ s/-/&#8209;/g; #substitute nonbreaking dash		
					
			$output_str.= td( {-align=>"center", -class=>"junction_$c->{$key}->{annotate_key}"}, $gene_string );
			$output_str.= td( {-align=>"center", -class=>"junction_$c->{$key}->{annotate_key}"}, $coord_string );
			$output_str.= td( {-class=>"junction_$c->{$key}->{annotate_key}"}, $info->{product} );
		}
		## /gene data
			
		$output_str.= end_Tr;
		
		$row_bg_color_index = ($row_bg_color_index+1)%2;
		$body_countdown--;
	}
	
	$output_str.= end_table;

}

sub html_alignment_file
{
	my ($settings, $interval) = @_;
	
	my $title = '';
	if (defined $interval->{deletion})
	{
		$title = "Deletion";
	}
	elsif (defined $interval->{hybrid})
	{
		$title = "Hybrid";
	}
	else
	{
		$title = "Mutation";
	}
	
	open HTML, ">$interval->{file_name}" or die "Could not open file: $interval->{file_name}";

	
	my $q = new CGI;

	print HTML
		start_html(
			-title => $title, 
			-head  => style({type => 'text/css'},$header_style_string),
	    );

	
	if (defined $interval->{deletion})
	{
		print HTML html_deletion_table_string([$interval->{deletion}], 0);
	}
	elsif (defined $interval->{hybrid})
	{
		print HTML html_junction_table_string($settings, [$interval->{hybrid}], undef, 0);
	}
	else
	{
		print HTML html_snp_table_string([$interval], 0);
	}
		
	print HTML p;
	my $ao = Breseq::AlignmentOutput->new;
	
	$interval->{insert_start} = 0 if (!defined $interval->{insert_start});
	$interval->{insert_end} = 0 if (!defined $interval->{insert_end});
	print HTML $ao->html_alignment($interval->{bam_path}, $interval->{fasta_path}, "$interval->{seq_id}:$interval->{start}.$interval->{insert_start}-$interval->{end}.$interval->{insert_end}", $interval);	
	
	print HTML end_html;
	close HTML;
}

sub save_text_deletion_file
{
	my ($deletion_file_name, $deletions_ref) = @_;

	open DEL, ">$deletion_file_name" or die "Could not open: $deletion_file_name";
	print DEL join("\t", 'seq_id', 'start', 'end', 'size', 'left_cov', 'left_inside_cov', 'right_inside_cov', 'right_cov', 'genes') . "\n";
	foreach my $d (@$deletions_ref)
	{
		print DEL join("\t", $d->{seq_id}, $d->{start}, $d->{end}, $d->{size}, $d->{left_unique_cov}, $d->{left_inside_unique_cov},
				$d->{right_inside_unique_cov}, $d->{right_unique_cov}, $d->{genes}) . "\n"; 
	}
	close DEL;
}

sub text_alignment_file
{
	my ($file_name, $title, $interval, $ref_seq) = @_;
	open TEXT, ">$file_name" or die "Could not open file: $file_name" or die "Could not open file: $file_name";
	
	if (defined $interval->{deletion})
	{
		### ADD deletion table at top of text alignment
	}
	elsif (defined $interval->{hybrid})
	{
		print TEXT text_junction_table_string([$interval->{hybrid}]);
	}
	else
	{
		print TEXT text_snp_table_string([$interval]);
	}
	
	##interval must have all reads loaded in preparation for producing alignment
	my $alignment = AlignmentMaker::text_alignment($interval, $ref_seq);
	
	print TEXT 	"$alignment->{aligned_ref_seq} reference\n" if (defined $alignment->{aligned_ref_seq});
	print TEXT  "$alignment->{aligned_change_seq}\n" if (defined $alignment->{aligned_change_seq});

	for my $line (@{$alignment->{lines}})
	{
		print TEXT "$line->{aligned_seq}  $line->{query}\n";
	}
	close TEXT;
}

sub read_genome_diff
{
	my ($file_name) = @_;
	
	## Create empty genome diff object.
	## Add mutations to it and then write file.
	my $gd = Breseq::GenomeDiff->new(-file_name=>$file_name);

	my $mutation_info;
	@{$mutation_info->{mutations}} = grep {$_->{type} eq 'WR'} $gd->list;
	
	foreach my $mut (@{$mutation_info->{mutations}})
	{
		$mut->{insert_start} = $mut->{insert_position};
		$mut->{insert_end} = $mut->{insert_position};
		$mut->{start} = $mut->{position};
		$mut->{end} = $mut->{position};
		$mut->{ref_seq}	= $mut->{ref_base},
		$mut->{new_seq} = $mut->{new_base},
		$mut->{total_coverage_string} = $mut->{tot_cov},
		$mut->{best_coverage_string} = $mut->{new_cov},
	}
	
	@{$mutation_info->{deletions}} = grep {$_->{type} eq 'MC'} $gd->list;	
	@{$mutation_info->{unknowns}}  = grep {$_->{type} eq 'UN'} $gd->list;
	
	foreach my $mut (@{$mutation_info->{unknowns}}, @{$mutation_info->{deletions}})
	{
		$mut->{size} = $mut->{end} - $mut->{start} + 1;
	}
	
	@{$mutation_info->{hybrids}} = grep {$_->{type} eq 'JC'} $gd->list;
	foreach my $item (@{$mutation_info->{hybrids}})
	{
		$item->{side_1}->{start} = $item->{side_1_position};
		$item->{side_1}->{end} = $item->{side_1_position};
		$item->{side_1}->{strand} = $item->{side_1_strand};
		$item->{side_1}->{seq_id} = $item->{side_1_seq_id};
		$item->{side_1}->{redundant} = $item->{side_1_redundant};
		$item->{side_1}->{overlap} = $item->{side_1_overlap};

		$item->{side_2}->{start} = $item->{side_2_position};
		$item->{side_2}->{end} = $item->{side_2_position};
		$item->{side_2}->{strand} = $item->{side_2_strand};
		$item->{side_2}->{seq_id} = $item->{side_2_seq_id};
		$item->{side_2}->{redundant} = $item->{side_2_redundant};
		$item->{side_2}->{overlap} = $item->{side_2_overlap};

		$item->{seq_id} = $item->{key};
		
		## regenerate the alignment overlap from the junction_key
		my $scj = Breseq::Shared::junction_name_split($item->{key});
		$item->{alignment_overlap} = $scj->{alignment_overlap};
		$item->{side_1}->{alignment_pos} = $scj->{side_1}->{position};
		$item->{side_2}->{alignment_pos} = $scj->{side_2}->{position};
		
		if ($item->{alignment_overlap} == 0)
		{
			$item->{start} = $item->{flanking_left};
			$item->{end} = $item->{flanking_left}+1;			
		}
		elsif ($item->{alignment_overlap} > 0)
		{
			$item->{start} = $item->{flanking_left}+1;
			$item->{end} = $item->{flanking_left}+$item->{alignment_overlap};
		}
		else ## ($item->{overlap} < 0)
		{
			$item->{start} = $item->{flanking_left}+1;
			$item->{end} = $item->{flanking_left}-$item->{alignment_overlap};
		}		
	}		
	return $mutation_info;
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

sub make_nonbreaking
{
	my ($text) = @_;
	$text =~ s/-/&#8209;/g; #substitute nonbreaking dash		
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

