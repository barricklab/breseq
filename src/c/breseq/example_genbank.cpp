// @JEB might want to use this as a template for loading Genbank files quickly...
//
// This program was prepared by the Regents of the University of California at 
// Los Alamos National Laboratory (the University) under contract No. W-7405-ENG-36 
// with the U.S. Department of Energy (DOE). All rights in the program are reserved 
// by the DOE and the University.  Permission is granted to the public to copy and 
// use this software without charge, provided that this Notice and any statement of 
// authorship are reproduced on all copies.  Neither the U.S. Government nor the 
// University makes any warranty, express or implied, or assumes any liability or 
// responsibility for the use of this software.

#include "stdafx.h"
#include "annotation.h"
#include "univout.h"

#include <sstream>

using namespace std;

// Local functions
int read_gbk_key(ifstream &fin);
bool read_locus_GBK(ifstream &fin);
void read_accession_GBK(ifstream &fin, SeqIdPtr &sip);
void read_version_GBK(ifstream &fin, SeqIdPtr &sip);
string read_source_GBK(ifstream &fin);
void read_sequence_GBK(ifstream &m_fin, DNA3Seq &m_seq, float &m_gc_content,
	float &m_a_content, float &m_t_content, float &m_g_content, float &m_c_content,
	float &m_other_content);

void read_base_count_GBK(ifstream &m_fin, DNA3Seq &m_seq, 
	float &m_gc_content, float &m_a_content, float &m_t_content, float &m_g_content, 
	float &m_c_content, float &m_other_content);

unsigned int count_bases_GBK(ifstream &m_fin, float &m_gc_content, 
	float &m_a_content, float &m_t_content, float &m_g_content, float &m_c_content,
	float &m_other_content);

void load_custom_color_gbk(GeneAnnotation &m_annot, const string &m_data);

int next_key_GBK(ifstream &m_fin, const bool &m_clear_line = true);

SeqIdPtr write_accession_GBK(SeqIdPtr &m_sip, const string &m_accession);
SeqIdPtr write_gi_GBK(SeqIdPtr &m_sip, const unsigned int &m_gi);

int parse_gene_GBK(ifstream &m_fin, GeneAnnotation &m_gene);
int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna);
int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna, GeneAnnotation &m_gene, bool &m_add_rna);
int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna);
int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna, GeneAnnotation &m_gene, bool &m_add_rna);
int parse_imp_GBK(ifstream &m_fin, GeneAnnotation &m_imp);
int parse_user_GBK(ifstream &m_fin, GeneAnnotation &m_sig);
int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds);
int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds, GeneAnnotation &m_gene, bool &m_add_gene);

int parse_field_GBK(ifstream &m_fin, pair<string, string> &m_field);

// Enumerate all possible GBK keys
enum{	GBK_EOF = 0,
		GBK_NO_KEY,
		GBK_UNKNOWN_KEY,
		GBK_LOCUS,
		GBK_ACCESSION,
		GBK_VERSION,
		GBK_SOURCE,
		GBK_FEATURES,
		GBK_ORIGIN,
		GBK_BASE_COUNT
};

// Enumerate all annotation keys
enum {
	GBK_ANNOT_END = 0,
	GBK_ANNOT_SOURCE,
	GBK_ANNOT_GENE,
	GBK_ANNOT_CDS,
	GBK_ANNOT_RNA,
	GBK_ANNOT_tRNA,
	GBK_ANNOT_IMP,
	GBK_ANNOT_USER,
	GBK_ANNOT_UNKNOWN,
	GBK_ANNOT_NONE
};

// Track the line number for error reporting!
extern unsigned long int line_number;

bool DNAMol::loadGBK(const std::string &m_filename, streampos &m_pos)
{
	// Read a Genebank flat file

	// First, open the file [in binary mode to allow use of
	// read()].
	ifstream fin(m_filename.c_str(), ios::binary);

	if(!fin){
		throw error_msg("Unable to open Genbank Flat File (gbk)");
	}

	// Are we reading from the body of this file?
	if(m_pos > 0){
		fin.seekg(m_pos);
	}
	else{
		// If we're reading from the head of the file, reset the line counter
		line_number = 1;
	}
	
	int key;
	string line;

	uout << "Parsing genbank records ..." << uout_endl;

	// Some defaults for GBK files (or until I find out how to 
	// parse these entries!)
	info_map[SOURCE] = "Unknown";

	while( (key = read_gbk_key(fin)) != GBK_EOF){
		switch(key){
			case GBK_NO_KEY:
				// Read and throw away the line
				getline(fin, line);
				line_number ++;
				break;
			case GBK_UNKNOWN_KEY:
				// Read and throw away the line
				getline(fin, line);
				line_number ++;
				break;
			case GBK_LOCUS:
				is_circular = read_locus_GBK(fin);
				break;
			case GBK_ACCESSION:
				// Load the NCBI accesion as a SeqIdPtr
				read_accession_GBK(fin, sip);
				break;
			case GBK_VERSION:
				// Load the NCBI accesion as a SeqIdPtr
				read_version_GBK(fin, sip);
				break;
			case GBK_SOURCE:
				info_map[TAXA_NAME] = read_source_GBK(fin);
				break;
			case GBK_FEATURES:
				loadGBKFeatures(fin);
				break;
			case GBK_ORIGIN:
				// The read_sequence function has the option to 
				// set the gc_content
				read_sequence_GBK(fin, seq, gc_content, a_content, t_content, 
					g_content, c_content, other_content);

				processGeneList(true /* Loading this data for the first time */);
				
				m_pos = fin.tellg();

				// All done. Is there more data to read?
				return (fin.eof() == false);
			case GBK_BASE_COUNT:
				read_base_count_GBK(fin, seq, gc_content, a_content, t_content, 
					g_content, c_content, other_content);

				break;

			default:
				throw error_msg("loadGBK: Unknown key encountered");
		};
	}

	return false;
}

void DNAMol::loadGBKFeatures(ifstream &m_fin)
{
	// Read and process all of the feature elements in a GBK file
	// Skip to the next line
	string buffer;

	getline(m_fin, buffer);
	line_number ++;
	
	int annot_key = next_key_GBK(m_fin);
	int last_annot_key = GBK_ANNOT_NONE;

	GeneAnnotation tmp_gene;

	while(annot_key != GBK_ANNOT_END){
		
		int cur_annot_key = annot_key;

		switch(annot_key){
			case GBK_ANNOT_END:
				return;
			case GBK_ANNOT_NONE:
				// Didn't find a key -- keep reading
				annot_key = next_key_GBK(m_fin);
				break;
			case GBK_ANNOT_SOURCE:
				// Skip the source feature for now. We'll need to parse
				// this feature to extract the taxon id).
				getline(m_fin, buffer);
				line_number ++;
				
				annot_key = next_key_GBK(m_fin);
				break;
			case GBK_ANNOT_GENE:
				annot_key = parse_gene_GBK(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case GBK_ANNOT_CDS:
				if(last_annot_key == GBK_ANNOT_GENE){
					bool add_gene = true;

					annot_key = parse_cds_GBK(m_fin, tmp_gene, gene_list.back(), add_gene);

					if(add_gene){
						// Save this CDS
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_cds_GBK(m_fin, tmp_gene);

					// Save this CDS
					gene_list.push_back(tmp_gene);
				}

				break;
			case GBK_ANNOT_RNA:
				if(last_annot_key == GBK_ANNOT_GENE){
					bool add_rna = true;

					annot_key = parse_rna_GBK(m_fin, tmp_gene, gene_list.back(), add_rna);

					if(add_rna){
						// Save this RNA
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_rna_GBK(m_fin, tmp_gene);

					// Save this RNA
					gene_list.push_back(tmp_gene);
				}

				break;
			case GBK_ANNOT_tRNA:
				if(last_annot_key == GBK_ANNOT_GENE){
					bool add_trna = true;

					annot_key = parse_trna_GBK(m_fin, tmp_gene, gene_list.back(), add_trna);

					if(add_trna){
						// Save this tRNA
						gene_list.push_back(tmp_gene);
					}
					// else {we copied the protein records into the last gene read}
				}
				else{
					annot_key = parse_trna_GBK(m_fin, tmp_gene);

					// Save this tRNA
					gene_list.push_back(tmp_gene);
				}

				break;
			case GBK_ANNOT_IMP:
				annot_key = parse_imp_GBK(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case GBK_ANNOT_USER:
				annot_key = parse_user_GBK(m_fin, tmp_gene);

				// Save this gene
				gene_list.push_back(tmp_gene);
				break;
			case GBK_ANNOT_UNKNOWN:
				// Do nothing for now
				break;
		};

		last_annot_key = cur_annot_key;
	}
}

int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds)
{
	SeqIdPtr sip = NULL;

	// Clear any existing info
	m_cds.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_cds.type(GeneAnnotation::CDS);

	m_cds.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if(seg_list.empty() == false){
		m_cds.segments(seg_list);
	}

	m_cds.start(range.first);
	m_cds.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 1024;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_cds.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_cds.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_cds.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			// Promote this CDS to a gene
			m_cds.type(GeneAnnotation::GENE);

			m_cds.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "EC_number"){
			m_cds.info(GeneAnnotation::EC, field.second);
		} else if(field.first == "protein_id"){
			// Set the accession
			sip = write_accession_GBK(sip, field.second);
		} else if(field.first == "db_xref"){
		
			unsigned int pos = field.second.find("GI:");

			if(pos != string::npos){
				pos += 3; /* strlen("GI:") */

				sip = write_gi_GBK(sip, 
					atoi(field.second.substr(pos, field.second.size() - pos).c_str()) );
			}
		} else if(field.first == "pseudo"){
			// Turn this record into a pseduo gene
			m_cds.type(GeneAnnotation::PSEUDO_GENE);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_cds, field.second);
		}
	}

	// Set the SeqId
	if(sip){
		m_cds.seqid(sip);

		sip = SeqIdSetFree(sip);
	}

	return annot_key;
}

int parse_cds_GBK(ifstream &m_fin, GeneAnnotation &m_cds, GeneAnnotation &m_gene, bool &m_add_gene)
{

	SeqIdPtr sip = NULL;

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	bool is_comp = read_range(m_fin, range, seg_list);

	// Does this CDS record match the last gene read?
	if((range.first == m_gene.start()) && (range.second == m_gene.stop())){
		// We have a match!
		m_add_gene = false;
	}
	else{
		// No match -- this CDS record does not correspond to the last gene read
		m_add_gene = true;

		// Clear any existing info
		m_cds.clear();
		
		// Define multiple annotation segments if needed
		if(seg_list.empty() == false){
			m_cds.segments(seg_list);
		}

		m_cds.start(range.first);
		m_cds.stop(range.second);
		m_cds.is_complement(is_comp);
		m_cds.type(GeneAnnotation::CDS);
	}

	GeneAnnotation &gene_ref = ((m_add_gene == true) ? m_cds : m_gene);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 1024;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			gene_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			gene_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			gene_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			gene_ref.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "EC_number"){
			gene_ref.info(GeneAnnotation::EC, field.second);
		} else if(field.first == "protein_id"){
			// Set the accession
			sip = write_accession_GBK(sip, field.second);
		} else if(field.first == "db_xref"){
			unsigned int pos = field.second.find("GI:");

			if(pos != string::npos){
				pos += 3; /* strlen("GI:") */

				sip = write_gi_GBK(sip, 
					atoi(field.second.substr(pos, field.second.size() - pos).c_str()) );
			}
		} else if(field.first == "pseudo"){
			// Turn this record into a pseduo gene
			gene_ref.type(GeneAnnotation::PSEUDO_GENE);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_cds, field.second);
		}
	}

	// Set the SeqId
	if(sip){
		gene_ref.seqid(sip);

		sip = SeqIdSetFree(sip);
	}

	return annot_key;
}


int parse_gene_GBK(ifstream &m_fin, GeneAnnotation &m_gene)
{
	// Clear any existing info
	m_gene.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_gene.type(GeneAnnotation::GENE);

	m_gene.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if(seg_list.empty() == false){
		m_gene.segments(seg_list);
	}

	m_gene.start(range.first);
	m_gene.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_gene.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_gene.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_gene.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_gene.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "pseudo"){
			// Turn this record into a pseudo gene
			m_gene.type(GeneAnnotation::PSEUDO_GENE);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_gene, field.second);
		}
	}

	return annot_key;
}

int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna)
{
	// Clear any existing info
	m_rna.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_rna.type(GeneAnnotation::RNA);

	m_rna.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if(seg_list.empty() == false){
		m_rna.segments(seg_list);
	}

	m_rna.start(range.first);
	m_rna.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_rna.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_rna.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_rna.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_rna.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_rna, field.second);
		}
	}

	return annot_key;
}

int parse_rna_GBK(ifstream &m_fin, GeneAnnotation &m_rna, GeneAnnotation &m_gene, bool &m_add_rna)
{
	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	bool is_comp = read_range(m_fin, range, seg_list);

	// Does this RNA record match the last gene read?
	if((range.first == m_gene.start()) && (range.second == m_gene.stop())){
		// We have a match!
		m_add_rna = false;
	}
	else{
		// No match -- this RNA record does not correspond to the last gene read
		m_add_rna = true;

		// Clear any existing info
		m_rna.clear();
		
		// Define multiple annotation segments if needed
		if(seg_list.empty() == false){
			m_rna.segments(seg_list);
		}

		m_rna.start(range.first);
		m_rna.stop(range.second);
		m_rna.is_complement(is_comp);
	}	

	GeneAnnotation &rna_ref = ((m_add_rna == true) ? m_rna : m_gene);

	rna_ref.type(GeneAnnotation::RNA);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			rna_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			rna_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			rna_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			rna_ref.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(rna_ref, field.second);
		}
	}

	return annot_key;
}

int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna)
{
	// Clear any existing info
	m_trna.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_trna.type(GeneAnnotation::tRNA);

	m_trna.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if(seg_list.empty() == false){
		m_trna.segments(seg_list);
	}

	m_trna.start(range.first);
	m_trna.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			m_trna.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_trna.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			m_trna.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_trna.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_trna, field.second);
		}
	}

	return annot_key;
}

int parse_trna_GBK(ifstream &m_fin, GeneAnnotation &m_trna, GeneAnnotation &m_gene, bool &m_add_trna)
{
	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	bool is_comp = read_range(m_fin, range, seg_list);

	// Does this RNA record match the last gene read?
	if((range.first == m_gene.start()) && (range.second == m_gene.stop())){
		// We have a match!
		m_add_trna = false;
	}
	else{
		// No match -- this tRNA record does not correspond to the last gene read
		m_add_trna = true;

		// Clear any existing info
		m_trna.clear();
		
		// Define multiple annotation segments if needed
		if(seg_list.empty() == false){
			m_trna.segments(seg_list);
		}

		m_trna.start(range.first);
		m_trna.stop(range.second);
		m_trna.is_complement(is_comp);
	}	

	GeneAnnotation &trna_ref = ((m_add_trna == true) ? m_trna : m_gene);

	trna_ref.type(GeneAnnotation::tRNA);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 32;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "gene"){
			trna_ref.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			trna_ref.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "note"){
			trna_ref.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			trna_ref.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(trna_ref, field.second);
		}
	}

	return annot_key;
}


int parse_imp_GBK(ifstream &m_fin, GeneAnnotation &m_imp)
{
	// Clear any existing info
	m_imp.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_imp.type(GeneAnnotation::IMP);

	m_imp.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if(seg_list.empty() == false){
		m_imp.segments(seg_list);
	}

	m_imp.start(range.first);
	m_imp.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 128;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "note"){
			m_imp.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "product"){
			m_imp.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "standard_name"){
			m_imp.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "db_xref"){
			m_imp.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_imp, field.second);
		}
	}

	return annot_key;
}

int parse_user_GBK(ifstream &m_fin, GeneAnnotation &m_sig)
{
	// Clear any existing info
	m_sig.clear();

	// Read the range of this annotation
	pair<unsigned int, unsigned int> range;
	list< pair<unsigned int, unsigned int> > seg_list;

	m_sig.type(GeneAnnotation::USER);

	m_sig.is_complement(read_range(m_fin, range, seg_list));

	// Define multiple annotation segments if needed
	if(seg_list.empty() == false){
		m_sig.segments(seg_list);
	}

	m_sig.start(range.first);
	m_sig.stop(range.second);

	int annot_key;

	pair<string, string> field;

	const unsigned int expected_size = 128;
	field.second.reserve(expected_size);

	while( (annot_key = parse_field_GBK(m_fin, field)) == GBK_ANNOT_NONE){
		if(field.first == "note"){
			m_sig.info(GeneAnnotation::NOTE , field.second);
		} else if(field.first == "gene"){
			m_sig.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "locus_tag"){
			m_sig.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "product"){
			m_sig.info(GeneAnnotation::PRODUCT, field.second);
		} else if(field.first == "standard_name"){
			m_sig.info(GeneAnnotation::LOCUS, field.second);
		} else if(field.first == "db_xref"){
			m_sig.info(GeneAnnotation::LOCUS_TAG, field.second);
		} else if(field.first == "custom_color"){
			
			// This is a Genomorama specific extension to the GBK standard
			load_custom_color_gbk(m_sig, field.second);
		}
	}

	return annot_key;
}

int parse_field_GBK(ifstream &m_fin, pair<string, string> &m_field)
{
	// Check for a possible annotation key
	const int annot_key = next_key_GBK(m_fin, false /* Don't clear the line */);

	if(annot_key != GBK_ANNOT_NONE){
		return annot_key;
	}

	// Read a field pair in one of two forms:
	//
	// Multiline:
	// /key="string stuff here
	//		over multiple lines"
	// Single Line
	// /key=number
	// Boolean
	// /key

	const unsigned int buffer_size = 96;
	char buffer[buffer_size];
	char *start_ptr;
	char *stop_ptr;
	char *ptr;
	
	// Count the number of matching '(' and ')'
	int paren_count = 0;

	m_fin.getline(buffer, buffer_size);
	line_number ++;
	
	// How many characters did we read?
	unsigned int len = m_fin.gcount();

	// Read the key from buffer
	start_ptr = strchr(buffer, '/');

	if(start_ptr == NULL){
		throw error_msg(":parse_field: Unable to find key start");
	}
	
	// Skip the '/' character
	start_ptr++;

	stop_ptr = strchr(start_ptr, '=');

	if(stop_ptr == NULL){
		// Note that keys do NOT have to have associated fields!
		stop_ptr = start_ptr;
		
		while((*stop_ptr != '\0') && !isspace(*stop_ptr)){
			stop_ptr++;
		}
		
		*stop_ptr = '\0';

		// Just return ...
		m_field.first = start_ptr;

		m_field.second = "";

		return annot_key;
	}

	// Save the key
	*stop_ptr = '\0';
	
	ptr = stop_ptr - 1;

	while(ptr != start_ptr){
		if(isspace(*ptr)){
			*ptr = '\0';

			ptr --;

			continue;
		}
		
		break;
	}
	
	m_field.first = start_ptr;

	// Read the field value
	start_ptr = stop_ptr + 1;

	while(isspace(*start_ptr)){
		start_ptr++;
	}

	if(*start_ptr == '('){
		paren_count ++;
	}

	// Do we have a single line field?
	if((paren_count == 0) && *start_ptr != '"'){
		// Yes -- this is a single line field
		stop_ptr = buffer + len - 2;

		while(isspace(*stop_ptr)){
			*stop_ptr = '\0';

			stop_ptr --;
		}
		
		m_field.second = start_ptr;

		return annot_key;
	}


	// We have a (possibly) multiline field.
	// skip the '"' character [but not the '(' character]
	if(paren_count == 0){
		start_ptr++;
	}

	// Clear the field variable
	m_field.second = "";

	stop_ptr = buffer + len - 2;

	while(true){
		
		while((*stop_ptr == '\r') || isspace(*stop_ptr)){
			*stop_ptr = '\0';
			stop_ptr --;
		}
		
		if(*stop_ptr == '"'){
			*stop_ptr = '\0';
			
			// all done
			m_field.second += start_ptr;

			return annot_key;
		}

		if((paren_count != 0) && (*stop_ptr == ')')){
			// Count the number of parens so far
			ptr = start_ptr + 1; // Skip the first '('

			while(*ptr != '\0'){
				if(*ptr == '('){
					paren_count++;
				}

				if(*ptr == ')'){
					paren_count--;
				}

				ptr++;
			}

			if(paren_count == 0){
				// all done
				m_field.second += start_ptr;

				return annot_key;
			}
		}

		// Save the current line
		m_field.second += start_ptr;

		// Add a space between multiple lines for easy reading
		m_field.second += ' ';

		// Read another line
		m_fin.getline(buffer, buffer_size);
		line_number ++;
		
		len = m_fin.gcount();
		
		// Empty lines are not allowed!
		if(len == 0){
			throw error_msg("Unexpected blank line or end of file encountered");
		}

		start_ptr = buffer;

		while(isspace(*start_ptr)){
			start_ptr++;
		}

		stop_ptr = buffer + len - 2;
	}

	return annot_key;
}

int next_key_GBK(ifstream &m_fin, const bool &m_clear_line /* = true */)
{
	const int buffer_size = 21;
	
	char buffer[buffer_size + 1];
	char *start_ptr = buffer;

	// Terminate the array
	buffer[buffer_size] = '\0';

	if(m_fin.read(buffer, buffer_size) == false){
		throw error_msg("Unable to read next annotation key");
	}
	
	// Find the start of the string
	while(isspace(*start_ptr)){
		start_ptr++;
	}

	// Is this an empty string?
	if(*start_ptr == '\0'){
		if(m_clear_line){
			// Throw away the rest of the line
			m_fin.ignore(80, '\n');
			line_number ++;
		}

		return GBK_ANNOT_NONE;
	}

	// Make the sting upper case
	for(char *ptr = start_ptr;*ptr != '\0';ptr++){
		*ptr = toupper(*ptr);
	}

	// Have we read into the base count section?
	if(strncmp(start_ptr, "BASE", 4 /*strlen("BASE")*/) == 0){
		// Rewind the stream by buffer_size characters
		m_fin.seekg(-buffer_size, ios::cur);

		return GBK_ANNOT_END;
	}

	if(strncmp(start_ptr, "ORIGIN", 6 /*strlen("ORIGIN")*/) == 0){
		// Rewind the stream by buffer_size characters
		m_fin.seekg(-buffer_size, ios::cur);

		return GBK_ANNOT_END;
	}

	if(strncmp(start_ptr, "CDS", 3 /*strlen("CDS")*/) == 0){
		return GBK_ANNOT_CDS;
	}

	if(strncmp(start_ptr, "SOURCE", 6 /*strlen("SOURCE")*/) == 0){
		return GBK_ANNOT_SOURCE;
	}

	if(strncmp(start_ptr, "GENE", 4 /*strlen("GENE")*/) == 0){
		return GBK_ANNOT_GENE;
	}

	if(strncmp(start_ptr, "TRNA", 4 /*strlen("TRNA")*/) == 0){
		return GBK_ANNOT_tRNA;
	}

	// Does the buffer contain the string "RNA"?
	if(strstr(start_ptr, "RNA") != NULL){
		return GBK_ANNOT_RNA;
	}
	
	if(strncmp(start_ptr, "USER", 4 /*strlen("USER")*/) == 0){
		return GBK_ANNOT_USER;
	}
	
	// Did not match this key
	return GBK_ANNOT_IMP;
}

void read_base_count_GBK(ifstream &m_fin, DNA3Seq &m_seq,
						 float &m_gc_content, float &m_a_content, float &m_t_content, float &m_g_content, 
						 float &m_c_content, float &m_other_content)
{
	// Read the number of bases to expect in the sequence
	const char buffer_size = 126;
	char buffer[buffer_size];
	char *start_ptr, *stop_ptr;

	m_fin.getline(buffer, buffer_size);
	line_number ++;
	
	///////////// Read the number of a's /////////////
	// 1) skip leading spaces
	start_ptr = buffer;

	while(isdigit(*start_ptr) == false){
		start_ptr++;
	}

	// 2) Find the end of the "a" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int a_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;

	///////////// Read the number of c's /////////////
	// 1) skip leading spaces
	while(isdigit(*start_ptr) == false){
		start_ptr++;
	}

	// 2) Find the end of the "c" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int c_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;

	//////////// Read the number of g's /////////////
	// 1) skip leading spaces
	while(isdigit(*start_ptr) == false){
		start_ptr++;
	}

	// 2) Find the end of the "g" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int g_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;
	
	//////////// Read the number of t's /////////////
	// 1) skip leading spaces
	while(isdigit(*start_ptr) == false){
		start_ptr++;
	}

	// 2) Find the end of the "t" counts
	for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
		// Do nothing!
	}
	
	*stop_ptr = '\0';

	const int t_count = atoi(start_ptr);

	start_ptr = stop_ptr + 1;

	//////////// Read the number of Others's /////////////
	// 1) skip leading spaces but allow for optional "others" count
	while((*start_ptr != '\0') && !isdigit(*start_ptr)){
		start_ptr++;
	}
	
	int other_count = 0;

	if(*start_ptr != '\0'){
		// 2) Find the end of the "others" counts
		for(stop_ptr = start_ptr;isdigit(*stop_ptr);stop_ptr++){
			// Do nothing!
		}
		
		*stop_ptr = '\0';

		other_count = atoi(start_ptr);
	}

	// Set the sequence size
	const int num_base = a_count + t_count + g_count + c_count + other_count;

	m_gc_content = ((float)(g_count + c_count))/(num_base);
	m_a_content = ((float)(a_count))/(num_base);
	m_t_content = ((float)(t_count))/(num_base);
	m_g_content = ((float)(g_count))/(num_base);
	m_c_content = ((float)(c_count))/(num_base);
	m_other_content = ((float)(other_count))/(num_base);

	m_seq = DNA3Seq(num_base);
}

// Read nucleotide sequence and return the gc content
void read_sequence_GBK(ifstream &m_fin, DNA3Seq &m_seq, float &m_gc_content,
					   float &m_a_content, float &m_t_content, float &m_g_content, float &m_c_content,
					   float &m_other_content)
{
	streampos pos = m_fin.tellg();

	// Have we read the number of bases in the sequence?
	if(m_seq.size() == 0){

		// We don't know the size of the DNA molecule
		// Count the number of bases and ALSO the GC content!
		unsigned int mol_size = count_bases_GBK(m_fin, m_gc_content,
			m_a_content, m_t_content, m_g_content, m_c_content,
			m_other_content);

		m_seq = DNA3Seq(mol_size);
	
		// Restore the state of the input file
		m_fin.clear();

		m_fin.seekg(pos);
	}

	// The buffer size is a parameter that must be tuned for the file IO of 
	// a given machine.
	const unsigned int buffer_size = 2046;
	char buffer[buffer_size];
	char *ptr;
	
	unsigned int base_count = 0;
	DNA3Seq::iterator iter = m_seq.begin();
	
	// First, throw away the line that contains "ORIGIN"
	m_fin.getline(buffer, buffer_size);
	line_number ++;
	
	// Track the current position
	pos += strlen(buffer);

	// Read until we hit a "/" symbol or reach the end of the file
	while(m_fin.eof() == false){
		m_fin.read(buffer, buffer_size);
		ptr = buffer;
		
		const unsigned int len = m_fin.gcount();

		for(unsigned int i = 0;i < len;i++, ptr++){
			
			// We expect upper case letters, so test for these first
			if( (*ptr >= 'A') && (*ptr <= 'Z') ){
				// Save this base
				if(*ptr == 'A'){
					iter.base(DNA3Seq::A);
				}
				else{
					if(*ptr == 'T'){
						iter.base(DNA3Seq::T);
					}
					else{
						if(*ptr == 'G'){
							iter.base(DNA3Seq::G);
						}
						else{
							if(*ptr == 'C'){
								iter.base(DNA3Seq::C);
							}
							else{
								iter.base(DNA3Seq::UNKNOWN);
							}
						}
					}
				}
								
				iter++;

				// Keep track of the number of bases read
				base_count++;
				
				continue;
			}
			
			// Is this a lower case letter?
			if( (*ptr >= 'a') && (*ptr <= 'z') ){
				// Save this base
				if(*ptr == 'a'){
					iter.base(DNA3Seq::A);
				}
				else{
					if(*ptr == 't'){
						iter.base(DNA3Seq::T);
					}
					else{
						if(*ptr == 'g'){
							iter.base(DNA3Seq::G);
						}
						else{
							if(*ptr == 'c'){
								iter.base(DNA3Seq::C);
							}
							else{
								iter.base(DNA3Seq::UNKNOWN);
							}
						}
					}
				}
								
				iter++;

				// Keep track of the number of bases read
				base_count++;
				
				continue;
			}
			
			if(*ptr == '/'){
				// all done!
				if(base_count != m_seq.size()){
					throw error_msg("Did not read enough bases");
				}
				
				pos += i;

				m_fin.clear();
				m_fin.seekg(pos);

				return;
			}
		}

		// Update the current position
		pos += len;
	}

	throw error_msg(":read_sequence: Could not find end-of-sequence terminator");
}

unsigned int count_bases_GBK(ifstream &m_fin, float &m_gc_content,
							 float &m_a_content, float &m_t_content, float &m_g_content, float &m_c_content,
							 float &m_other_content)
{
	// The buffer size is a parameter that must be tuned for the file IO of 
	// a given machine.
	const unsigned int buffer_size = 2048;
	char buffer[buffer_size];
	char *ptr;
	
	unsigned int base_count = 0;
	
	// Zero the GC content
	m_gc_content = 0.0f;
	m_a_content = 0.0f;
	m_t_content = 0.0f;
	m_g_content = 0.0f;
	m_c_content = 0.0f;
	m_other_content = 0.0f;

	// First, throw away the line that contains "ORIGIN"
	m_fin.getline(buffer, buffer_size);
	
	// We're going to rewind, so don't increment the line number counter
	
	// Read until we hit a "/" symbol or reach the end of the file
	while(m_fin.eof() == false){
		m_fin.read(buffer, buffer_size);
		ptr = buffer;
		
		const unsigned int len = m_fin.gcount();

		for(unsigned int i = 0;i < len;i++, ptr++){

			// We expect upper case letters, so test for these first
			if( (*ptr >= 'A') && (*ptr <= 'Z') ){
				if( *ptr == 'A' ){
					m_a_content++;
					base_count ++;
					continue;
				}

				if( *ptr == 'T' ){
					m_t_content++;
					base_count ++;
					continue;
				}
				
				if( *ptr == 'G' ){
					m_gc_content ++;
					m_g_content++;
					base_count ++;
					continue;
				}
				
				if( *ptr == 'C' ){
					m_gc_content ++;
					m_c_content++;
					base_count ++;
					continue;
				}
				
				
				m_other_content++;
				base_count++;
				continue;
			}
			
			if(*ptr == '/'){
				// all done!
				if(base_count > 0){
					m_gc_content /= base_count;
					m_a_content /= base_count;
					m_t_content /= base_count;
					m_g_content /= base_count;
					m_c_content /= base_count;
					m_other_content /= base_count;
				}

				return base_count;
			}
			
			// Test for lower case letters
			if( (*ptr >= 'a') && (*ptr <= 'z') ){
				if( *ptr == 'a' ){
					m_a_content++;
					base_count ++;
					continue;
				}

				if( *ptr == 't' ){
					m_t_content++;
					base_count ++;
					continue;
				}
				
				if( *ptr == 'g' ){
					m_gc_content ++;
					m_g_content++;
					base_count ++;
					continue;
				}
				
				if( *ptr == 'c' ){
					m_gc_content ++;
					m_c_content++;
					base_count ++;
					continue;
				}
				
				// Is this a letter?
				m_other_content++;
				base_count++;
				continue;
			}
					
			#ifdef OPTIMIZE
			 
			switch(buffer[i]){
				case '0':
				case '1':
				case '2':
				case '3':
				case '4':
				case '5':
				case '6':
				case '7':
				case '8':
				case '9':
				case ' ':
				case '\t':
				case '\n':
				case '\r':
					// Skip digits and white space
					break;
				case '/':
					// all done!

					if(base_count > 0){
						m_gc_content /= base_count;
						m_a_content /= base_count;
						m_t_content /= base_count;
						m_g_content /= base_count;
						m_c_content /= base_count;
						m_other_content /= base_count;
					}

					return base_count;
				case 'a':
				case 'A':
					m_a_content++;
					base_count ++;
					break;
				case 't':
				case 'T':
					m_t_content++;
					base_count ++;
					break;
				case 'g':
				case 'G':
					m_g_content++;
					m_gc_content++;
					base_count ++;
					break;
				case 'c':
				case 'C':
					m_c_content++;
					m_gc_content++;
					base_count ++;
					break;
				default:
					m_other_content++;
					base_count++;
					break;
			};
			
			#endif // OPTIMIZE
		}
	}

	throw error_msg(":count_bases: Could not find end-of-sequence terminator");

	return 0;
}

string read_source_GBK(ifstream &fin)
{
	string taxa;
	string buffer;

	getline(fin, buffer);
	line_number ++;

	stringstream ss(buffer);

	while(ss >> buffer){
		if(taxa.empty()){
			taxa = buffer;
		}
		else{
			taxa += " " + buffer;
		}
	}

	return taxa;
}

void read_version_GBK(ifstream &fin, SeqIdPtr &sip)
{
	string buffer;

	getline(fin, buffer);
	line_number ++;

	stringstream ss(buffer);

	buffer = "";

	ss >> buffer;

	if(buffer.empty()){
		// No version info to read
		return;
	}

	// Is this a gi or accession?
	unsigned int pos = buffer.find(":");

	if(pos == string::npos){
		// Accession
		sip = write_accession_GBK(sip, buffer);
	}
	else{
		// Gi
		sip = write_gi_GBK(sip, atoi(buffer.c_str() + pos + 1) );
	}

	buffer = "";

	ss >> buffer;

	if(buffer.empty()){
		// No version info to read
		return;
	}

	// Is this a gi or accession?
	pos = buffer.find(":");

	if(pos == string::npos){
		// Accession
		sip = write_accession_GBK(sip, buffer);
	}
	else{
		// Gi
		sip = write_gi_GBK(sip, atoi(buffer.c_str() + pos + 1) );
	}
}

void read_accession_GBK(ifstream &fin, SeqIdPtr &sip)
{
	string accession;

	fin >> accession;

	// Convert this string into a valid SeqIdPtr.
	// Overwrite any existing accession entries.
	sip = write_accession_GBK(sip, accession);
}

// Write an accession to the given SeqIdPtr. If an accession entry
// already exists, this code will overwrite it! The updated SedIdPtr is
// returned and the input pointer is invalid.
SeqIdPtr write_accession_GBK(SeqIdPtr &m_sip, const string &m_accession)
{
	SeqIdPtr sip = NULL;
	SeqIdPtr tmp_new = NULL;
	SeqIdPtr tmp_old = NULL;

	// Is this a properly formatted accession?
	if(m_accession.find('|') != string::npos){
		// Yes
		sip = SeqIdParse ((char*)m_accession.c_str());
	}
	else{
		// No
		sip = SeqIdParse( (char*)string("gb|" + m_accession).c_str() );
	}

	if(m_sip == NULL){
		return sip;
	}

	if(sip == NULL){
		throw error_msg(":write_accession: Unable to parse SeqId");
	}

	tmp_new = sip;
	tmp_old = m_sip;

	while(tmp_old != NULL){
		if(tmp_old->choice != SEQID_GENBANK){
			tmp_new->next = SeqIdDup(tmp_old);
			tmp_new = tmp_new->next;
		}

		tmp_old = tmp_old->next;
	}

	// Free the old ptr
	m_sip = SeqIdSetFree(m_sip);

	return sip;
}

// Write a gi to the given SeqIdPtr. If a gi entry
// already exists, this code will overwrite it! The updated SedIdPtr is
// returned and the input pointer is invalid.
SeqIdPtr write_gi_GBK(SeqIdPtr &m_sip, const unsigned int &m_gi)
{
	SeqIdPtr sip = NULL;
	SeqIdPtr tmp_new = NULL;
	SeqIdPtr tmp_old = NULL;

	sip = ValNodeNew(NULL);
	sip->choice = SEQID_GI;
	sip->data.intvalue = m_gi;

	if(m_sip == NULL){
		return sip;
	}

	tmp_new = sip;
	tmp_old = m_sip;

	while(tmp_old != NULL){
		if(tmp_old->choice != SEQID_GI){
			tmp_new->next = SeqIdDup(tmp_old);
			tmp_new = tmp_new->next;
		}

		tmp_old = tmp_old->next;
	}

	// Free the old ptr
	m_sip = SeqIdSetFree(m_sip);

	return sip;
}

bool read_locus_GBK(ifstream &fin)
{
	string line;

	getline(fin, line);
	line_number ++;

	string::iterator iter;

	for(iter = line.begin();iter != line.end(); iter++){
		*iter = tolower(*iter);
	}

	return (string::npos != line.find("circular"));
}

int read_gbk_key(ifstream &fin)
{
	const int buffer_size = 12;
	
	char buffer[buffer_size + 1];

	// Terminate the array
	buffer[buffer_size] = '\0';

	if(fin.read(buffer, buffer_size) == false){
		return GBK_EOF;
	}

	stringstream ssin(buffer);
	string key;

	ssin >> key;

	if(key.empty()){
		return GBK_NO_KEY;
	}

	string::iterator iter;
	

	for(iter = key.begin();iter != key.end();iter++){
		*iter = toupper(*iter);
	}

	if(key == "LOCUS"){
		return GBK_LOCUS;
	}

	if(key == "ACCESSION"){
		return GBK_ACCESSION;
	}

	if(key == "VERSION"){
		return GBK_VERSION;
	}

	if(key == "SOURCE"){
		return GBK_SOURCE;
	}

	if(key == "FEATURES"){
		return GBK_FEATURES;
	}

	if(key == "ORIGIN"){
		return GBK_ORIGIN;
	}
	
	if(key == "BASE"){
		return GBK_BASE_COUNT;
	}

	// Did not match this key
	return GBK_UNKNOWN_KEY;
}

void load_custom_color_gbk(GeneAnnotation &m_annot, const string &m_data)
{
	// Parse the /custom_color="FORMAT:XX, YY, ZZ, ..."
	// Where the only allowed format is currently RGBA and each of the four colors
	// is an integer between [0, 255]
	 
	unsigned int pos = m_data.find("RGBA:");

	if(pos != string::npos){
	
		pos += 5; /* strlen("RGBA:") */
		
		const unsigned int len = m_data.size();
		
		// Skip any leading spaces
		while( (pos < len) && !isdigit(m_data[pos]) ){
			pos++;
		}
		
		// Read the red channel
		unsigned int start = pos;
		
		while( (pos < len) && isdigit(m_data[pos]) ){
			pos++;
		}
		
		const int red = atoi( m_data.substr(start, pos - start).c_str() );
		
		// Skip any leading spaces and comma
		while( (pos < len) && !isdigit(m_data[pos]) ){
			pos++;
		}
		
		start = pos;
		
		// Read the green channel
		while( (pos < len) && isdigit(m_data[pos]) ){
			pos++;
		}
		
		const int green = atoi( m_data.substr(start, pos - start).c_str() );
		
		// Skip any leading spaces and comma
		while( (pos < len) && !isdigit(m_data[pos]) ){
			pos++;
		}
		
		start = pos;
		
		// Read the blue channel
		while( (pos < len) && isdigit(m_data[pos]) ){
			pos++;
		}
		
		const int blue = atoi( m_data.substr(start, pos - start).c_str() );
		
		// Skip any leading spaces and comma
		while( (pos < len) && !isdigit(m_data[pos]) ){
			pos++;
		}
		
		start = pos;
		
		// Read the alpha channel
		while( (pos < len) && isdigit(m_data[pos]) ){
			pos++;
		}
		
		const int alpha = atoi( m_data.substr(start, pos - start).c_str() );
		
		// Validate the color values
		if( (red < 0) || (red > 255) ){
			throw ":load_custom_color_gbk: red channel out of bounds";
		}
		
		if( (green < 0) || (green > 255) ){
			throw ":load_custom_color_gbk: green channel out of bounds";
		}
		
		if( (blue < 0) || (blue > 255) ){
			throw ":load_custom_color_gbk: blue channel out of bounds";
		}
		
		if( (alpha < 0) || (alpha > 255) ){
			throw ":load_custom_color_gbk: alpha channel out of bounds";
		}
		
		float color_norm = 1.0f/255.0f;
		
		m_annot.color(red*color_norm, green*color_norm, blue*color_norm);
		m_annot.color_alpha(alpha*color_norm);
	}

}
