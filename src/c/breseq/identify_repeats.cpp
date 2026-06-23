/*****************************************************************************

AUTHORS

  Jeffrey E. Barrick <jeffrey.e.barrick@gmail.com>
  David B. Knoester

LICENSE AND COPYRIGHT

  Copyright (c) 2008-2010 Michigan State University
  Copyright (c) 2011-2022 The University of Texas at Austin

  breseq is free software; you can redistribute it and/or modify it under the
  terms the GNU General Public License as published by the Free Software
  Foundation; either version 1, or (at your option) any later version.

*****************************************************************************/

#include "libbreseq/identify_repeats.h"

#include "libbreseq/reference_sequence.h"

using namespace std;

namespace breseq {

// Splits one CSV line on commas, respecting double-quoted fields (which may
// not themselves contain embedded double-quotes -- write_repeats_csv never
// produces those, so we don't need to handle quote-escaping here).
static vector<string> split_csv_line(const string& line)
{
  vector<string> fields;
  string field;
  bool in_quotes = false;

  for (size_t i = 0; i < line.size(); i++) {
    char c = line[i];
    if (c == '"') {
      in_quotes = !in_quotes;
    } else if ((c == ',') && !in_quotes) {
      fields.push_back(field);
      field.clear();
    } else {
      field += c;
    }
  }
  fields.push_back(field);

  return fields;
}

vector<cRepeatMatch> parse_show_coords_file(const string& coords_file_name)
{
  vector<cRepeatMatch> result;

  ifstream in(coords_file_name.c_str());
  assert(in);

  string line;
  while (getline(in, line)) {
    if (line.size() == 0) continue;

    vector<string> f = split_on_whitespace(line);
    ASSERT(f.size() == 11, "Expected 11 columns in show-coords output (run with flags '-r -T -l -H'), found "
           + to_string(f.size()) + " on line:\n" + line);

    cRepeatMatch m;
    m.start_1 = from_string<uint32_t>(f[0]);
    m.end_1   = from_string<uint32_t>(f[1]);

    // The query side (S2/E2) is reported with start>end when it matched the
    // reverse complement of the reference side -- normalize that into strand_2.
    uint32_t raw_start_2 = from_string<uint32_t>(f[2]);
    uint32_t raw_end_2   = from_string<uint32_t>(f[3]);
    if (raw_start_2 <= raw_end_2) {
      m.strand_2 = 1;
      m.start_2  = raw_start_2;
      m.end_2    = raw_end_2;
    } else {
      m.strand_2 = -1;
      m.start_2  = raw_end_2;
      m.end_2    = raw_start_2;
    }

    m.length_1 = from_string<uint32_t>(f[4]);
    m.length_2 = from_string<uint32_t>(f[5]);
    m.percent_identity = from_string<double>(f[6]);
    // f[7]/f[8] are the total reference/query sequence lengths -- not needed here.
    m.seq_id_1 = f[9];
    m.seq_id_2 = f[10];

    result.push_back(m);
  }

  return result;
}

// Identifies a single region for graph-building purposes in group_repeats().
struct cRegionKey {
  string   seq_id;
  uint32_t start;
  uint32_t end;

  bool operator<(const cRegionKey& other) const {
    if (seq_id != other.seq_id) return seq_id < other.seq_id;
    if (start != other.start) return start < other.start;
    return end < other.end;
  }
};

vector<cRepeatRegion> group_repeats(const vector<cRepeatMatch>& matches)
{
  // Build the undirected match graph: one node per distinct region, edges
  // connecting the two regions of each match.
  map<cRegionKey, set<cRegionKey> > adjacency;
  map<cRegionKey, uint32_t> region_length;

  for (vector<cRepeatMatch>::const_iterator it = matches.begin(); it != matches.end(); it++) {
    const cRepeatMatch& m = *it;
    cRegionKey k1; k1.seq_id = m.seq_id_1; k1.start = m.start_1; k1.end = m.end_1;
    cRegionKey k2; k2.seq_id = m.seq_id_2; k2.start = m.start_2; k2.end = m.end_2;

    adjacency[k1].insert(k2);
    adjacency[k2].insert(k1);
    region_length[k1] = m.length_1;
    region_length[k2] = m.length_2;
  }

  // map<> already iterates its keys in sorted (seq_id, start, end) order, so this
  // gives deterministic node visitation order and, in turn, deterministic group IDs.
  vector<cRegionKey> all_nodes;
  for (map<cRegionKey, set<cRegionKey> >::const_iterator it = adjacency.begin(); it != adjacency.end(); it++)
    all_nodes.push_back(it->first);

  vector<cRepeatRegion> result;
  set<cRegionKey> visited;
  uint32_t next_group_id = 1;

  for (size_t i = 0; i < all_nodes.size(); i++) {
    if (visited.count(all_nodes[i])) continue;

    // Breadth-first traversal to find this node's connected component.
    vector<cRegionKey> component;
    vector<cRegionKey> to_visit;
    to_visit.push_back(all_nodes[i]);
    visited.insert(all_nodes[i]);
    while (!to_visit.empty()) {
      cRegionKey node = to_visit.back();
      to_visit.pop_back();
      component.push_back(node);

      const set<cRegionKey>& neighbors = adjacency[node];
      for (set<cRegionKey>::const_iterator nit = neighbors.begin(); nit != neighbors.end(); nit++) {
        if (!visited.count(*nit)) {
          visited.insert(*nit);
          to_visit.push_back(*nit);
        }
      }
    }
    sort(component.begin(), component.end());

    // A group requires every pair within the component to be a match (a clique).
    bool is_clique = true;
    for (size_t a = 0; is_clique && (a < component.size()); a++) {
      const set<cRegionKey>& neighbors = adjacency[component[a]];
      for (size_t b = 0; b < component.size(); b++) {
        if (a == b) continue;
        if (!neighbors.count(component[b])) { is_clique = false; break; }
      }
    }

    if (is_clique) {
      double min_identity = 100.0;
      for (vector<cRepeatMatch>::const_iterator it = matches.begin(); it != matches.end(); it++) {
        cRegionKey k1; k1.seq_id = it->seq_id_1; k1.start = it->start_1; k1.end = it->end_1;
        cRegionKey k2; k2.seq_id = it->seq_id_2; k2.start = it->start_2; k2.end = it->end_2;
        if (binary_search(component.begin(), component.end(), k1) && binary_search(component.begin(), component.end(), k2))
          min_identity = min(min_identity, it->percent_identity);
      }

      uint32_t group_id = next_group_id++;
      for (size_t c = 0; c < component.size(); c++) {
        cRepeatRegion r;
        r.group_id = group_id;
        r.group_size = component.size();
        r.seq_id = component[c].seq_id;
        r.start = component[c].start;
        r.end = component[c].end;
        r.length = region_length[component[c]];
        r.min_percent_identity = min_identity;
        result.push_back(r);
      }
    } else {
      // The connected component isn't a clean clique (not every pair matched
      // above threshold) -- leave each region as its own ungrouped singleton
      // rather than guess at sub-groups. See group_repeats()'s header comment.
      cerr << "  NOTE: " << component.size() << " regions near " << component[0].seq_id << ":" << component[0].start
           << " form a connected cluster of repeats, but not every pair matched above the threshold." << endl
           << "        They were left ungrouped (each assigned its own group ID)." << endl;

      for (size_t c = 0; c < component.size(); c++) {
        cRepeatRegion r;
        r.group_id = next_group_id++;
        r.group_size = 1;
        r.seq_id = component[c].seq_id;
        r.start = component[c].start;
        r.end = component[c].end;
        r.length = region_length[component[c]];
        r.min_percent_identity = 0;
        result.push_back(r);
      }
    }
  }

  return result;
}

void write_repeats_csv(const vector<cRepeatRegion>& regions, const string& file_name)
{
  string path = path_to_dirname(file_name);
  create_path(path);

  ofstream out(file_name.c_str());
  out << join(make_vector<string>
              ("group_id")("group_size")
              ("seq_id")("start")("end")("length")
              ("min_percent_identity"),
              ",") << endl;

  for (vector<cRepeatRegion>::const_iterator it = regions.begin(); it != regions.end(); it++) {
    const cRepeatRegion& r = *it;
    out << join(make_vector<string>
                (to_string(r.group_id))
                (to_string(r.group_size))
                (double_quote(r.seq_id))
                (to_string(r.start))
                (to_string(r.end))
                (to_string(r.length))
                (r.group_size > 1 ? to_string(r.min_percent_identity, 2) : string("")),
                ",") << endl;
  }
}

vector<cRepeatRegion> read_repeats_csv(const string& file_name)
{
  vector<cRepeatRegion> regions;

  ifstream in(file_name.c_str());
  assert(in);

  string line;
  bool header_skipped = false;
  while (getline(in, line)) {
    if (line.size() == 0) continue;
    if (!header_skipped) {
      header_skipped = true;
      continue;
    }

    vector<string> f = split_csv_line(line);
    ASSERT(f.size() == 7, "Expected 7 columns in repeats CSV file, found "
           + to_string(f.size()) + " on line:\n" + line);

    cRepeatRegion r;
    r.group_id   = from_string<uint32_t>(f[0]);
    r.group_size = from_string<uint32_t>(f[1]);
    r.seq_id     = f[2];
    r.start      = from_string<uint32_t>(f[3]);
    r.end        = from_string<uint32_t>(f[4]);
    r.length     = from_string<uint32_t>(f[5]);
    r.min_percent_identity = f[6].empty() ? 0 : from_string<double>(f[6]);

    regions.push_back(r);
  }

  return regions;
}

void identify_repeats(
                      const vector<string>& reference_file_names,
                      const string& output_file_name,
                      const uint32_t minimum_length,
                      const double minimum_identity
                      )
{
  ASSERT(SYSTEM_CAPTURE("which nucmer", true).size() > 0,
         "Required executable \"nucmer\" not found in PATH.\nInstall MUMmer, e.g. 'conda install -c bioconda mummer4'.");
  ASSERT(SYSTEM_CAPTURE("which show-coords", true).size() > 0,
         "Required executable \"show-coords\" not found in PATH.\nInstall MUMmer, e.g. 'conda install -c bioconda mummer4'.");

  cReferenceSequences ref_seq_info;
  ref_seq_info.LoadFiles(reference_file_names);

  // Working files live next to the output CSV.
  string out_dir = path_to_dirname(output_file_name);
  create_path(out_dir);
  string stem        = out_dir + "/" + path_to_filename(output_file_name) + ".identify_repeats";
  string fasta_path   = stem + ".fna";
  string delta_path   = stem + ".delta";
  string coords_path  = stem + ".coords";

  ref_seq_info.WriteFASTA(fasta_path);

  // Self-vs-self alignment. --maxmatch (rather than the default --mum/--mumreference)
  // is required to recover every pairwise combination among repeat families with
  // three or more copies, since the default modes only report matches that are
  // unique in the reference and/or query.
  SYSTEM("nucmer --maxmatch -p " + double_quote(stem) + " " + double_quote(fasta_path) + " " + double_quote(fasta_path),
         false, false, false);

  SYSTEM("show-coords -r -T -l -H " + double_quote(delta_path) + " > " + double_quote(coords_path),
         false, false, false);

  vector<cRepeatMatch> raw_matches = parse_show_coords_file(coords_path);

  vector<cRepeatMatch> filtered;
  for (vector<cRepeatMatch>::iterator it = raw_matches.begin(); it != raw_matches.end(); it++) {
    cRepeatMatch& m = *it;

    // Skip the trivial match of an entire region against itself (e.g. the whole-
    // sequence self-diagonal).
    if ((m.seq_id_1 == m.seq_id_2) && (m.start_1 == m.start_2) && (m.end_1 == m.end_2))
      continue;

    if ((m.length_1 < minimum_length) || (m.length_2 < minimum_length))
      continue;

    if (m.percent_identity < minimum_identity)
      continue;

    // A self-vs-self alignment reports every off-diagonal repeat pair twice (once
    // with each side as the "reference"). Keep only one canonical ordering.
    if (make_pair(m.seq_id_1, make_pair(m.start_1, m.end_1))
        > make_pair(m.seq_id_2, make_pair(m.start_2, m.end_2)))
      continue;

    filtered.push_back(m);
  }

  vector<cRepeatRegion> grouped = group_repeats(filtered);
  write_repeats_csv(grouped, output_file_name);

  remove_file(fasta_path);
  remove_file(delta_path);
  remove_file(coords_path);
}

}
