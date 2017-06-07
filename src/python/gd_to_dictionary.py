#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Thu June  1 16:56:12 2017
read in all contents of all specified .gd file into python dictionary


improvements needed:

@author: ded
"""

import sys
import re
import argparse
from collections import OrderedDict
import os


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""
    # Taken from https://gist.github.com/brantfaircloth/1252339 to allow relative path in argparse for gd directory identification
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

parser = argparse.ArgumentParser(description="read in gd file names, generate dictionary")
parser.add_argument("-g", "--gd", help="direcotry contianing gd files to compare", action=FullPaths)
parser.add_argument('gdnames', nargs='*', help="filenames of .gd files to compare. CAN NOT BE USED WITH -g option")
args = parser.parse_args()

# assertions to make sure correct options or option sets are selected, and any relevent warnings
assert not (args.gd is None and args.gdnames == []), "gd files or directory containing gd files must be specified"
assert args.gd is None or args.gdnames == [], "gd files or directoy must be specified, not both"

# generate files_to_compare list of gd files to be compared regardless of how they were entered.
if args.gd:
    os.chdir(args.gd)
    files_to_compare = [x for x in os.listdir(".") if x.endswith(".gd")]
else:
    files_to_compare = [x for x in args.gdnames if x.endswith(".gd")]
    assert len(files_to_compare) == len(args.gdnames), "Non Genome Diff file provided: %s" %([x for x in args.gdnames if not x.endswith(".gd")])


master_gd_dict = {}
named_fields = []

# generate master dictionary
for files in files_to_compare:
    sample = files.split("/")[-1].replace(".gd", "")  # remove path information from sample name
    assert files.endswith(".gd"), "non-gd file found. This should never trigger as previous code should have eliminated this file: %s" % files

    assert sample not in master_gd_dict, "Duplicated gd name. %s read in despite already being in dictionary from one of the other .gds listed:\n%s." % (sample, files_to_compare)
    master_gd_dict[sample] = OrderedDict()  # ordered to ensure headers stored at front

    with open(files, "r") as F:
        for line in F:
            line = line.rstrip().split('\t')
            commented_out = False

            # Header Lines  TODO: decide if this needs to be stored as a dictionary of items, instead of single list of 'header'
            if re.match("^#=", line[0]):
                assert len(line) <= 2, "The following line begins with header notation, but contains %i fields rather than the expected 2 (or 1 in the case of version line) meaning it is not formatted correctly.\n%s" % (len(line), line)

                if len(line) == 1:  # Version line
                    assert line[0] == "#=GENOME_DIFF 1.0", "The following line lacks tabs, and is not the #=GENOME_DIFF 1.0 line:\n%s" % line
                    assert "header" not in master_gd_dict[sample], "Multiple version lines detected, possible formatting error:\n%s" % line
                    master_gd_dict[sample]["header"] = [line]  # Defined here to ensure is first line, else block will throw KeyError if version line not detected previously

                try:  # All other header/metadata lines
                    master_gd_dict[sample]["header"].append(line)
                except KeyError:
                    assert False, "The following line was encountered before the version line. Please ensure file starts with #=GENOME_DIFF 1.0. Note space, NOT tab between F and 1.\n%s" % line

                continue  # move to next line to avoid final block regarding extra key=value fields

            # Identify lines that have been commented out
            if re.match("^#", line[0]):
                line[0] = line[0].replace("#", "")  # Remove comment designation so line can be stored using existing lines, then reapply later
                commented_out = True  # Needed to replace comment status after storage

            # Verify ID is unique and store 3 core elements to new dictionary entry
            assert line[1] not in master_gd_dict[sample], "The ID %i is duplicated in %s" % (line[1], files)
            master_gd_dict[sample][line[1]] = OrderedDict()  # OrderedDict used such that iterating over keys will naturally reproduce original order
            master_gd_dict[sample][line[1]]["type"] = line[0]
            master_gd_dict[sample][line[1]]["id"] = line[1]
            master_gd_dict[sample][line[1]]["parent-ids"] = line[2]

            # Validation lines
            if re.match("^[A-Z]{4}$", line[0]):

                # No fields common to all Validation entries.

                if line[0] == "NOTE":  # Listed first as is most common
                    master_gd_dict[sample][line[1]]["note"] = line[3]

                elif line[0] in ["CURA", "FPOS"]:  # share 100% required fields
                    master_gd_dict[sample][line[1]]["expert"] = line[3]

                elif line[0] in ["TSEQ", "PFLP", "RFLP", "PFGE", "MASK"]:
                    master_gd_dict[sample][line[1]]["seq_id"] = line[3]  # shared among all 5

                    if line[0] == "MASK":
                        master_gd_dict[sample][line[1]]["position"] = line[4]
                        master_gd_dict[sample][line[1]]["size"] = line[5]

                    elif line[0] == "PFGE":
                        master_gd_dict[sample][line[1]]["restriction_enzyme"] = line[4]

                    else:  # TSEQ, PFLP, and RFLP share additional 4 common fields
                        master_gd_dict[sample][line[1]]["primer1_start"] = line[4]
                        master_gd_dict[sample][line[1]]["primer1_end"] = line[5]
                        master_gd_dict[sample][line[1]]["primer2_start"] = line[6]
                        master_gd_dict[sample][line[1]]["primer2_end"] = line[7]

                        if line[0] == "RFLP":
                            master_gd_dict[sample][line[1]]["enzyme"] = line[8]
                else:
                    assert line[0] == "PHYL", "The following 4 letter code mutation not recognized, please check formatting:\n%s\nFrom:\n%s" % (line[0], line)
                    master_gd_dict[sample][line[1]]["gd"] = line[3]

            # Mutation Line
            elif re.match("^[A-Z]{3}$", line[0]):

                # Store 2 core elements of each mutation to dictionary
                master_gd_dict[sample][line[1]]["seq_id"] = line[3]
                master_gd_dict[sample][line[1]]["position"] = line[4]

                # Add other required elements of each mutation to dictionary
                if line[0] in ["SNP", "INS"]:  # Both SNP and INS share 100% of required fields
                    master_gd_dict[sample][line[1]]["new_seq"] = line[5]

                elif line[0] in ["DEL", "AMP", "SUB", "CON", "INV"]:
                    master_gd_dict[sample][line[1]]["size"] = line[5]  # Required for all, only additional required for DEL, INV

                    if line[0] == "AMP":
                        master_gd_dict[sample][line[1]]["new_copy_number"] = line[6]

                    elif line[0] == "SUB":
                        master_gd_dict[sample][line[1]]["new_seq"] = line[6]

                    elif line[0] == "CON":
                        master_gd_dict[sample][line[1]]["region"] = line[6]

                else:
                    assert line[0] == "MOB", "The following 3 letter code mutation not recognized, please check formatting:\n%s\nFrom:\n%s" % (line[0], line)
                    master_gd_dict[sample][line[1]]["repeat_name"] = line[5]
                    master_gd_dict[sample][line[1]]["strand"] = line[6]
                    master_gd_dict[sample][line[1]]["duplicated_size"] = line[7]

            # Evidence lines
            elif re.match("^[A-Z]{2}$", line[0]):

                # No fields common to all evidence.

                # Add other required fields to each evidence entry in dictionary
                if line[0] in ["RA", "MC", "UN"]:
                    master_gd_dict[sample][line[1]]["seq_id"] = line[3]  # shared by all evidence save JC

                    if line[0] in ["MC", "UN"]:  # MC & UN share start and end, UN requires no additional
                        master_gd_dict[sample][line[1]]["start"] = line[4]
                        master_gd_dict[sample][line[1]]["end"] = line[5]

                        if line[0] == "MC":
                            master_gd_dict[sample][line[1]]["start_range"] = line[6]
                            master_gd_dict[sample][line[1]]["end_range"] = line[7]

                    else:  # RA
                        master_gd_dict[sample][line[1]]["position"] = line[4]
                        master_gd_dict[sample][line[1]]["insert_position"] = line[5]
                        master_gd_dict[sample][line[1]]["ref_base"] = line[6]
                        master_gd_dict[sample][line[1]]["new_base"] = line[7]

                else:
                    assert line[0] == "JC", "2 letter code found not matching 'RA, MC, UN, JC':\n%s" % line
                    master_gd_dict[sample][line[1]]["side_1_seq_id"] = line[3]
                    master_gd_dict[sample][line[1]]["side_1_position"] = line[4]
                    master_gd_dict[sample][line[1]]["side_1_strand"] = line[5]
                    master_gd_dict[sample][line[1]]["side_2_seq_id"] = line[6]
                    master_gd_dict[sample][line[1]]["side_2_position"] = line[7]
                    master_gd_dict[sample][line[1]]["side_2_strand"] = line[8]
                    master_gd_dict[sample][line[1]]["overlap"] = line[9]

            else:
                assert False, "The following line was not added to dictionary in %s:\n%s" % (files, line)

            for entry in master_gd_dict[sample][line[1]]:  # TODO decide if this is better declared rather than looping for each entry
                if entry not in named_fields:
                    named_fields.append(entry)

            # Add additional key=value pairs to dictionary with key, and value as keys and values for dictionary entries.
            if len(line) > len(master_gd_dict[sample][line[1]]):
                for entry in line[len(master_gd_dict[sample][line[1]]):]:
                    assert len(entry.split("=")) == 2, "Non key=value pair found for mutation:\n%s\n\nSpecifically: %s" % (line, entry)
                    assert entry.split("=")[0] not in named_fields, "Key in key=value pair, duplicates explicitly named field, will need to recode this 'somehow'\n%s\n\nSpecifically Duplicated: %s" % (line, entry.split("=")[0])
                    assert entry.split("=")[0] not in master_gd_dict[sample][line[1]], "Duplicated key from key=value pairs in mutation:\n%s\n\nSpecifically: %s" % (line, entry.split("=")[0])
                    master_gd_dict[sample][line[1]][entry.split("=")[0]] = entry.split("=")[1]

            if commented_out:
                master_gd_dict[sample][line[1]]["type"] = "#" + line[0]  # Replace comment as appropriate
                # will need to deal with in subsequent functions if filtering or sorting
                # Alternatively, could re-engineer to have additional key saying 'commented' or 'ignored' or something

