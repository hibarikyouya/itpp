#!/usr/bin/python2.7

# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

import os, sys
import odjo

# a function that opens and reads a file
def file2list(filename):
    f = open(filename)
    l = f.read().splitlines()
    f.close()
    return l


def main():

    if len(sys.argv) == 2 and not os.path.isfile(sys.argv[1]):
        """If a PFAM id is given as argument, then run the search on this
        id."""
        # creates a domain object with the pfam id given as argument
        domain = odjo.domain(sys.argv[1])
        # gets the pdb ids that correspond to the pfam domain
        pdb_ids = domain.get_pdb_ids()
        # creates an interaction object/class with these informations
        # "whatever" because the last argument is useless
        inter = odjo.interaction(domain.pfam_id, pdb_ids, "whatever")
        # gets domains annotations for all pdb that correspond to the pfam
        # domain of interest
        annot = inter.get_annotations()
        # searches for interacting domains
        analysis = inter.analysis_2(annot)
        # calculate the frequencies for each pair of interacting dommains
        freq = inter.frequencies(analysis)
        # creates a graph based on the frequencies
        odjo.graph().draw(freq)
        # creates a matrix based on the frequencies
        odjo.matrixcsv().frequencies_matrix(freq)
    elif len(sys.argv) == 3 and not os.path.isfile(sys.argv[1]) and \
            sys.argv[2] == '--ipfam':
        """If a PFAM id is given as argument with the --ipfam option, then
        run the search on this id, but the search will be limited to the
        domains recorded in the iPFAM database."""
        domain = odjo.domain(sys.argv[1])
        # gets the pdb ids that correspond to the pfam domain
        pdb_ids = domain.get_pdb_ids()
        # via iPFAM, gets the domains that interact with the pfam domain
        interacting_domains = domain.get_interactions()
        # creates an interaction object/class with these informations
        inter = odjo.interaction(domain.pfam_id, pdb_ids, interacting_domains)
        # gets domains annotations for all pdb that correspond to the pfam
        # domain of interest
        annot = inter.get_annotations()
        # searches for interacting domains within the iPFAM list
        analysis = inter.analysis(annot)
        # calculate the frequencies for each pair of interacting dommains
        freq = inter.frequencies(analysis)
        # creates a graph based on the frequencies
        odjo.graph().draw(freq)
        # creates a matrix based on the frequencies
        odjo.matrixcsv().frequencies_matrix(freq)
    elif len(sys.argv) == 2 and os.path.isfile(sys.argv[1]):
        """If the argument is a file containing a list of PDB id, then run the
        program of all this PDB ids list without the iPFAM stuff."""
        # parses the file given as argument
        pdb_ids = file2list(sys.argv[1])
        # creates a structure object based on the pdb id of the list, but it
        # could have been whatever...
        struct = odjo.structure(pdb_ids[0])
        # gets domains annotations for all pdb ids of the list
        annotations = struct.get_domains_from_list(pdb_ids)
        # creates an interaction object/class
        # the arguments are useless
        inter = odjo.interaction("whatever", "whatever", "whatever")
        # searches for interacting domains
        analysis = inter.analysis_3(annotations)
        # calculate the frequencies for each pair of interacting dommains
        freq = inter.frequencies_from_list(analysis)
        # creates a graph based on the frequencies
        odjo.graph().draw_list(freq)
        # creates a matrix based on the frequencies
        odjo.matrixcsv().frequencies_matrix(freq)
    else:
        """In all other cases, prints the usage."""
        print "Usage:"
        print ""
        print "For an unique PFAM identifier:"
        print "\t./odjo.py <PFAM identifier>"
        print ""
        print "For example, you can try:"
        print "\t./odjo.py PF00931 (35 structures, fast)"
        print "\t./odjo.py PF02180 (92 structures, longer)"
        print "\t./odjo.py PF00452 (147 structures, longer)"
        print ""
        print "For an unique PFAM identifier with iPFAM:"
        print "\t./odjo.py <PFAM identifier> --ipfam"
        print ""
        print "For a list of PDB identifiers (one per line):"
        print "\t./odjo.py <file with list of PDB identifier>"
        print ""
        print "For example, you can try:"
        print "\t./odjo.py list.txt"
        sys.exit()


if __name__ == '__main__':
    main()

