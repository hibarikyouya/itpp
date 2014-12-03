# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

""" structure class """

# common imports
import re
import urllib2
import xml.etree.ElementTree as et
from itertools import combinations

# biopyton imports
from Bio.PDB import PDBList
from Bio.PDB import PDBParser

class structure:

    """
    This class was made to contain methods for the PDB structure files.

    @param pdb_id: the PDB identifier for the structure
    @type pdb_id: string
    """

    def __init__(self, pdb_id):
        self.pdb_id = pdb_id
        self.rcsb_resturl = 'http://www.rcsb.org/pdb/rest/'

    def get_pfam_ids(self):
        """Returns pfam domain ids from a pbd structure id."""
        # uses the rcsb rest interface
        request = self.rcsb_resturl+'hmmer?structureId='+self.pdb_id
        print "Retrieves pfam informations for "+self.pdb_id+"..."
        socket = urllib2.urlopen(request)
        # the rcsb server returns a xml file which is parsed with ElementTree
        xml = socket.read()
        hmmer3 = et.fromstring(xml)
        # iterates on each pfamHit and pushes informartions in a list
        accession_nums = []
        for pfamhit in hmmer3:
            accession_num = re.sub("\..*", "", pfamhit.attrib['pfamAcc'])
            accession_nums.append(accession_num)
        return accession_nums

    def get_domains_from_list(self, pdb_ids):
        """Calls get_domains for a list of pdb ids, and returns a list."""
        annotations = []
        for pdb_id in pdb_ids:
            annotation = {}
            annotation['pdb_id'] = pdb_id
            struct = structure(pdb_id)
            domains = struct.get_domains()
            annotation['domains'] = domains
            domain_combinations = self.get_domain_combinations(domains)
            annotation['combinations'] = domain_combinations
            annotations.append(annotation)
        return annotations

    def get_domain_combinations(self, domains):
        """Returns all possible domains pairs combinations from a list."""
        domain_ids = []
        for domain in domains:
            domain_id = domain['id']
            domain_ids.append(domain_id)
        # uses itertools
        domain_combinations = combinations(domain_ids, 2)
        return domain_combinations


    def get_domains(self):
        """Returns pfam domain annotations for a pdb structure id."""
        # uses the rcsb rest interface
        request = self.rcsb_resturl+'hmmer?structureId='+self.pdb_id
        print "Connexion to rcsb for "+self.pdb_id+"..."
        socket = urllib2.urlopen(request)
        # the rcsb server returns a xml file which is parsed with ElementTree
        xml = socket.read()
        hmmer3 = et.fromstring(xml)
        # iterates on each pfamHit and pushes informartions in a list
        domains = []
        for pfamhit in hmmer3:
            domain = {}
            domain['id'] = re.sub("\..*", "", pfamhit.attrib['pfamAcc'])
            domain['name'] = pfamhit.attrib['pfamName']
            domain['chain'] = pfamhit.attrib['chainId']
            domain['start'] = pfamhit.attrib['pdbResNumStart']
            domain['end'] = pfamhit.attrib['pdbResNumEnd']
            domains.append(domain)
        return domains

    def get_pdb_file(self, directory):
        """Retrieves a pdb file."""
        pdb_list = PDBList()
        pdb_list.retrieve_pdb_file(self.pdb_id, pdir=directory)
    
    def get_model(self, pdb_id, filename):
        """Returns the first model of a pdb file."""
        parser = PDBParser(PERMISSIVE=1)
        structure = parser.get_structure(filename, pdb_id)
        model = structure[0]
        return model

    def get_residues(self, model, domain):
        """Returns the list of the residues of a domain."""
        # extracts domain informations
        chain_id = domain['chain']
        start = domain['start']
        end = domain['end']
        chain = model[chain_id]
        # builds the list of residues
        residues = []
        residues_id = range(int(start), int(end)+1)
        for residue in chain.get_list():
            if residue.get_id()[1] in residues_id:
                residues.append(residue)
        return residues

    def serial_numbers(self, atoms):
        """Returns a list of serial numbers for an atoms list."""
        serial_numbers = []
        for atom in atoms:
            serial_number = atom.get_serial_number()
            serial_numbers.append(serial_number)
        return serial_numbers

    def atoms2residues(self, atoms):
        """Returns a list of residues from a list of atoms."""
        residues = []
        for atom in atoms:
            residue = atom.get_parent()
            if residue not in residues:
                residues.append(residue)
        return residues

