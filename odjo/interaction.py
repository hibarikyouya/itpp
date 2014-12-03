# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

""" interaction class """

# local imports
from .structure import structure
from .domain import domain
from .representation import representation

# regular imports
import re, os
from Bio.PDB import PDBList
from Bio.PDB import PDBParser
from Bio.PDB import Selection
from Bio.PDB import NeighborSearch

class interaction:

    """
    This class was originally made to process informations for a PFAM domain.

    @param pfam_id: the identifier of the PFAM domain
    @type pfam_id: string

    @param pdb_ids: the list of PDB structure that contains the PFAM domain
    @type pdb_ids: list

    @param pfam_ids: the list of PFAM domains that interacts with the domain,
    this list comes from the iPFAM database
    @type pfam_ids: list

    """

    def __init__(self, pfam_id, pdb_ids, pfam_ids):
        self.pfam_id = pfam_id
        self.pdb_ids = pdb_ids
        # a pfam domain ids list for iPFAM
        self.pfam_ids = pfam_ids
        self.directory = 'pdb_and_png'

    def get_pdb_files(self):
        """Retrieves all pdb files corresponding the domains of interest."""
        pdb_list = PDBList()
        for pdb in self.pdb_ids:
            # put files in the directory pointed by the self.directory variable
            pdb_list.retrieve_pdb_file(pdb, pdir=self.directory)

    def get_filtered_pdb_ids(self):
        """Returns the list of pdb ids corresponding to an interaction."""
        ids = []
        # get pfam ids for each pdb file
        for pdb_id in self.pdb_ids:
            acc_nums = structure(pdb_id).get_pfam_ids()
            # check if a pfam id is in the known domains list
            for acc_num in acc_nums:
                if acc_num in self.pfam_ids and pdb_id not in ids:
                    ids.append(pdb_id)
        return ids

    def get_annotations(self):
        """Returns a list of pfam domains annotations for a list of PDB."""
        pdbs = []
        for pdb_id in self.pdb_ids:
            # a function from the structure class
            domains = structure(pdb_id).get_domains()
            entry = {}
            entry['pdb_id'] = pdb_id
            entry['domains'] = domains
            pdbs.append(entry)
        annotations = {}
        annotations['pfam_id'] = self.pfam_id
        annotations['pdbs'] = pdbs
        return annotations

    def analysis(self, annotations):
        """Returns an interaction analysis dict based on iPFAM."""
        analysis = {}
        analysis['interactants'] = []
        pfam_id = annotations['pfam_id']
        pdbs = annotations['pdbs']
        for pdb in pdbs:
            pdb_id = pdb['pdb_id']
            domains = pdb['domains']
            # just extracts the informations for the domain passed as an
            # argument
            main_domain = self._get_domain_of_interest(domains)
            main_name = main_domain['name']
            analysis['name'] = main_name
            for domain in domains:
                name = domain['name']
                # checks the iPFAM list
                if domain['id'] in self.pfam_ids:
                    filename = self._id2filename(pdb_id)
                    rel_path = self.directory+'/'+filename
                    if not os.path.isfile(rel_path):
                        # downloads the PDB file
                        structure(pdb_id).get_pdb_file(self.directory)
                    # searches for interactions between this two domains and
                    # in this PDB structure
                    inter = self.interaction(rel_path, pdb_id, main_domain,
                            domain)
                    if inter:
                        # export an annotated structure in a png, but we have
                        # problems with that : pymol does the job for a few
                        # structures but bugs for a lot of structures.
                        # Please uncomment if you wanna test
                        #representation().export_png(rel_path, pdb_id, inter)
                        analysis['interactants'].append(name)
        return analysis
                    
    def analysis_2(self, annotations):
        """Returns an interaction analysis dict without a priori."""
        analysis = {}
        analysis['interactants'] = []
        pfam_id = annotations['pfam_id']
        pdbs = annotations['pdbs']
        for pdb in pdbs:
            pdb_id = pdb['pdb_id']
            domains = pdb['domains']
            # just extracts the informations for the domain passed as an
            # argument
            main_domain = self._get_domain_of_interest(domains)
            main_name = main_domain['name']
            analysis['name'] = main_name
            for domain in domains:
                # to not test the main domain with itself
                if domain['id'] != main_domain['id']:
                    name = domain['name']
                    filename = self._id2filename(pdb_id)
                    rel_path = self.directory+'/'+filename
                    if not os.path.isfile(rel_path):
                        # downloads the PDB file
                        structure(pdb_id).get_pdb_file(self.directory)
                    # searches for interactions between this two domains and
                    # in this PDB structure
                    inter = self.interaction(rel_path, pdb_id, main_domain,
                            domain)
                    if inter:
                        # export an annotated structure in a png, but we have
                        # problems with that : pymol does the job for a few
                        # structures but bugs for a lot of structures.
                        # Please uncomment if you wanna test
                        #representation().export_png(rel_path, pdb_id, inter)
                        analysis['interactants'].append(name)
        return analysis

    def analysis_3(self, annotations):
        analysis = []
        for annotation in annotations:
            pdb_id = annotation['pdb_id']
            domains = annotation['domains']
            # converts a PDB id to a filename
            filename = self._id2filename(pdb_id)
            rel_path = self.directory+'/'+filename
            if not os.path.isfile(rel_path):
                # downloads the PDB file
                structure(pdb_id).get_pdb_file(self.directory)
            # all possible domains combinations (2) for an annotation
            combs = annotation['combinations']
            for comb in combs:
                # extracts domain annotations for an id of the combination
                domain_1 = self._get_domain_from_id(comb[0], domains)
                domain_2 = self._get_domain_from_id(comb[1], domains)
                name_1 = domain_1['name']
                name_2 = domain_2['name']
                # searches for interactions between this two domains and
                # in this PDB structure
                inter = self.interaction(rel_path, pdb_id, domain_1, domain_2)
                if inter:
                    # export an annotated structure in a png, but we have
                    # problems with that : pymol does the job for a few
                    # structures but bugs for a lot of structures.
                    # Please uncomment if you wanna test
                    #representation().export_png(rel_path, pdb_id, inter)
                    # appends the interaction in a list of entries
                    self._append_interactants(name_1, name_2, analysis)
                    self._append_interactants(name_2, name_1, analysis)
        return analysis

    def interaction(self, pdb_id, filename, domain_1, domain_2):
        """Returns a dict with informations (atoms, residues...) if two domains
        interact with each other, and returns False if not."""
        print "Searching for interactions in "+pdb_id+"..."
        # creates a strucuture object/class to extract atoms of the two domains
        model = structure(pdb_id).get_model(pdb_id, filename)
        residues_1 = structure(pdb_id).get_residues(model, domain_1)
        residues_2 = structure(pdb_id).get_residues(model, domain_2)
        atoms_1 = Selection.unfold_entities(residues_1, 'A')
        atoms_2 = Selection.unfold_entities(residues_2, 'A')
        # gets the serial numbers of the atoms
        numbers_1 = structure(pdb_id).serial_numbers(atoms_1)
        numbers_2 = structure(pdb_id).serial_numbers(atoms_2)
        # the search starts here !
        atoms = Selection.unfold_entities(model, 'A')
        nsearch = NeighborSearch(atoms)
        interacting_atoms_1 = []
        interacting_atoms_2 = []
        for atom in atoms:
            if atom.get_serial_number() in numbers_1:
                point = atom.get_coord()
                # This is how we detect an interaction, we put 5 angstroms
                # here.
                # This is the simplest method we can use, and we're not sure
                # that it is correct.
                # Originally we have planned to go further by doing a surface
                # and accesssion analysis, but we had no time.
                # We hope we can talk about that during the talk.
                neighbors = nsearch.search(point, 5)
                for neighbor in neighbors:
                    if neighbor.get_serial_number() in numbers_2:
                        interacting_atoms_2.append(neighbor)
                        if atom not in interacting_atoms_1:
                            interacting_atoms_1.append(atom)
        # returns a dict with all residues and atoms
        if len(interacting_atoms_2) > 0:
            infos = {}
            infos['1'] = {}
            infos['2'] = {}
            # just get the parent residues for the list of atoms
            interacting_residues_1 = structure(pdb_id).atoms2residues(
                    interacting_atoms_1)
            interacting_residues_2 = structure(pdb_id).atoms2residues(
                    interacting_atoms_2)
            infos['1']['atoms'] = interacting_atoms_1
            infos['2']['atoms'] = interacting_atoms_2
            infos['1']['residues'] = interacting_residues_1
            infos['2']['residues'] = interacting_residues_2
            return infos
        else: return False

    def frequencies(self, analysis):
        """Returns the frequencies of interaction pairs based on the number
        of occurences of each domain in the interacting domain list."""
        frequencies = []
        interactants = analysis['interactants']
        # uniquify the list
        unique_interactants = list(set(interactants))
        for interactant in unique_interactants:
            occurence_num = interactants.count(interactant)
            frequency = float(occurence_num) / float(len(interactants))
            entry = {}
            entry['domain_1'] = analysis['name']
            entry['domain_2'] = interactant
            entry['frequency'] = str(round(frequency, 2))
            frequencies.append(entry)
        return frequencies

    def frequencies_from_list(self, analysis):
        """Calls the frequencies() function for a list of interactions."""
        frequencies = []
        for entry in analysis:
            entry_freq = self.frequencies(entry)
            frequencies += entry_freq
        return frequencies


    # private functions
    
    def _ids2filenames(self, pdb_ids):
        """Returns a list of filenames from a pdb id list."""
        ids = []
        for pdb_id in pdb_ids:
            pdb_id = pdb_id.lower()
            pdb_id = re.sub('^', 'pdb', pdb_id)
            pdb_id = re.sub('$', '.ent', pdb_id)
            ids.append(pdb_id)
        return ids

    def _id2filename(self, pdb_id):
        """Returns a list of filenames from a pdb id list."""
        filename = pdb_id.lower()
        filename = re.sub('^', 'pdb', filename)
        filename = re.sub('$', '.ent', filename)
        return filename

    def _get_domain_of_interest(self, domains):
        """Returns the domain of interest from a list of domains."""
        for domain in domains:
            if domain['id'] == self.pfam_id:
                return domain

    def _get_domain_from_id(self, pfam_id, domains):
        """Returns the domains list corresponding to a PFAM id."""
        for domain in domains:
            if domain['id'] == pfam_id:
                return domain

    def _append_interactants(self, name_1, name_2, analysis):
        """Takes two domains and a list of dict that contains all interacting
        domains for each domain. This function creates the dict if necessary
        and adds the interacting domains."""
        # the entry for a domain already exists
        if any(entry['name'] == name_1 for entry in analysis):
            for entry in analysis:
                if entry['name'] == name_1:
                    entry['interactants'].append(name_2)
        # the entry for a domain doesn't exist
        else:
            new_entry = {}
            new_entry['name'] = name_1
            new_entry['interactants'] = []
            new_entry['interactants'].append(name_2)
            analysis.append(new_entry)

