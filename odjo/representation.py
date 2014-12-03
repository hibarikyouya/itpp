# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

""" representation class """

import re
import pymol
from Bio.PDB import PDBParser

class representation:

    """
    This class contains methods to export an image based on an annotated PDB
    structure (pymol).
    """

    def export_png(self, rel_path, pdb_id, infos):
        """Uses pymol to export a png representing a PDB structure. The image
        displays the chains that contains the domains (cartoon) and the
        interacting residues (surface). THE CHOICE HAS BEEN MADE TO NOT
        HIGHLIGHT THE DOMAINS. """

        # extracts informations for the chain and the residues
        domain_1 = self._extract_info(infos, '1')
        domain_2 = self._extract_info(infos, '2')

        # gets informations
        chain_1 = domain_1['chain']
        chain_2 = domain_2['chain']
        residues_1 = '+'.join(domain_1['ids'])
        residues_2 = '+'.join(domain_2['ids'])

        # build the select command for pymol
        chain_1_obj = 'chain '+chain_1+' and (not resi '+residues_1+')'
        chain_2_obj = 'chain '+chain_2+' and (not resi '+residues_2+')'
        residues_1_obj = 'chain '+chain_1+' and resi '+residues_1
        residues_2_obj = 'chain '+chain_1+' and resi '+residues_2

        # run pymol without its gui
        pymol.pymol_argv = ['pymol','-qc']
        pymol.finish_launching()

        # the file
        pymol.cmd.load(rel_path, pdb_id)

        pymol.cmd.disable('all')
        pymol.cmd.enable(pdb_id)
        # background
        pymol.cmd.bg_color(color='white')
        pymol.cmd.hide('all')
        pymol.cmd.deselect()
        # objects
        pymol.cmd.create('chain_1', chain_1_obj)
        pymol.cmd.create('chain_2', chain_2_obj)
        pymol.cmd.create('res_1', residues_1_obj)
        pymol.cmd.create('res_2', residues_2_obj)
        # colors
        pymol.cmd.color('lightblue', 'chain_1')
        pymol.cmd.color('nitrogen', 'res_1')
        pymol.cmd.color('yelloworange', 'chain_2')
        pymol.cmd.color('red', 'res_2')
        # representation type
        pymol.cmd.show('cartoon', 'chain_1')
        pymol.cmd.show('cartoon', 'chain_2')
        pymol.cmd.show('surface', 'res_1')
        pymol.cmd.show('surface', 'res_2')
        # zoom
        pymol.cmd.create('final', 'chain_1 or chain_2 or res_1 or res_2')
        pymol.cmd.zoom('final')
        # export
        pymol.cmd.png(re.sub('.ent$', '.png', rel_path))

        # quit
        pymol.cmd.quit()

    # private functions

    def _extract_info(self, infos, domain):
        """Returns the chain and the residues corresponding to a domain in a
        PDB file."""
        ids = []
        residues = infos[domain]['residues']
        # the residues are necessarily in the same chain, so :
        chain = residues[0].get_parent().get_id()
        for residue in residues:
            number = str(residue.get_id()[1])
            if number not in ids:
                ids.append(number)
        return {'chain': chain, 'ids': ids}

