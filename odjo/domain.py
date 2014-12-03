# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

""" domain class """

import csv
import urllib2
from lxml import etree

class domain:

    """
    This class was made to contain methods for the PFAM domains.

    @param pfam_id: the PFAM identifier for the domain
    @type pfam_id: string
    """

    def __init__(self, pfam_id):
        self.xml_dir = 'xml'
        # pfam domain accession number
        self.pfam_id = pfam_id
        # rcsb
        self.rcsb_resturl = 'http://www.rcsb.org/pdb/rest/'
        # ipfam variables
        self.ipfam_url = 'ftp://selab.janelia.org/pub/ipfam/Current_Release/'
        self.hetero_filename = 'heterodomain_interaction.csv'
        self.homo_filename = 'heterodomain_interaction.csv'

    def get_pdb_ids(self):
        """Returns a list of pdb id codes for a pfam access number."""
        # builds the a query in xml format
        tree = etree.parse(self.xml_dir+'/pfamtopdb_query.xml')
        query_field = tree.xpath('pfamID')
        query_field[0].text = self.pfam_id
        xml_query = etree.tostring(tree, encoding='UTF-8',
                xml_declaration=True)
        # request to rcsb server
        print "Connexion to rcsb to get pbds for "+self.pfam_id+"..."
        req = urllib2.Request(self.rcsb_resturl+'search', data=xml_query)
        socket = urllib2.urlopen(req)
        result = socket.read()
        # build the list
        pdb_ids = result.split('\n')
        pdb_ids = filter(None, pdb_ids)
        return pdb_ids


    def get_interactions(self):
        """Returns a list of interacting domains based on the ipfam
        database."""
        # the returned ids list
        ids = []
        # opens sockets to get csv files
        print "Connexion to ipfam..."
        hetero_socket = urllib2.urlopen(self.ipfam_url+self.hetero_filename)
        homo_socket = urllib2.urlopen(self.ipfam_url+self.hetero_filename)
        # parses csv files
        table = csv.reader(hetero_socket, delimiter='\t')
        for row in table:
            if row[0] == self.pfam_id:
                ids.append(row[2])
        if len(ids) > 0:
            return ids
        else:
            print "This domain has no known interactions"

