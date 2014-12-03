# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

""" matrixcsv class """

import csv

class matrixcsv:

    """
    This class contains methods for CSV files I/O.
    """

    def frequencies_matrix(self, frequencies_list):
        """Writes a csv file for frequencies analysis."""
        # we take a copy
        frequencies = frequencies_list[:]
        header = self._domains_from_frequencies(frequencies)
        row_names = header[:]
        with open('matrix.csv', 'wb') as csvfile:
            writer = csv.writer(csvfile, quoting=csv.QUOTE_ALL)
            writer.writerow(['']+header)
            # checks for all domains couples
            for domain_1 in row_names:
                current_row = []
                for domain_2 in header:
                    # the frequency exists
                    if any(entry['domain_1'] == domain_1 and \
                            entry['domain_2'] == domain_2 \
                            for entry in frequencies):
                        for entry in frequencies:
                            if entry['domain_1'] == domain_1 and \
                                    entry['domain_2'] == domain_2:
                                current_row.append(float(entry['frequency']))
                    # the frequency doesn't exist so put 0
                    else:
                        current_row.append(int(0))
                writer.writerow([domain_1]+current_row)

    # private functions

    def _domains_from_frequencies(self, frequencies):
        """Returns all domains of a frequencies analysis."""
        domains = []
        for frequency in frequencies:
            # we want to have a uniquified list
            if frequency['domain_1'] not in domains:
                domains.append(frequency['domain_1'])
            if frequency['domain_2'] not in domains:
                domains.append(frequency['domain_2'])
        return domains

