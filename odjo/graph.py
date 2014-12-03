# M2BIBS ITPP project
# Ariane Odjo & Nourdine Bah

""" graph class """

import pydot

class graph:

    """
    This class contains methods for the graphs draw with the interaction
    frequencies.

    @param interactions: a list of dict and each contains {domain_1:x,
    domain_2:y, frequency:f}
    @type interaction: list
    """

    def draw_list(self, interactions):
        """Draws an oriented graph for an interaction frequencies analysis."""
        # creates graph
        g = pydot.Dot(graph_type='digraph')
        # the nodes are the domains and the edges are annotated with the
        # frequencies
        for entry in interactions:
            node1 = pydot.Node(name=entry['domain_1'])
            node2 = pydot.Node(name=entry['domain_2'])
            edge = pydot.Edge(node1, node2, label=entry['frequency'],
                    fontsize="12")
            g.add_node(node1)
            g.add_node(node2)
            g.add_edge(edge)
        # write png
        g.write_png('graph.png')

    def draw(self, interactions):
        """Draws a non-oriented graph for an interaction frequencies
        analysis."""
        # creates graph
        g = pydot.Dot(graph_type='graph')
        # the nodes are the domains and the edges are annotated with the
        # frequencies
        for entry in interactions:
            node1 = pydot.Node(name=entry['domain_1'])
            node2 = pydot.Node(name=entry['domain_2'])
            edge = pydot.Edge(node1, node2, label=entry['frequency'],
                    fontsize="12")
            g.add_node(node1)
            g.add_node(node2)
            g.add_edge(edge)
        # write png
        g.write_png('graph.png')

