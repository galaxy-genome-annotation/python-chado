"""
Contains possible interactions with the Chado Features
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from chado.client import Client
from future import standard_library

standard_library.install_aliases()


class FeatureClient(Client):
    """
    Access to the chado features
    """

    def load_fasta(self, fasta, analysis, organism, sequence_type):
        """
        Load features from a fasta file

        :type fasta: str
        :param fasta: Path to the Fasta file to load

        :type analysis: int
        :param analysis: Analysis ID

        :type organism: int
        :param organism: Organism ID

        :type sequence_type: str
        :param sequence_type: Sequence type

        :rtype: None
        :return: None
        """
        raise Exception("Not implemented")

    def load_gff(self, gff, analysis, organism, sequence_type):
        """
        Load features from a gff file

        :type gff: str
        :param gff: Path to the Fasta file to load

        :type analysis: int
        :param analysis: Analysis ID

        :type organism: int
        :param organism: Organism ID

        :type sequence_type: str
        :param sequence_type: Sequence type

        :rtype: None
        :return: None
        """
        raise Exception("Not implemented")
