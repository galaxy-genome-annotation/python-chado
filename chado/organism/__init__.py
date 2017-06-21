"""
Contains possible interactions with the Apollo Organisms Module
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from chado.client import Client
from chado.models import *


class OrganismClient(Client):
    """
    Access to the chado organism table
    """

    def add_organism(self, genus, common, abbr, species=None, comment=None):
        """
        Add a new organism to the Chado database

        :type genus: str
        :param genus: The genus of the organism

        :type common: str
        :param common: The common name of the organism

        :type abbr: str
        :param abbr: The abbreviation of the organism

        :type species: str
        :param species: The species of the organism

        :type comment: str
        :param comment: A comment / description

        :rtype: dict
        :return: Organism information
        """

        # check if the organism exists
        res = self.session.query(Organism).filter_by(common_name=common)

        if (res.count() > 0):
            raise Exception("Found a preexisting organism with the same attributes in the database")

        org = Organism()
        org.abbreviation = abbr
        org.genus = genus
        org.species = species
        org.common_name = common
        org.comment = comment
        self.session.add(org)
        self.session.commit()
        return {
            'abbreviation': org.abbreviation,
            'comment': org.comment,
            'common_name': org.common_name,
            'genus': org.genus,
            'organism_id': org.organism_id,
            'species': org.species,
        }

    def get_organisms(self, genus=None, species=None, common=None, abbr=None, comment=None):
        """
        Get all organisms

        :type genus: str
        :param genus: genus filter

        :type species: str
        :param species: species filter

        :type common: str
        :param common: common filter

        :type abbr: str
        :param abbr: abbr filter

        :type comment: str
        :param comment: comment filter

        :rtype: list of dict
        :return: Organisms information
        """

        # check if the organism exists
        res = self.session.query(Organism)
        if genus:
            res = res.filter_by(genus=genus)
        if species:
            res = res.filter_by(species=species)
        if common:
            res = res.filter_by(common_name=common)
        if abbr:
            res = res.filter_by(abbreviation=abbr)
        if comment:
            res = res.filter_by(comment=comment)

        data = []
        for org in res:
            data.append({
                'organism_id': org.organism_id,
                'genus': org.genus,
                'species': org.species,
                'abbreviation': org.abbreviation,
                'common_name': org.common_name,
                'comment': org.comment,
            })
        return data

    def delete_all_organisms(self, confirm=False):
        """
        Delete all organisms

        :type confirm: bool
        :param confirm: Confirm that you really do want to delete ALL of the organisms.

        :rtype: None
        :return: None
        """

        # check if the organism exists
        res = self.session.query(Organism).delete()
        self.session.commit()
        return res
