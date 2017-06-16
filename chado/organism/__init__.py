"""
Contains possible interactions with the Apollo Organisms Module
"""
from chado.client import Client
from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import mapper, sessionmaker

class Organism:
    pass


class OrganismClient(Client):

    def __init__(self, *args, **kwargs):
        super(self.__class__, self).__init__(*args, **kwargs)
        organism = Table('organism', self.metadata, autoload=True)
        mapper(Organism, organism)


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

        :rtype: dict
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
