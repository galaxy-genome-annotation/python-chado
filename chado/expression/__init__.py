"""
Contains possible interactions with the Chado base for expressions and biomaterials
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from chado.client import Client

from future import standard_library

standard_library.install_aliases()


class ExpressionClient(Client):
    """
    Interact with expressions
    """


#WIP separator)
#    def add_expression(self, organism, analysis, match_type, file_path,
#                       separator="tab", biomaterial_provider=None, array_design=None, assay_id=None,
#                       acquisition_id=None, quantification_id=None, file_extension=None,
#                       start_regex=None, stop_regex=None, use_column=False):

#        return ""


#To test

    def get_biomaterials(self, provider_id="", biomaterial_id="", organism_id="", dbxref_id=""):
        """
        List biomaterials in the database

        :type organism_id: str
        :param organism_id: Limit query to the selected organism

        :type biomaterial_id: str
        :param biomaterial_id: Limit query to the selected biomaterial

        :type provider_id: str
        :param provider_id: Limit query to the selected provider

        :type dbxref_id: str
        :param dbxref_id: Limit query to the selected ref

        :rtype: list
        :return: List of biomaterials
        """
        res = self.session.query(self.model.biomaterial)

        if biomaterial_id:
            res = res.filter_by(biomaterial_id=biomaterial_id)
        if provider_id:
            res = res.filter_by(biosourceprovider_id=provider_id)
        if organism_id:
            res = res.filter_by(taxon_id=organism_id)
        if dbxref_id:
            res = res.filter_by(dbxref_id=dbxref_id)

        data = []
        for biomat in res:
            data.append({
                'provider_id': biomat.biosourceprovider_id,
                'description': biomat.description,
                'biomaterial_id': biomat.biomaterial_id,
                'organism_id': biomat.taxon_id,
                'dbxref_id': biomat.dbxref_id,
                'name': biomat.name,
            })
#        return data
        return(res[0].items)


# WIP
    def add_biomaterial(self, organism_id, biomaterial_name, organism_name="", provider="", accession="", sra_accession="", bioproject_accession="", attributes={}):
        """
        Add a new biomaterial to the database

        :type organism_id: str
        :param organism_id: The id of the associated organism

        :type biomaterial_name: str
        :param biomaterial_name: Biomaterial name

        :type organism_name: str
        :param organism_name: Orgnaism name (can be different from the organism ID provided)

        :type provider: str
        :param provider: Biomaterial provider name

        :type accession: str
        :param accession: Biomaterial accession

        :type sra_accession: str
        :param sra_accession: SRA acession

        :type bioproject_accession: str
        :param bioproject_accession: Bioproject accession

        :type attributes: dict
        :param attributes: Custom attributes (In JSON dict form)

        :rtype: dict
        :return: Job information
        """

        biomat = self.model.biomaterial()
        return ""
