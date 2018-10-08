"""
Contains possible interactions with the Chado base for expressions and biomaterials
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import chado
import json

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

    def get_biomaterials(self, provider_id="", biomaterial_id="", organism_id="", biomaterial_name=""):
        """
        List biomaterials in the database

        :type organism_id: str
        :param organism_id: Limit query to the selected organism

        :type biomaterial_id: str
        :param biomaterial_id: Limit query to the selected biomaterial id

        :type provider_id: str
        :param provider_id: Limit query to the selected provider

        :type biomaterial_name: str
        :param biomaterial_name: Limit query to the selected biomaterial name


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
        if biomaterial_name:
            res = res.filter_by(name=biomaterial_name)

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
        return data


# WIP
    def add_biomaterial(self, biomaterial_name, organism_id,
                        description="", analysis_id="", biomaterial_provider="", biosample_accession="", sra_accession="",
                        bioproject_accession="", attributes={}):
        """
        Add a new biomaterial to the database

        :type biomaterial_name: str
        :param biomaterial_name: Biomaterial name

        :type organism_id: str
        :param organism_id: The id of the associated organism

        :type description: str
        :param description: Description of the biomaterial

        :type analysis_id: str
        :param analysis_id: Analysis ID

        :type biomaterial_provider: str
        :param biomaterial_provider: Biomaterial provider name

        :type biosample_accession: str
        :param biosample_accession: Biosample accession number

        :type sra_accession: str
        :param sra_accession: SRA accession number

        :type bioproject_accession: str
        :param bioproject_accession: Bioproject accession number

        :type attributes: str
        :param attributes: Custom attributes (In JSON dict form)

        :rtype: str
        :return: Biomaterial

        """
        biosourceprovider_id = None
        dbxref_id = None
        json_attributes = []

        if attributes:
            json_attributes = json.loads(attributes)

        # Add provider into table if not existing
        if biomaterial_provider:
            biosourceprovider_id = self._create_biomaterial_contact(biomaterial_provider)

        # Create biomaterial
        biomaterial_id = self._create_biomaterial(biomaterial_name, analysis_id, organism_id, biosourceprovider_id, dbxref_id, description)

        # Create DB for accession type, add to dbxref table, and then add to biomaterial_dbxref

        if sra_accession:
            dbxref_id = self._register_accession('sra', 'NCBI SRA', '', sra_accession)
            self._add_to_biomaterial_dbxref(biomaterial_id, dbxref_id)
        if bioproject_accession:
            dbxref_id = self._register_accession('bioproject', 'NCBI BioProject', '', bioproject_accession)
            self._add_to_biomaterial_dbxref(self, biomaterial_id, dbxref_id)
        if biosample_accession:
            dbxref_id = self._register_accession('biosample', 'NCBI BioSample', '', biosample_accession)
            self._add_to_biomaterial_dbxref(biomaterial_id, dbxref_id)

        # Load attributes
        if json_attributes:
            self._add_biomaterial_props(biomaterial_id, json_attributes)

        return self.get_biomaterials(biomaterial_id=biomaterial_id)[0]


    def _register_accession(self, url_name, db_name, db_description, accession):
        # Create DB if not existing, get id and update otherwise
        url = "http://www.ncbi.nlm.nih.gov/"
        urlprefix = "http://www.ncbi.nlm.nih.gov/" + url_name + "/"
        db_id = ""
        res = self.session.query(self.model.db).filter_by(name=db_name)
        if res.count():
            db_id = res.one().db_id
        if not db_id:
            database = self.model.db()
            database.name = db_name
            database.description = description
            database.urlprefix = urlprefix
            database.url = url
            self.session.add(database)
            self.session.commit()
            db_id = database.db_id
        # Update DB
        else:
            self.session.query(self.model.db).filter_by(db_id=db_id).update({
                'name' : db_name,
                'description' : db_description,
                'urlprefix': urlprefix,
                'url' : url
            })
            self.session.commit()
        # Register in dbxref

        dbxref_id = ""
        res = self.session.query(self.model.dbxref).filter_by(accession=accession, db_id=db_id)
        if res.count():
            dbxref_id = res.one().dbxref_id
        # Create if not existing,  no need for update
        if not dbxref_id:
            database = self.model.dbxref()
            database.accession = accession
            database.db_id = db_id
            self.session.add(database)
            self.session.commit()
            dbxref_id = database.dbxref_id

        return dbxref_id

    def _create_biomaterial_contact(self, contact_name):

        contact_id = ""
        res = self.session.query(self.model.contact).filter_by(name=contact_name)
        if res.count():
            contact_id = res.one().contact_id

        if not contact_id:
            contact = self.model.contact()
            contact.name = contact_name
            self.session.add(contact)
            self.session.commit()
            contact_id = contact.contact_id

        return contact_id

    def _add_to_biomaterial_dbxref(self, biomaterial_id, dbxref_id):

        biomaterial_dbxref_id = ""
        res = self.session.query(self.model.biomaterial_dbxref).filter_by(biomaterial_id=biomaterial_id, dbxref_id=dbxref_id)
        if res.count():
            biomaterial_dbxref_id = res.one().biomaterial_dbxref_id
        # Create if not existing,  no need for update
        if not biomaterial_dbxref_id:
            database = self.model.biomaterial_dbxref()
            database.biomaterial_id = biomaterial_id
            database.dbxref_id = dbxref_id
            self.session.add(database)
            self.session.commit()
            biomaterial_dbxref_id = database.biomaterial_dbxref_id

        return biomaterial_dbxref_id

    def _create_biomaterial(self, biomaterial_name, analysis_id, organism_id, biosourceprovider_id, dbxref_id, description):

        # Check if biomaterial exist
        res_biomaterial = self.session.query(self.model.biomaterial).filter_by(name=biomaterial_name, taxon_id=organism_id)
        biomaterial_id = ""
        if res_biomaterial.count():
            print("Biomaterial already exists for this organism. Will update")
            biomaterial_id = res_biomaterial.one().biomaterial_id

        if analysis_id:
            analysis_name = ""
            res_analysis = self.session.query(self.model.analysis).filter_by(analysis_id=analysis_id)
            if res_analysis.count():
                analysis_name = res_analysis.one().name
            else:
                print("Analysis with not found : will ignore")

        if (description == "" and analysis_name):
            description = 'This biomaterial: ' + biomaterial_name + ', was created for the analysis: ' + analysis_name

        if not biomaterial_id:
            biomat = self.model.biomaterial()
            biomat.name = biomaterial_name
            biomat.description = description
            biomat.taxon_id = organism_id
            biomat.biosourceprovider_id = biosourceprovider_id
            biomat.dbxref_id = dbxref_id
            self.session.add(biomat)
            self.session.commit()
            biomaterial_id = biomat.biomaterial_id
        else:
            self.session.query(self.model.biomaterial).filter_by(biomaterial_id=biomaterial_id).update({
                'description': description,
                'biosourceprovider_id': biosourceprovider_id,
                'dbxref_id': dbxref_id
            })
            self.session.commit()

        return biomaterial_id

    def _add_biomaterial_props(self, biomaterial_id, prop_dict):

        # Insert in Controlled Vocabulary table if not existing
        for key, value in prop_dict.items():
            try:
                propterm = self.ci.get_cvterm_id(key, 'biomaterial_property')
            # 'Internal DB?'
            except chado.RecordNotFoundError:
                propterm = self.ci.create_cvterm(key, 'biomaterial_property', 'internal')

            # Cache management go here (one day)
            # We need unicity for (biomaterial_id, term_id, rank)
            rank = 0
            res = self.session.query(self.model.biomaterialprop).filter_by(biomaterial_id=biomaterial_id, type_id=propterm)
            # Make sure to increment the rank if existing
            if res.count():
                rank = res.count()
            # Add to DB
            prop = self.model.biomaterialprop()
            prop.type_id=propterm
            prop.biomaterial_id=biomaterial_id
            prop.value = value
            prop.rank = rank
            self.session.add(prop)

        self.session.commit()
