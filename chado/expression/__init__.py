"""
Contains possible interactions with the Chado base for expressions and biomaterials
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import csv
import json

import chado
from chado.client import Client

from chakin.io import warn

from future import standard_library

standard_library.install_aliases()


class ExpressionClient(Client):
    """
    Interact with expressions
    """

    def get_biomaterials(self, provider_id="", biomaterial_id="", organism_id="", biomaterial_name="", analysis_id=""):
        """
        List biomaterials in the database

        :type organism_id: int
        :param organism_id: Limit query to the selected organism

        :type biomaterial_id: int
        :param biomaterial_id: Limit query to the selected biomaterial id

        :type provider_id: int
        :param provider_id: Limit query to the selected provider

        :type biomaterial_name: str
        :param biomaterial_name: Limit query to the selected biomaterial name

        :type analysis_id: int
        :param analysis_id: Limit query to the selected analysis_id

        :rtype: list
        :return: List of biomaterials
        """

        if analysis_id:
            res = self.session.query(self.model.biomaterial) \
                .join(self.model.assay_biomaterial, self.model.biomaterial.biomaterial_id == self.model.assay_biomaterial.biomaterial_id) \
                .join(self.model.assay, self.model.assay_biomaterial.assay_id == self.model.assay.assay_id) \
                .join(self.model.acquisition, self.model.assay.assay_id == self.model.acquisition.assay_id) \
                .join(self.model.quantification, self.model.acquisition.acquisition_id == self.model.quantification.acquisition_id) \
                .filter(self.model.quantification.analysis_id == analysis_id)
        else:
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

    def add_biomaterial(self, biomaterial_name, organism_id,
                        description="", biomaterial_provider="", biosample_accession="", sra_accession="",
                        bioproject_accession="", attributes={}):
        """
        Add a new biomaterial to the database

        :type biomaterial_name: str
        :param biomaterial_name: Biomaterial name

        :type organism_id: int
        :param organism_id: The id of the associated organism

        :type description: str
        :param description: Description of the biomaterial

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

        :rtype: dict
        :return: Biomaterial details

        """
        biosourceprovider_id = None
        dbxref_id = None
        json_attributes = []
        # Set to none because it only change the description, and does not add it to DB. Misleading
        analysis_id = None

        if attributes:
            if isinstance(attributes, str):
                json_attributes = json.loads(attributes)
            else:
                json_attributes = attributes

        # Add provider into table if not existing
        if biomaterial_provider:
            biosourceprovider_id = self._create_biomaterial_contact(biomaterial_provider)

        # Create biomaterial
        biomaterial_id = self._create_biomaterial(biomaterial_name, organism_id,
                                                  analysis_id=analysis_id, biosourceprovider_id=biosourceprovider_id,
                                                  dbxref_id=dbxref_id, description=description)

        # Create DB for accession type, add to dbxref table, and then add to biomaterial_dbxref

        if sra_accession:
            dbxref_id = self._register_accession('sra', 'NCBI SRA', '', sra_accession)
            self._add_to_biomaterial_dbxref(biomaterial_id, dbxref_id)
        if bioproject_accession:
            dbxref_id = self._register_accession('bioproject', 'NCBI BioProject', '', bioproject_accession)
            self._add_to_biomaterial_dbxref(biomaterial_id, dbxref_id)
        if biosample_accession:
            dbxref_id = self._register_accession('biosample', 'NCBI BioSample', '', biosample_accession)
            self._add_to_biomaterial_dbxref(biomaterial_id, dbxref_id)

        # Load attributes
        if json_attributes:
            self._add_biomaterial_props(biomaterial_id, json_attributes)

        self.session.commit()

        return self.get_biomaterials(biomaterial_id=biomaterial_id)[0]

    def add_expression(self, organism_id, analysis_id, file_path, separator="\t", unit=None, query_type="mRNA", match_on_name=False, re_name=None, skip_missing=False):
        """
        Add an expression matrix file to the database

        :type organism_id: int
        :param organism_id: The id of the associated organism

        :type analysis_id: int
        :param analysis_id: The id of the associated analysis

        :type file_path: str
        :param file_path: File path

        :type separator: str
        :param separator: Separating character in the matrix file (ex : ','). Default character is tab.

        :type unit: str
        :param unit: The units associated with the loaded values (ie, FPKM, RPKM, raw counts)

        :type query_type: str
        :param query_type: The feature type (e.g. \'gene\', \'mRNA\', 'polypeptide', \'contig\') of the query. It must be a valid Sequence Ontology term.

        :type match_on_name: bool
        :param match_on_name: Match features using their name instead of their uniquename

        :type re_name: str
        :param re_name: Regular expression to extract the feature name from the input file (first capturing group will be used).

        :type skip_missing: bool
        :param skip_missing: Skip lines with unknown features or GO id instead of aborting everything.

        :rtype: str
        :return: Number of expression data loaded

        """
        self._reset_cache()
        seqterm = self.ci.get_cvterm_id(query_type, 'sequence')
        self._init_feature_cache(organism_id, seqterm, match_on_name)
        self._init_analysisfeature_cache(analysis_id)

        uniq_name = self._get_unique_name(organism_id, analysis_id)
        contact_id = self._create_generic_contact()
        arraydesign_id = self._create_generic_arraydesign(contact_id)

        results = self._process_matrix_file(file_path, separator)

        # Create all biomaterials, get quantification ID
        quant_list = []
        for biomaterial in results["biomaterial_list"]:
            quant_list.append(self._expression_create_biomaterial_structure(biomaterial, organism_id, analysis_id, arraydesign_id, contact_id, uniq_name, unit))

        # Manage features
        num = 0
        for index, feature in enumerate(results["feature_list"]):
            self._manage_feature_expression(feature, results["data"][index], quant_list, organism_id, analysis_id, arraydesign_id, query_type, re_name, skip_missing)
            num += len(results["data"][index])

        self.session.commit()

        self._reset_cache()

        return "%s expression values added" % num

    def delete_biomaterials(self, names=None, ids=None, organism_id="", analysis_id=""):
        """
        Will delete biomaterials based on selector. Only one selector will be used.

        :type names: str
        :param names: JSON list of biomaterial names to delete.

        :type ids: str
        :param ids: JSON list of biomaterial ids to delete.

        :type organism_id: int
        :param organism_id: Delete all biomaterial associated with this organism id.

        :type analysis_id: int
        :param analysis_id: Delete all biomaterial associated with this analysis id.

        :rtype: str
        :return: Number of deleted biomaterials

        """
        if not (names or ids or organism_id or analysis_id):
            return("Please select one filter")

        if names:
            if isinstance(names, str):
                names = json.loads(names)

        if ids:
            if isinstance(ids, str):
                ids = json.loads(ids)

        if names:
            res = self.session.query(self.model.biomaterial).filter(self.model.biomaterial.name.in_(names))
            if not res.count():
                return("No biomaterials with given names were found")

        elif ids:
            res = self.session.query(self.model.biomaterial).filter(self.model.biomaterial.biomaterial_id.in_(ids))
            if not res.count():
                return ("No biomaterials with these ids were found")

        elif organism_id:
            res = self.session.query(self.model.biomaterial).filter_by(taxon_id=organism_id)
            if not res.count():
                return("No biomaterials associated with this organism were found")

        elif analysis_id:
            res = self.session.query(self.model.biomaterial.biomaterial_id) \
                .join(self.model.assay_biomaterial, self.model.biomaterial.biomaterial_id == self.model.assay_biomaterial.biomaterial_id) \
                .join(self.model.assay, self.model.assay_biomaterial.assay_id == self.model.assay.assay_id) \
                .join(self.model.acquisition, self.model.assay.assay_id == self.model.acquisition.assay_id) \
                .join(self.model.quantification, self.model.acquisition.acquisition_id == self.model.quantification.acquisition_id) \
                .filter(self.model.quantification.analysis_id == analysis_id)
            if not res.count():
                return("No biomaterials associated with this analysis ID were found")
            else:
                print(str(res.count()) + " biomaterials associated with this analysis ID were found")
                ids = [biomat[0] for biomat in res]
                res = self.session.query(self.model.biomaterial).filter(self.model.biomaterial.biomaterial_id.in_(ids))

        number = str(res.delete(synchronize_session=False))
        self.session.commit()
        return number + " biomaterial(s) deleted"

    def delete_all_biomaterials(self, confirm=False):
        """
        Delete all biomaterials

        :type confirm: bool
        :param confirm: Confirm that you really do want to delete ALL of the biomaterials.

        :rtype: None
        :return: None

        """
        if confirm:
            res = self.session.query(self.model.biomaterial).delete(synchronize_session=False)
            self.session.commit()
            return res
        else:
            print("No --confirm option set. Not doing anything")

    def _create_generic_contact(self):
        return self._create_biomaterial_contact("Not provided", "Caution: This is a generic contact created by the expression module. Delete with caution.")

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
            database.description = db_description
            database.urlprefix = urlprefix
            database.url = url
            self.session.add(database)
            self.session.flush()
            self.session.refresh(database)
            db_id = database.db_id
        # Update DB
        else:
            self.session.query(self.model.db).filter_by(db_id=db_id).update({
                'name': db_name,
                'description': db_description,
                'urlprefix': urlprefix,
                'url': url
            })

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
            self.session.flush()
            self.session.refresh(database)
            dbxref_id = database.dbxref_id

        return dbxref_id

    def _create_biomaterial_contact(self, contact_name, description=""):
        contact_id = ""
        res = self.session.query(self.model.contact).filter_by(name=contact_name)
        if res.count():
            contact_id = res.one().contact_id

        if not contact_id:
            contact = self.model.contact()
            contact.name = contact_name
            contact.description = description
            self.session.add(contact)
            self.session.flush()
            self.session.refresh(contact)
            contact_id = contact.contact_id
        else:
            if not res.one().description == description:
                res.one().description = description
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
            self.session.flush()
            self.session.refresh(database)
            biomaterial_dbxref_id = database.biomaterial_dbxref_id

        return biomaterial_dbxref_id

    def _create_biomaterial(self, biomaterial_name, organism_id, analysis_id=None,
                            biosourceprovider_id=None, dbxref_id=None, description=None):

        # Check if biomaterial exist
        res_biomaterial = self.session.query(self.model.biomaterial).filter_by(name=biomaterial_name)
        biomaterial_id = ""
        if res_biomaterial.count():
            biomaterial_id = res_biomaterial.one().biomaterial_id
            # Do not update if not set and existing in DB
            if not description:
                description = res_biomaterial.one().description
            if not dbxref_id:
                dbxref_id = res_biomaterial.one().dbxref_id
            if not biosourceprovider_id:
                res_biomaterial.one().biosourceprovider_id

        analysis_name = ""
        if analysis_id:
            res_analysis = self.session.query(self.model.analysis).filter_by(analysis_id=analysis_id)
            if res_analysis.count():
                analysis_name = res_analysis.one().name
            else:
                warn("Analysis not found: will ignore")

        if (not description and analysis_name):
            description = 'This biomaterial: ' + biomaterial_name + ', was created for the analysis: ' + analysis_name

        if not biomaterial_id:
            biomat = self.model.biomaterial()
            biomat.name = biomaterial_name
            biomat.description = description
            biomat.taxon_id = organism_id
            biomat.biosourceprovider_id = biosourceprovider_id
            biomat.dbxref_id = dbxref_id
            self.session.add(biomat)
            self.session.flush()
            self.session.refresh(biomat)
            biomaterial_id = biomat.biomaterial_id
        else:
            self.session.query(self.model.biomaterial).filter_by(biomaterial_id=biomaterial_id).update({
                'description': description,
                'biosourceprovider_id': biosourceprovider_id,
                'dbxref_id': dbxref_id
            })

        return biomaterial_id

    def _add_biomaterial_props(self, biomaterial_id, prop_dict):

        # Insert in Controlled Vocabulary table if not existing
        for key, value in prop_dict.items():
            try:
                propterm = self.ci.get_cvterm_id(key, 'biomaterial_property')
            except chado.RecordNotFoundError:
                propterm = self.ci.create_cvterm(key, 'biomaterial_property', 'tripal')

            # Cache management go here (one day)
            # We need unicity for (biomaterial_id, term_id, rank)
            rank = 0
            res = self.session.query(self.model.biomaterialprop).filter_by(biomaterial_id=biomaterial_id, type_id=propterm)
            # Make sure to increment the rank if existing
            if res.count():
                rank = res.count()
            # Add to DB
            prop = self.model.biomaterialprop()
            prop.type_id = propterm
            prop.biomaterial_id = biomaterial_id
            prop.value = value
            prop.rank = rank
            self.session.add(prop)

    def _create_generic_arraydesign(self, contact_id):
        res = self.session.query(self.model.arraydesign).filter_by(name='Not provided', manufacturer_id=contact_id)
        array_id = ""
        if res.count():
            array_id = res.one().arraydesign_id

        if not array_id:
            array = self.model.arraydesign()
            array.name = 'Not provided'
            array.description = 'Caution: This is a generic arraydesign created by the expression module. Delete with caution.'
            array.manufacturer_id = contact_id
            array.platformtype_id = 1
            self.session.add(array)
            self.session.flush()
            self.session.refresh(array)
            array_id = array.arraydesign_id

        return array_id

    def _create_generic_assay(self, contact_id, arraydesign_id, uniq_name):
        res = self.session.query(self.model.assay).filter_by(name=uniq_name)
        assay_id = ""
        if res.count():
            assay_id = res.one().assay_id

        if not assay_id:
            assay = self.model.assay()
            assay.name = uniq_name
            assay.arraydesign_id = arraydesign_id
            assay.operator_id = contact_id
            self.session.add(assay)
            self.session.flush()
            self.session.refresh(assay)
            assay_id = assay.assay_id

        return assay_id

    def _create_generic_acquisition(self, assay_id, uniq_name):
        res = self.session.query(self.model.acquisition).filter_by(name=uniq_name)
        acqui_id = ""
        if res.count():
            acqui_id = res.one().acquisition_id

        if not acqui_id:
            acqui = self.model.acquisition()
            acqui.name = uniq_name
            acqui.assay_id = assay_id
            self.session.add(acqui)
            self.session.flush()
            self.session.refresh(acqui)
            acqui_id = acqui.acquisition_id

        return acqui_id

    def _create_generic_quantification(self, acquisition_id, analysis_id, uniq_name):
        res = self.session.query(self.model.quantification).filter_by(name=uniq_name)
        quantification_id = ""
        if res.count():
            quantification_id = res.one().quantification_id

        if not quantification_id:
            quantification = self.model.quantification()
            quantification.name = uniq_name
            quantification.acquisition_id = acquisition_id
            quantification.analysis_id = analysis_id
            self.session.add(quantification)
            self.session.flush()
            self.session.refresh(quantification)
            quantification_id = quantification.quantification_id

        return quantification_id

    def _create_generic_channel(self):
        res = self.session.query(self.model.channel).filter_by(name='Not provided')
        channel_id = ""
        if res.count():
            channel_id = res.one().channel_id

        if not channel_id:
            channel = self.model.channel()
            channel.name = 'Not provided'
            channel.definition = 'Caution: This is a generic channel created by the expression module'
            self.session.add(channel)
            self.session.flush()
            self.session.refresh(channel)
            channel_id = channel.channel_id

        return channel_id

    def _create_assay_biomaterial(self, assay_id, biomaterial_id, channel_id):
        res = self.session.query(self.model.assay_biomaterial).filter_by(rank=1, assay_id=assay_id, biomaterial_id=biomaterial_id, channel_id=channel_id)
        assay_biomaterial_id = ""
        if res.count():
            assay_biomaterial_id = res.one().assay_biomaterial_id

        if not assay_biomaterial_id:
            assay_biomaterial = self.model.assay_biomaterial()
            assay_biomaterial.rank = 1
            assay_biomaterial.assay_id = assay_id
            assay_biomaterial.biomaterial_id = biomaterial_id
            assay_biomaterial.channel_id = channel_id
            self.session.add(assay_biomaterial)
            self.session.flush()
            self.session.refresh(assay_biomaterial)
            assay_biomaterial_id = assay_biomaterial.assay_biomaterial_id

        return assay_biomaterial_id

    def _get_unique_name(self, organism_id, analysis_id):
        # Generate 'unique' name for couple (organism, analysis)
        uniq_name = " from "
        res = self.session.query(self.model.organism).filter_by(organism_id=organism_id)
        if not res.count():
            raise Exception("No organism found with id " + organism_id)
        uniq_name += res.one().common_name
        uniq_name += "  for  "
        res = self.session.query(self.model.analysis).filter_by(analysis_id=analysis_id)
        if not res.count():
            raise Exception("No analysis found with id " + analysis_id)
        uniq_name += res.one().name

        return uniq_name

    def _process_matrix_file(self, file_path, separator):
        # Return a dict with 'biomaterial_list', 'data', and 'feature_list'
        data = []
        feature_list = []
        try:
            with open(file_path) as f:
                reader = csv.reader(f, delimiter=separator)
                # Get headers (biomat list)
                biomaterial_full_list = next(reader)
                # Remove first col (=transcript ids)
                biomaterial_full_list = biomaterial_full_list[1:]
                # Cleanup empty strings
                biomaterial_list = [x for x in biomaterial_full_list if x]
                expected_len = len(biomaterial_list)
                for line in reader:
                    # Get feature name
                    feature_list.append(line.pop(0))
                    float_line = [float(i) for i in line]
                    if not len(float_line) == expected_len:
                        raise Exception("Different number of expressions values and biomaterials for feature '%s'" % feature_list[-1])
                    data.append(float_line)

            if not len(biomaterial_list) == len(set(biomaterial_list)):
                raise Exception("Duplicates found in Biomaterials")

            if not len(feature_list) == len(set(feature_list)):
                raise Exception("Duplicates found in Features")

            return {'biomaterial_list': biomaterial_list, 'feature_list': feature_list, 'data': data}

        except IOError:
            raise Exception("Could not read file at path '%s'" % file_path)

    def _expression_create_biomaterial_structure(self, biomaterial, organism_id, analysis_id, arraydesign_id, contact_id, uniq_name, unit=None):

        # Will create the biomaterial with minimum info and set up DBs
        biomaterial_id = self._create_biomaterial(biomaterial, organism_id, analysis_id=analysis_id, biosourceprovider_id=contact_id)
        uniq_name = biomaterial + uniq_name

        # Create default values for the following
        assay_id = self._create_generic_assay(contact_id, arraydesign_id, uniq_name)
        acquisition_id = self._create_generic_acquisition(assay_id, uniq_name)
        quantification_id = self._create_generic_quantification(acquisition_id, analysis_id, uniq_name)

        # Create generic channel (required for assay_biomaterial table)
        channel_id = self._create_generic_channel()
        self._create_assay_biomaterial(assay_id, biomaterial_id, channel_id)
        # Try to add the unit to cvterm list
        if (unit):
            try:
                cvterm_id = self.ci.get_cvterm_id('unit_of_measure', 'sep')
            except chado.RecordNotFoundError:
                definition = """A unit of measure is a quantity which is a standard of measurement for some dimension.
                                For example, the Meter is a Unit O fMeasure for the dimension of length, as is the Inch.
                                There is no intrinsic property of a UnitOfMeasure that makes it primitive or fundamental;
                                rather, a system of units (e.g. Systeme International Unit) defines a set of orthogonal dimensions and assigns units for each. [ SUMO:unit of measure ]"""
                cvterm_id = self.ci.create_cvterm('unit_of_measure', 'sep', 'sep', term_definition=definition, accession='00056')

            # Add unit
            res = self.session.query(self.model.quantificationprop).filter_by(quantification_id=quantification_id, type_id=cvterm_id, value=unit)
            if not res.count():
                quantificationprop = self.model.quantificationprop()
                quantificationprop.quantification_id = quantification_id
                quantificationprop.type_id = cvterm_id
                quantificationprop.value = unit
                self.session.add(quantificationprop)
                self.session.flush()
                self.session.refresh(quantificationprop)

        return quantification_id

    def _create_expression_element(self, feature_id, arraydesign_id):
        res = self.session.query(self.model.element).filter_by(arraydesign_id=arraydesign_id, feature_id=feature_id)
        element_id = ""
        if res.count():
            element_id = res.one().element_id
        if not element_id:
            element = self.model.element()
            element.arraydesign_id = arraydesign_id
            element.feature_id = feature_id
            self.session.add(element)
            self.session.flush()
            self.session.refresh(element)
            element_id = element.element_id

        return element_id

    def _set_elementresult(self, element_id, quantification_id, signal):
        res = self.session.query(self.model.elementresult).filter_by(element_id=element_id, quantification_id=quantification_id)
        elementresult_id = ""
        if res.count():
            elementresult_id = res.one().elementresult_id
        if not elementresult_id:
            elementresult = self.model.elementresult()
            elementresult.element_id = element_id
            elementresult.quantification_id = quantification_id
            elementresult.signal = signal
            self.session.add(elementresult)
            self.session.flush()
            self.session.refresh(elementresult)
            elementresult_id = elementresult.elementresult_id
        else:
            res.one().signal = signal

        return elementresult_id

    def _manage_feature_expression(self, feature_name, feature_expression_list, quantification_list, organism_id, analysis_id, arraydesign_id, query_type, re_name, skip_missing):
        # Get the feature ID
        feature_id = self._match_feature(feature_name, re_name, query_type, organism_id, skip_missing)

        element_id = self._create_expression_element(feature_id, arraydesign_id)
        self._add_analysis_feature(feature_id, analysis_id)

        # Iterate over quantification (one per biomaterial), and expression value for the selected feature
        for index, quantification_id in enumerate(quantification_list):
            self._set_elementresult(element_id, quantification_id, feature_expression_list[index])

    def _setup():
        """
        if not hasattr(self.model, 'biomaterial'):
            biomaterial_table = Table(
                'biomaterial', self.metadata,
                Column(biomaterial_id, Integer, primary_key=True),
                Column(),
                Column(),
                Column(),
                Column(),
                Column(),
            )
        """
        return None
