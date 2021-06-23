from . import ChadoTestCase, ci


class ExpressionTest(ChadoTestCase):

    def test_add_biomaterial(self):

        # Setup testing data
        biomaterial_name = "test_biomat"
        biomaterial_provider = "Helena Ventseslav"
        biomaterial_description = "This is a description"
        biosample_accession = 'biosample-foobar'
        sra_accession = "sra-bar"
        bioproject_accession = 'bioproject-foo'
        attrs = {
            'foo': 'bar',
            'oof': 'arb'
        }
        org = self._create_fake_org()

        biomat_ret = self.ci.expression.add_biomaterial(biomaterial_name, org["organism_id"], description=biomaterial_description, biomaterial_provider=biomaterial_provider, biosample_accession=biosample_accession, sra_accession=sra_accession, bioproject_accession=bioproject_accession, attributes=attrs)

        biomat = self.ci.expression.get_biomaterials(biomaterial_name=biomaterial_name)[0]

        assert biomat_ret == biomat

        assert biomat, "Issue : Biomaterial not created"
        assert biomat["name"] == biomaterial_name, "Biomaterial name issue : " + biomat["name"]
        assert biomat["provider_id"], "Bioprovider not set! "
        assert biomat["description"] == biomaterial_description, "Description issue : " + biomat["description"]

        res = self.ci.session.query(self.ci.model.biomaterialprop).filter_by(biomaterial_id=biomat['biomaterial_id'])

        props = {prop.type_id: prop.value for prop in res}
        for at in attrs:
            type_id = self.ci.get_cvterm_id(at, 'biomaterial_property')
            assert type_id in props

            assert props[type_id] == attrs[at]

        dbs = self.ci.session.query(self.ci.model.db.db_id, self.ci.model.db.name, self.ci.model.db.description) \
            .filter((self.ci.model.db.name == 'NCBI SRA') | (self.ci.model.db.name == 'NCBI BioProject') | (self.ci.model.db.name == 'NCBI BioSample'))

        dbs = {db.name: db.db_id for db in dbs}

        biomat = self.ci.session.query(self.ci.model.biomaterial).filter_by(biomaterial_id=biomat['biomaterial_id']).one()

        xrefs = {dbx.dbxref.db_id: dbx.dbxref.accession for dbx in biomat.biomaterial_dbxref_collection}

        assert len(dbs) == 3
        assert len(xrefs) == 3

        assert xrefs[dbs['NCBI SRA']] == "sra-bar"
        assert xrefs[dbs['NCBI BioProject']] == "bioproject-foo"
        assert xrefs[dbs['NCBI BioSample']] == "biosample-foobar"

    def test_add_expression(self):

        # Setup testing data
        # Expression file
        expression_file_path = "./test-data/expression.matrix"

        org = self._create_fake_org()
        an = self._create_fake_an()

        # Feature file (fasta)
        feature_file_path = "./test-data/proteins.fa"
        self.ci.feature.load_fasta(fasta=feature_file_path, analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='mRNA')
        self.ci.expression.add_expression(org['organism_id'], an['analysis_id'], expression_file_path, separator="\t")

        biomat_list = self.ci.expression.get_biomaterials()
        assert len(biomat_list) == 8, "Unexpected number of biomaterials created"

        elements = self.ci.session.query(self.ci.model.element.arraydesign_id, self.ci.model.elementresult.quantification_id, self.ci.model.elementresult.signal) \
            .join(self.ci.model.elementresult, self.ci.model.element.element_id == self.ci.model.elementresult.element_id)

        assert elements.count() == 24

    def test_add_expression_title(self):

        # Setup testing data
        # Expression file
        expression_file_path = "./test-data/expression_title.matrix"

        org = self._create_fake_org()
        an = self._create_fake_an()

        # Feature file (fasta)
        feature_file_path = "./test-data/proteins.fa"
        self.ci.feature.load_fasta(fasta=feature_file_path, analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='mRNA')
        self.ci.expression.add_expression(org['organism_id'], an['analysis_id'], expression_file_path, separator="\t")

        biomat_list = self.ci.expression.get_biomaterials()
        assert len(biomat_list) == 8, "Unexpected number of biomaterials created"

        elements = self.ci.session.query(self.ci.model.element.arraydesign_id, self.ci.model.elementresult.quantification_id, self.ci.model.elementresult.signal) \
            .join(self.ci.model.elementresult, self.ci.model.element.element_id == self.ci.model.elementresult.element_id)

        assert elements.count() == 24

    def setUp(self):

        self.ci = ci
        ci.feature.delete_features()
        ci.organism.delete_organisms()
        ci.analysis.delete_analyses()
        ci.expression.delete_all_biomaterials(confirm=True)

        ci.session.commit()

    def tearDown(self):
        ci.feature.delete_features()
        ci.organism.delete_organisms()
        ci.analysis.delete_analyses()
        ci.expression.delete_all_biomaterials(confirm=True)

        ci.session.commit()
