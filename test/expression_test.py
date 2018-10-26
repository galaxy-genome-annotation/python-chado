from . import ChadoTestCase, ci


class ExpressionTest(ChadoTestCase):

    def test_add_biomaterial(self):

        # Setup testing data
        biomaterial_name = "test_biomat"
        biomaterial_provider = "Helena Ventseslav"
        biomaterial_description = "This is a description"
        org = self._create_fake_org()

        biomat = self.ci.expression.add_biomaterial(biomaterial_name, org["organism_id"], description=biomaterial_description, biomaterial_provider=biomaterial_provider)

        assert biomat, "Issue : Biomaterial not created"
        assert biomat["name"] == biomaterial_name, "Biomaterial name issue : " + biomat["name"]
        assert biomat["provider_id"], "Bioprovider not set! "
        assert biomat["description"] == biomaterial_description, "Description issue : " + biomat["description"]

    def test_add_expression(self):

        # Setup testing data
        # Expression file
        expression_file_path = "./test-data/expression.matrix"

        org = self._create_fake_org()
        an = self._create_fake_an()

        # Feature file (fasta)
        feature_file_path = "./test-data/proteins.fa"
        self.ci.feature.load_fasta(fasta=feature_file_path, analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.expression.add_expression(org['organism_id'], an['analysis_id'], expression_file_path, separator="\t")

        biomat_list = self.ci.expression.get_biomaterials()
        assert len(biomat_list) == 8, "Unexpected number of biomaterials created"

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
