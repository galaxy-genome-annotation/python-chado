from nose.tools import raises

from . import ChadoTestCase, ci_reflect_tripal


class GoTest(ChadoTestCase):

    def test_load_go_simple(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='polypeptide')

        an_go = self._create_fake_an('GO')
        self.ci.load.go(input="./test-data/go.gaf", analysis_id=an_go['analysis_id'], organism_id=org['organism_id'])

        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV', analysis_id=an['analysis_id'])

        f = feats[0]
        contigterm = self.ci.get_cvterm_id('polypeptide', 'sequence')
        assert f['dbxref_id'] is None, "fasta features properly created"
        assert f['organism_id'] == org['organism_id'], "fasta features properly created"
        assert f['analysis_id'] == an['analysis_id'], "fasta features properly created"
        assert f['name'] == 'Q02123|VNBP_POPMV', "fasta features properly created"
        assert f['uniquename'] == 'Q02123|VNBP_POPMV', "fasta features properly created"
        assert f['residues'] == 'MVNMRKVLALMQVFRERYDHKCDFNFCDIAVSIVCRSELDFINEPGLSNYAKRRRARRLGRCVRCFRVNPGFYFTKRCDGITCVPGISWNYDVEDYIKRGRVTGDRETPSTFHGYGYPVGHKT', "fasta features properly created"
        assert f['seqlen'] == 123, "fasta features properly created"
        assert f['md5checksum'] == 'a10c50557506954bf61efe1f997aa8d3', "fasta features properly created"
        assert f['type_id'] == contigterm, "fasta features properly created"
        assert f['is_analysis'] is False, "fasta features properly created"
        assert f['is_obsolete'] is False, "fasta features properly created"

        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV', analysis_id=an_go['analysis_id'])

        assert len(feats) == 1, "feature associated with go analysis too"

        cvterms = self.ci.feature.get_feature_cvterms(feats[0]['feature_id'])
        assert len(cvterms) == 1, "1 cvterm association"

        term = cvterms[0]
        assert term['db_name'] == "GO", "associated to correct cvterm"
        assert term['dbxref_accession'] == "0005852", "associated to correct cvterm"
        assert term['cv_name'] == "cellular_component", "associated to correct cvterm"

        ans = self.ci.feature.get_feature_analyses(feats[0]['feature_id'])
        assert len(ans) == 2, "2 analysis association"

        ansass = ans[1]['analysisfeatureprop'][0]
        print(ans[1]['analysisfeatureprop'])
        assert ansass['value'] == "GO:0005852", "analysisfeatureprop associated to correct cvterm"
        assert ansass['db_name'] == "GO", "analysisfeatureprop associated to correct cvterm"
        assert ansass['dbxref_accession'] == "0005852", "analysisfeatureprop associated to correct cvterm"
        assert ansass['cv_name'] == "cellular_component", "analysisfeatureprop associated to correct cvterm"

    def test_load_go_tsv_re(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='polypeptide')

        an_go = self._create_fake_an('GO')
        self.ci.load.go(input="./test-data/go.tsv", analysis_id=an_go['analysis_id'], organism_id=org['organism_id'], name_column=1, go_column=2, match_on_name=True, re_name='(.*)FOOBAR')

        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV', analysis_id=an['analysis_id'])

        f = feats[0]
        contigterm = self.ci.get_cvterm_id('polypeptide', 'sequence')
        assert f['dbxref_id'] is None, "fasta features properly created"
        assert f['organism_id'] == org['organism_id'], "fasta features properly created"
        assert f['analysis_id'] == an['analysis_id'], "fasta features properly created"
        assert f['name'] == 'Q02123|VNBP_POPMV', "fasta features properly created"
        assert f['uniquename'] == 'Q02123|VNBP_POPMV', "fasta features properly created"
        assert f['residues'] == 'MVNMRKVLALMQVFRERYDHKCDFNFCDIAVSIVCRSELDFINEPGLSNYAKRRRARRLGRCVRCFRVNPGFYFTKRCDGITCVPGISWNYDVEDYIKRGRVTGDRETPSTFHGYGYPVGHKT', "fasta features properly created"
        assert f['seqlen'] == 123, "fasta features properly created"
        assert f['md5checksum'] == 'a10c50557506954bf61efe1f997aa8d3', "fasta features properly created"
        assert f['type_id'] == contigterm, "fasta features properly created"
        assert f['is_analysis'] is False, "fasta features properly created"
        assert f['is_obsolete'] is False, "fasta features properly created"

        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV', analysis_id=an_go['analysis_id'])

        assert len(feats) == 1, "feature associated with go analysis too"

        cvterms = self.ci.feature.get_feature_cvterms(feats[0]['feature_id'])
        assert len(cvterms) == 1, "1 cvterm association"

        term = cvterms[0]
        assert term['db_name'] == "GO", "associated to correct cvterm"
        assert term['dbxref_accession'] == "0005852", "associated to correct cvterm"
        assert term['cv_name'] == "cellular_component", "associated to correct cvterm"

        ans = self.ci.feature.get_feature_analyses(feats[0]['feature_id'])
        assert len(ans) == 2, "2 analysis association"

        ansass = ans[1]['analysisfeatureprop'][0]
        print(ans[1]['analysisfeatureprop'])
        assert ansass['value'] == "GO:0005852", "analysisfeatureprop associated to correct cvterm"
        assert ansass['db_name'] == "GO", "analysisfeatureprop associated to correct cvterm"
        assert ansass['dbxref_accession'] == "0005852", "analysisfeatureprop associated to correct cvterm"
        assert ansass['cv_name'] == "cellular_component", "analysisfeatureprop associated to correct cvterm"

    @raises(Exception)
    def test_load_go_bad_feat(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='polypeptide')

        an_go = self._create_fake_an('GO')
        self.ci.load.go(input="./test-data/go.tsv", analysis_id=an_go['analysis_id'], organism_id=org['organism_id'], name_column=1, go_column=2)

    @raises(Exception)
    def test_load_go_bad_query(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])

        an_go = self._create_fake_an('GO')
        self.ci.load.go(input="./test-data/go.gaf", analysis_id=an_go['analysis_id'], organism_id=org['organism_id'], query_type='foobar')

    def setUp(self):
        self.ci = ci_reflect_tripal
        self.ci.organism.delete_organisms()
        self.ci.analysis.delete_analyses()
        self.ci.feature.delete_features()

        self.ci.session.commit()

    def tearDown(self):
        self.ci.organism.delete_organisms()
        self.ci.analysis.delete_analyses()
        self.ci.feature.delete_features()

        self.ci.session.commit()
