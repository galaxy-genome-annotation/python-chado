import unittest

from . import ci


class FeatureTest(unittest.TestCase):

    def test_get_features(self):

        org = self._create_fake_org()
        org2 = self._create_fake_org('another')
        an = self._create_fake_an()
        an2 = self._create_fake_an('another')

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        feats = self.ci.feature.get_features(organism_id=org['organism_id'])
        assert len(feats) == 21, "fasta features properly created once"

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an2['analysis_id'], organism_id=org2['organism_id'])
        feats = self.ci.feature.get_features()
        assert len(feats) == 42, "fasta features properly created twice"
        feats = self.ci.feature.get_features(organism_id=org['organism_id'])
        assert len(feats) == 21, "fasta features properly fetched by organism"

        feats = self.ci.feature.get_features(analysis_id=an['analysis_id'])
        assert len(feats) == 21, "fasta features properly fetched by analysis_id"

        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')
        f = feats[0]
        contigterm = self.ci.get_cvterm_id('contig', 'sequence')
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

    def test_delete_features(self):

        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])

        deleted_feats = self.ci.feature.delete_features(uniquename='Q02123|VNBP_POPMV')
        feats = self.ci.feature.get_features()
        assert len(feats) == 20, "single fasta features properly deleted"
        assert deleted_feats == 1, "single fasta features properly deleted"

        deleted_feats = self.ci.feature.delete_features(analysis_id=an['analysis_id'])
        feats = self.ci.feature.get_features()
        assert len(feats) == 0, "fasta features properly deleted using analysis_id"
        assert deleted_feats == 20, "fasta features properly deleted using analysis_id"

    def _create_fake_org(self, uniqid=''):
        genus = "Testus" + uniqid
        common = "Testorg" + uniqid
        abbr = "Ttesta" + uniqid
        species = "testa" + uniqid
        comment = "A test org sta" + uniqid

        return self.ci.organism.add_organism(genus=genus, common=common, abbr=abbr, species=species, comment=comment)

    def _create_fake_an(self, uniqid=''):
        name = "analysis x" + uniqid
        program = "Magic" + uniqid
        programversion = "1.0" + uniqid
        algorithm = "mind" + uniqid
        sourcename = "src" + uniqid
        sourceversion = "2.1beta" + uniqid
        sourceuri = "http://example.org/"
        date_executed = "2018-02-03"

        return self.ci.analysis.add_analysis(name=name, program=program, programversion=programversion, algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri, date_executed=date_executed)

    def _del_dbxref(self):
        self.ci.session.query(self.ci.model.dbxref).filter(
            self.ci.model.dbxref.db_id == 1,
            (self.ci.model.dbxref.accession.like('VNBP%') | self.ci.model.dbxref.accession.like('%VIRU'))
        ).delete(synchronize_session='fetch')

    def setUp(self):
        self.ci = ci
        self.ci.organism.delete_organisms()
        self.ci.analysis.delete_analyses()
        self.ci.feature.delete_features()

        # Make sure dbxref are deleted too
        self._del_dbxref()

        self.ci.session.commit()

    def tearDown(self):
        self.ci.organism.delete_organisms()
        self.ci.analysis.delete_analyses()
        self.ci.feature.delete_features()

        # Make sure dbxref are deleted too
        self._del_dbxref()

        self.ci.session.commit()
