import unittest

from nose.tools import raises

from . import ci


class FastaTest(unittest.TestCase):

    def test_load_fasta_simple(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
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

    def test_load_fasta_seqtype(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        f = feats[0]
        contigterm = self.ci.get_cvterm_id('supercontig', 'sequence')
        assert f['type_id'] == contigterm, "correct sequence_type"

    def test_load_fasta_rename(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_name='([a-zA-Z0-9]+)\|[a-zA-Z0-9_]+')
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123')

        f = feats[0]
        assert f['name'] == 'Q02123', "fasta re_name"
        assert f['uniquename'] == 'Q02123|VNBP_POPMV', "fasta re_uniquename"

    def test_load_fasta_reuname(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_uniquename='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)')
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        f = feats[0]
        assert f['name'] == 'Q02123|VNBP_POPMV', "fasta re_name"
        assert f['uniquename'] == 'VNBP_POPMV', "fasta re_uniquename"

    def test_load_fasta_wrong_rename(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_name='([a-zA-Z]+)\|[a-zA-Z0-9_]+', re_uniquename='([a-zA-Z]+)\|[a-zA-Z0-9_]+')
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        f = feats[0]
        assert f['name'] == 'Q02123|VNBP_POPMV', "fasta wrong re_name"
        assert f['uniquename'] == 'Q02123|VNBP_POPMV', "fasta wrong re_uniquename"

    def test_load_fasta_update(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.load_fasta(fasta="./test-data/proteins_alt.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], update=True)
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        f = feats[0]
        assert f['residues'] == 'MVXXXDPR', "fasta features properly updated"
        assert f['seqlen'] == 8, "fasta features properly updated"
        assert f['md5checksum'] == '2eb1f62c8a6af8662f4ad1e155543779', "fasta features properly updated"

    @raises(Exception)
    def test_load_fasta_update_fail(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.load_fasta(fasta="./test-data/proteins_alt.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])

    @raises(Exception)
    def test_load_fasta_duplicate(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins_dup.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])

    def test_load_fasta_update_match_on_name_fail(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_uniquename='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)')
        self.ci.feature.load_fasta(fasta="./test-data/proteins_alt.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], update=True)
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], uniquename='Q02123|VNBP_POPMV')

        f = feats[0]
        assert f['residues'] == 'MVXXXDPR', "fasta features properly updated"
        assert f['seqlen'] == 8, "fasta features properly updated"
        assert f['md5checksum'] == '2eb1f62c8a6af8662f4ad1e155543779', "fasta features properly updated"

        feats = self.ci.feature.get_features(organism_id=org['organism_id'], uniquename='VNBP_POPMV')
        f = feats[0]
        assert f['seqlen'] == 123, "fasta features properly updated"

    def test_load_fasta_update_match_on_name(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_uniquename='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)')
        self.ci.feature.load_fasta(fasta="./test-data/proteins_alt.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], update=True, match_on_name=True)
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        f = feats[0]
        assert f['residues'] == 'MVXXXDPR', "fasta features properly updated"
        assert f['seqlen'] == 8, "fasta features properly updated"
        assert f['md5checksum'] == '2eb1f62c8a6af8662f4ad1e155543779', "fasta features properly updated"

    def test_load_fasta_db(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        # get a test db
        db = self.ci.session.query(self.ci.model.db).filter_by(name='null')[0]

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], db=db.db_id, re_db_accession='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)')
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        f = feats[0]
        assert f['dbxref_id'] > 0, "dbxref created"

        dbxref = self.ci.session.query(self.ci.model.dbxref).filter_by(dbxref_id=f['dbxref_id'])[0]
        assert dbxref.accession == 'VNBP_POPMV', "dbxref created correctly"

    def test_load_fasta_relation(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_name='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)', re_uniquename='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)', sequence_type='supercontig')
        self.ci.feature.load_fasta(fasta="./test-data/proteins_alt.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], re_name='([a-zA-Z0-9]+)\|[a-zA-Z0-9_]+', re_uniquename='([a-zA-Z0-9]+)\|[a-zA-Z0-9_]+', rel_type='part_of', re_parent='[a-zA-Z0-9]+\|([a-zA-Z0-9_]+)', parent_type='supercontig')

        child = self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123')[0]
        parent = self.ci.feature.get_features(organism_id=org['organism_id'], name='VNBP_POPMV')[0]

        rship = self.ci.session.query(self.ci.model.feature_relationship).filter_by(subject_id=child['feature_id']).all()[0]
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')

        assert rship.subject_id == child['feature_id'], "fasta features properly loaded with relationship"
        assert rship.object_id == parent['feature_id'], "fasta features properly loaded with relationship"
        assert rship.type_id == partofterm, "fasta features properly loaded with relationship"

    def test_load_fasta_remove_by_org(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        assert len(self.ci.feature.get_features()) == 21, "features are loaded"

        self.ci.organism.delete_organisms(organism_id=org['organism_id'])

        assert len(self.ci.feature.get_features()) == 0, "features removed when removing organism"

    def test_load_fasta_remove_by_an(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')

        assert len(self.ci.feature.get_features()) == 21, "features are loaded"

        self.ci.analysis.delete_analyses(analysis_id=an['analysis_id'])

        assert len(self.ci.feature.get_features()) == 0, "features removed when removing analysis"

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
