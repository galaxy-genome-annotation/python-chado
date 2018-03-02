import unittest

from . import ci


class FeaturePropTest(unittest.TestCase):

    def test_load_props(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        # there's a contig loaded by fasta and a supercontig in gff
        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        self.ci.feature.load_featureprops(tab_file="./test-data/props.tsv", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], prop_type='Foobar')

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015633m.g") \
            .one()

        # Check gene prop
        assert len(gene_f.featureprop_collection) == 1, "prop correctly loaded"
        assert gene_f.featureprop_collection[0].value == 'It\'s so cool', "prop correctly loaded"
        assert gene_f.featureprop_collection[0].rank == 0, "prop correctly loaded"
        assert gene_f.featureprop_collection[0].cvterm.name == "Foobar", "prop correctly loaded"
        assert gene_f.featureprop_collection[0].cvterm.cv.name == "feature_property", "prop correctly loaded"

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
