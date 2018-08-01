from . import ChadoTestCase, ci


class FeaturePropTest(ChadoTestCase):

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
