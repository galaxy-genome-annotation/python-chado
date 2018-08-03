from . import ChadoTestCase, ci, ci_no_reflect


class AnalysisTest(ChadoTestCase):

    def test_add_analysis(self):

        name = "analysis x"
        program = "Magic"
        programversion = "1.0"
        algorithm = "mind"
        sourcename = "src"
        sourceversion = "2.1beta"
        sourceuri = "http://example.org/"
        description = "Bla bla bla"
        date_executed = "2018-02-03"

        ana = ci_no_reflect.analysis.add_analysis(name=name, program=program, programversion=programversion, algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri, description=description, date_executed=date_executed)

        assert ana["analysis_id"] > 0, "ana properly created"

        ana = ci_no_reflect.analysis.get_analyses(name=name)[0]

        assert ana["name"] == name, "analysis properly created"
        assert ana["program"] == program, "analysis properly created"
        assert ana["programversion"] == programversion, "analysis properly created"
        assert ana["algorithm"] == algorithm, "analysis properly created"
        assert ana["sourcename"] == sourcename, "analysis properly created"
        assert ana["sourceversion"] == sourceversion, "analysis properly created"
        assert ana["sourceuri"] == sourceuri, "analysis properly created"
        assert ana["description"] == description, "analysis properly created"
        assert ana["timeexecuted"] == '2018-02-03 00:00:00', "analysis properly created"

        ci_no_reflect.analysis.delete_analyses(name=name)

        ans = ci_no_reflect.analysis.get_analyses(name=name)

        assert len(ans) == 0, "analysis properly deleted"

    def test_delete_cascade(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])

        ci_no_reflect.analysis.delete_analyses(analysis_id=an['analysis_id'])

        feats = ci.feature.get_features(organism_id=org['organism_id'])

        assert(len(feats) == 0)

    def test_delete_cascade_featureprops(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        props = ci.session.query(ci.model.featureprop).count()
        assert(props == 0)

        ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        props = ci.session.query(ci.model.featureprop).count()
        assert(props == 174)

        locs = ci.session.query(ci.model.featureloc).count()
        assert(locs == 182)

        rels = ci.session.query(ci.model.feature_relationship).count()
        assert(rels == 180)

        ci_no_reflect.analysis.delete_analyses(analysis_id=an_gff['analysis_id'])

        feats = len(ci.feature.get_features(organism_id=org['organism_id']))
        assert(feats == 1)  # The genome in first analysis

        props = ci.session.query(ci.model.featureprop).count()
        assert(props == 0)

        locs = ci.session.query(ci.model.featureloc).count()
        assert(locs == 1)

        rels = ci.session.query(ci.model.feature_relationship).count()
        assert(rels == 0)

    def test_delete_cascade_multiple(self):
        # When deleting analysis A don't delete features that are part of analysis A and B
        org = self._create_fake_org()
        an = self._create_fake_an()
        an2 = self._create_fake_an("second")

        ci.feature.load_fasta(fasta="./test-data/proteins.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])

        feat = ci.feature.get_features(organism_id=org['organism_id'], name='Q02123|VNBP_POPMV')[0]

        af = ci.model.analysisfeature()
        af.analysis_id = an2['analysis_id']
        af.feature_id = feat['feature_id']
        ci.session.add(af)
        ci.session.commit()

        ci_no_reflect.analysis.delete_analyses(analysis_id=an['analysis_id'])

        afs = ci.session.query(ci.model.analysisfeature).filter_by(feature_id=feat['feature_id'], analysis_id=an['analysis_id'])
        assert(afs.count() == 0)

        afs = ci.session.query(ci.model.analysisfeature).filter_by(feature_id=feat['feature_id'], analysis_id=an2['analysis_id'])
        assert(afs.count() == 1)

    def setUp(self):
        ci.feature.delete_features()
        ci_no_reflect.organism.delete_organisms()
        ci_no_reflect.analysis.delete_analyses()
        ci.session.query(ci.model.featureprop).delete()
        ci.session.query(ci.model.featureloc).delete()
        ci.session.query(ci.model.feature_relationship).delete()

        ci.session.commit()

        self.ci = ci_no_reflect

    def tearDown(self):
        ci.feature.delete_features()
        ci_no_reflect.organism.delete_organisms()
        ci_no_reflect.analysis.delete_analyses()
        ci.session.query(ci.model.featureprop).delete()
        ci.session.query(ci.model.featureloc).delete()
        ci.session.query(ci.model.feature_relationship).delete()

        ci.session.commit()
