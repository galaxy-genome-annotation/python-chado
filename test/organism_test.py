from . import ChadoTestCase, ci, ci_no_reflect


class OrganismTest(ChadoTestCase):

    def test_add_orgs(self):

        genus = "Testus"
        common = "Testorg"
        abbr = "Ttesta"
        species = "testa"
        comment = "A test org"

        org = self.ci.organism.add_organism(genus=genus, common=common, abbr=abbr, species=species, comment=comment)

        assert org["organism_id"] > 0, "org properly created"

        org = self.ci.organism.get_organisms(abbr=abbr)[0]

        assert org["genus"] == genus, "org properly created"
        assert org["common_name"] == common, "org properly created"
        assert org["abbreviation"] == abbr, "org properly created"
        assert org["species"] == species, "org properly created"
        assert org["comment"] == comment, "org properly created"

        self.ci.organism.delete_organisms(genus=genus)

        orgs = self.ci.organism.get_organisms(genus=genus)

        assert len(orgs) == 0, "orgs properly deleted"

    def test_delete_cascade(self):
        org = self._create_fake_org()

        ci.feature.load_fasta(fasta="./test-data/proteins.fa", organism_id=org['organism_id'])

        ci_no_reflect.organism.delete_organisms(genus=org["genus"])

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

        ci_no_reflect.organism.delete_organisms(genus=org["genus"])

        feats = len(ci.feature.get_features(organism_id=org['organism_id']))
        assert(feats == 0)

        props = ci.session.query(ci.model.featureprop).count()
        assert(props == 0)

        locs = ci.session.query(ci.model.featureloc).count()
        assert(locs == 0)

        rels = ci.session.query(ci.model.feature_relationship).count()
        assert(rels == 0)

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
