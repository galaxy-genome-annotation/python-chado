from . import ChadoTestCase, ci


class LoadTest(ChadoTestCase):

    def test_add_blast(self):

        # Setup testing data
        org = self._create_fake_org()
        an = self._create_fake_an()
        blast_file_path = "./test-data/blastx.xml"
        an_blast = self._create_fake_an('BLAST')
        an_blast_id = an_blast['analysis_id']
        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)
        test = self.ci.load.blast(an_blast_id, blast_file_path, blastdb="swissprot:display", query_type="mRNA")
        assert test['inserted'] == 1, test
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], uniquename='PAC:18136217', analysis_id=an['analysis_id'])
        assert feats != [], "Feature PAC:18136217 was not created"
        feat_id = feats[0]['feature_id']

        res = self.ci.session.query(self.ci.model.analysisfeatureprop) \
            .join(self.ci.model.analysisfeature, self.ci.model.analysisfeature.analysisfeature_id == self.ci.model.analysisfeatureprop.analysisfeature_id) \
            .filter(self.ci.model.analysisfeature.feature_id == feat_id, self.ci.model.analysisfeature.analysis_id == an_blast_id)

        assert res.count(), "No result in analysisfeatureprop table for this feature and analysis"
        assert res.count() == 1, "More than one result in analysisfeatureprop table for this feature and analysis"

        cvterm_id = self.ci.get_cvterm_id("analysis_blast_output_iteration_hits", "tripal")
        res = res.filter(self.ci.model.analysisfeatureprop.type_id == cvterm_id)
        assert res.count(), "Cvterm is not matching analysis_blast_output_iteration_hits"

    def test_add_interpro(self):

        # Setup testing data
        org = self._create_fake_org()
        an = self._create_fake_an()
        interpro_file_path = "./test-data/interproscan.xml"
        an_interpro = self._create_fake_an('INTERPRO')
        an_interpro_id = an_interpro['analysis_id']
        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)
        feats = self.ci.feature.get_features(organism_id=org['organism_id'], uniquename='PAC:18136217', analysis_id=an['analysis_id'])
        assert len(feats) != 0, "Feature PAC:18136217 was not created"
        feat_id = feats[0]['feature_id']

        cv_terms = self.ci.feature.get_feature_cvterms(feat_id)
        assert cv_terms == [], "CV term list is not empty before adding terms!"

        self.ci.load.interpro(an_interpro_id, interpro_file_path, parse_go=True, query_type='mRNA')

        cv_terms = self.ci.feature.get_feature_cvterms(feat_id)
        assert len(cv_terms) == 8, "Number of loaded CV term is not 8!"
        interpro_terms = [d for d in cv_terms if d['db_name'] == "INTERPRO"]
        assert len(interpro_terms) == 5, "Number of INTERPRO cvterms is not 5!"
        go_terms = [d for d in cv_terms if d['db_name'] == "GO"]
        assert len(go_terms) == 3, "Number of GO cvterms is not 3!"

        test_terms = [d for d in cv_terms if d['cvterm_name'] == "Kinase-like_dom_sf"]
        assert len(test_terms) == 1, "Term Kinase-like_dom_sf is not loaded!"
        test_term = test_terms[0]
        assert test_term['dbxref_accession'] == "IPR011009", "Wrong DBXREF accession number"
        assert test_term['db_name'] == "INTERPRO", "Wrong DB name"
        assert test_term['cvterm_definition'] == "Protein kinase-like domain superfamily", "Wrong cvterm definition"

        res = self.ci.session.query(self.ci.model.analysisfeatureprop) \
            .join(self.ci.model.analysisfeature, self.ci.model.analysisfeature.analysisfeature_id == self.ci.model.analysisfeatureprop.analysisfeature_id) \
            .filter(self.ci.model.analysisfeature.feature_id == feat_id, self.ci.model.analysisfeature.analysis_id == an_interpro_id)
        assert res.count(), "No result in analysisfeatureprop table for this feature and analysis"
        assert res.count() == 9, "Number of results not matching in analysisfeatureprop table for this feature and analysis : " + res.count()

        cvterm_id = self.ci.get_cvterm_id("analysis_interpro_xmloutput_hit", "tripal")
        res = res.filter(self.ci.model.analysisfeatureprop.type_id == cvterm_id)
        assert res.count(), "Cvterm analysis_interpro_xmloutput_hit not found in table"

    def setUp(self):

        self.ci = ci
        ci.feature.delete_features()
        ci.organism.delete_organisms()
        ci.analysis.delete_analyses()

        ci.session.commit()

    def tearDown(self):
        ci.feature.delete_features()
        ci.organism.delete_organisms()
        ci.analysis.delete_analyses()

        ci.session.commit()
