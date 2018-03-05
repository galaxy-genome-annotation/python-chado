import unittest

from nose.tools import raises

from . import ci


class GFFTest(unittest.TestCase):

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

    def test_load_gff(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015632m.g") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        # Check gene feature
        assert gene_f.dbxref_id is None, "gff>gene loaded correctly"
        assert gene_f.organism_id == org['organism_id'], "gff>gene loaded correctly"
        assert gene_f.name == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.uniquename == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.residues is None, "gff>gene loaded correctly"
        assert gene_f.seqlen is None, "gff>gene loaded correctly"
        assert gene_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>gene loaded correctly"
        assert gene_f.type_id == geneterm, "gff>gene loaded correctly"
        assert gene_f.is_analysis is False, "gff>gene loaded correctly"
        assert gene_f.is_obsolete is False, "gff>gene loaded correctly"

        # Check gene loc
        assert len(gene_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert gene_f.featureloc_collection[0].fmin == 4058459, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].fmax == 4062210, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmin_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmax_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].strand == 1, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].phase is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].residue_info is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].locgroup == 0, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].rank == 0, "gff>gene located correctly"

        src_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=gene_f.featureloc_collection[0].srcfeature_id) \
            .one()

        assert src_f.uniquename == "scaffold00001", "gff>gene loaded correctly"
        scaff1_id = src_f.feature_id

        # Check gene aliases
        exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
        syns = {synf.synonym.name: synf.synonym.type_id for synf in gene_f.feature_synonym_collection}

        assert len(syns) == 2, "gff>gene aliases loaded correctly"
        assert 'some-synonym' in syns, "gff>gene aliases loaded correctly"
        assert 'another synonym' in syns, "gff>gene aliases loaded correctly"
        assert syns['some-synonym'] == exactterm, "gff>gene aliases loaded correctly"
        assert syns['another synonym'] == exactterm, "gff>gene aliases loaded correctly"

        # Check gene dbxref
        dbs = self.ci.session.query(self.ci.model.db.db_id, self.ci.model.db.name, self.ci.model.db.description) \
            .filter((self.ci.model.db.name == 'GO') | (self.ci.model.db.name == 'FOOBAR') | (self.ci.model.db.name == 'FOOBARXX') | (self.ci.model.db.name == 'GFF_source'))
        for db in dbs:
            if db.name == "FOOBAR":
                assert db.description == "Added automatically by the GFF loader", "gff>gene dbxrefs db loaded correctly"

        dbs = {db.name: db.db_id for db in dbs}

        assert len(dbs) == 4, "gff>gene dbxrefs db loaded correctly"

        xrefs = {dbx.dbxref.accession: dbx.dbxref.db_id for dbx in gene_f.feature_dbxref_collection}

        assert len(xrefs) == 3, "gff>gene dbxrefs loaded correctly"
        assert '0061611' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert '6528B' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert 'phytozome6' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert xrefs['0061611'] == dbs['GO'], "gff>gene dbxrefs loaded correctly"
        assert xrefs['6528B'] == dbs['FOOBAR'], "gff>gene dbxrefs loaded correctly"
        assert xrefs['phytozome6'] == dbs['GFF_source'], "gff>gene dbxrefs loaded correctly"

        # Check gene featureprop
        expected = [
            'Gap___BLABLA___0',
            'Gap___BLOBLO___1',
            'Note___that\'s fantastic___0',
            'Note___really___1',
            'Poutrelle___test___1',
            'Poutrelle___lapinou___0',
        ]

        assert len(gene_f.featureprop_collection) == 6, "gff>gene loaded correctly"

        for prop in gene_f.featureprop_collection:
            assert prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank) in expected, "gff>gene loaded correctly"
            expected.remove(prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank))

        # Check mrna
        rna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136219") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        rnaterm = self.ci.get_cvterm_id('mRNA', 'sequence')

        # Check mRNA feature
        assert rna_f.dbxref_id is None, "gff>mRNA loaded correctly"
        assert rna_f.organism_id == org['organism_id'], "gff>mRNA loaded correctly"
        assert rna_f.name == "orange1.1g015615m", "gff>mRNA loaded correctly"
        assert rna_f.uniquename == "PAC:18136219", "gff>mRNA loaded correctly"
        assert rna_f.residues is None, "gff>mRNA loaded correctly"
        assert rna_f.seqlen is None, "gff>mRNA loaded correctly"
        assert rna_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>mRNA loaded correctly"
        assert rna_f.type_id == rnaterm, "gff>mRNA loaded correctly"
        assert rna_f.is_analysis is False, "gff>mRNA loaded correctly"
        assert rna_f.is_obsolete is False, "gff>mRNA loaded correctly"

        # Check mRNA loc
        assert len(rna_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert rna_f.featureloc_collection[0].fmin == 4058759, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].fmax == 4062210, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].is_fmin_partial is False, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].is_fmax_partial is False, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].strand == 1, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].phase is None, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].residue_info is None, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].locgroup == 0, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].rank == 0, "gff>mRNA located correctly"

        assert scaff1_id == rna_f.featureloc_collection[0].srcfeature_id, "gff>mRNA loaded correctly"

        # Check mRNA aliases
        exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
        syns = {synf.synonym.name: synf.synonym.type_id for synf in rna_f.feature_synonym_collection}

        assert len(syns) == 2, "gff>mRNA aliases loaded correctly"
        assert 'some-synonym' in syns, "gff>mRNA aliases loaded correctly"
        assert 'another synonym' in syns, "gff>mRNA aliases loaded correctly"
        assert syns['some-synonym'] == exactterm, "gff>mRNA aliases loaded correctly"
        assert syns['another synonym'] == exactterm, "gff>mRNA aliases loaded correctly"

        # Check mRNA dbxref
        xrefs = {dbx.dbxref.accession: dbx.dbxref.db_id for dbx in rna_f.feature_dbxref_collection}

        assert len(xrefs) == 3, "gff>mRNA dbxrefs loaded correctly"
        assert '0061621' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert '6528A' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert 'phytozome6' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['0061621'] == dbs['GO'], "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['6528A'] == dbs['FOOBARXX'], "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['phytozome6'] == dbs['GFF_source'], "gff>mRNA dbxrefs loaded correctly"

        # Check relationships
        parents = {x.object_id: x.type_id for x in rna_f.subject_in_relationships}
        assert len(parents) == 1, "mRNA relationships"
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')
        assert rna_f.subject_in_relationships[0].type_id == partofterm, "mRNA relationships"

        derivesfromterm = self.ci.get_cvterm_id('derives_from', 'sequence')
        peps = [x.subject_id for x in rna_f.object_in_relationships if x.type_id == derivesfromterm]
        assert len(peps) == 1, "mRNA relationships, single peptide"

        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=peps[0]) \
            .one()

        # Check pep feature
        pepterm = self.ci.get_cvterm_id('polypeptide', 'sequence')
        assert pep_f.dbxref_id is None, "gff>pep loaded correctly"
        assert pep_f.organism_id == org['organism_id'], "gff>pep loaded correctly"
        assert pep_f.name == "orange1.1g015615m", "gff>pep loaded correctly"
        assert pep_f.uniquename == "PAC:18136219-protein", "gff>pep loaded correctly"
        assert pep_f.residues is None, "gff>pep loaded correctly"
        assert pep_f.seqlen is None, "gff>pep loaded correctly"
        assert pep_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>pep loaded correctly"
        assert pep_f.type_id == pepterm, "gff>pep loaded correctly"
        assert pep_f.is_analysis is False, "gff>pep loaded correctly"
        assert pep_f.is_obsolete is False, "gff>pep loaded correctly"

        # Check pep loc
        assert len(pep_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].fmin == 4059234, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].fmax == 4061905, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].is_fmin_partial is False, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].is_fmax_partial is False, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].strand == 1, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].phase is None, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].residue_info is None, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].locgroup == 0, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].rank == 0, "gff>pep located correctly"

        assert scaff1_id == pep_f.featureloc_collection[0].srcfeature_id, "gff>pep loaded correctly"

        children = {x.subject_id: x for x in rna_f.object_in_relationships if x.type_id != derivesfromterm}
        assert len(children) == 15, "mRNA relationships, single peptide"

        cdsterm = self.ci.get_cvterm_id('CDS', 'sequence')
        exonterm = self.ci.get_cvterm_id('exon', 'sequence')
        utr3term = self.ci.get_cvterm_id('three_prime_UTR', 'sequence')
        utr5term = self.ci.get_cvterm_id('five_prime_UTR', 'sequence')
        for c in children:
            assert children[c].type_id == partofterm, "subsubfeatures"

            if children[c].subject.type_id == utr3term:
                subsub_f = children[c].subject

            assert children[c].subject.type_id in (cdsterm, exonterm, utr3term, utr5term), "subsubfeatures"

        # Check a subsubfeature
        assert subsub_f.dbxref_id is None, "gff>utr loaded correctly"
        assert subsub_f.organism_id == org['organism_id'], "gff>utr loaded correctly"
        assert subsub_f.name.endswith("-three_prime_UTR-scaffold00001:4061905..4062210"), "gff>utr loaded correctly"
        assert subsub_f.uniquename.endswith("-three_prime_UTR-scaffold00001:4061905..4062210"), "gff>utr loaded correctly"
        assert subsub_f.residues is None, "gff>utr loaded correctly"
        assert subsub_f.seqlen is None, "gff>utr loaded correctly"
        assert subsub_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>utr loaded correctly"
        assert subsub_f.type_id == utr3term, "gff>utr loaded correctly"
        assert subsub_f.is_analysis is False, "gff>utr loaded correctly"
        assert subsub_f.is_obsolete is False, "gff>utr loaded correctly"

        assert len(subsub_f.featureloc_collection) == 1, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].fmin == 4061905, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].fmax == 4062210, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].is_fmin_partial is False, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].is_fmax_partial is False, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].strand == 1, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].phase is None, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].residue_info is None, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].locgroup == 0, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].rank == 0, "gff>utr located correctly"

        # Check utr with 2 parents
        confused_child_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename='an_utr_with_two_parents') \
            .all()

        assert len(confused_child_f) == 1, "1 utr with 2 parents"

        confused_rels = confused_child_f[0].subject_in_relationships

        assert len(confused_rels) == 2, "1 utr with 2 parents"

        for r in confused_rels:
            assert (r.object.uniquename == 'PAC:18136239') or (r.object.uniquename == 'PAC:18136238'), "1 utr with 2 parents"

        # Check Derives_from
        derivesfrom = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename='some_special_cds') \
            .all()

        assert len(derivesfrom) == 1, "derives_from"

        derivesfrom_rels = derivesfrom[0].subject_in_relationships

        assert len(derivesfrom_rels) == 2, "derives_from"

        for r in derivesfrom_rels:
            assert (r.object.uniquename == 'PAC:18136217') or (r.object.uniquename == 'PAC:18136225'), "derives_from"

        terms = {cvt.cvterm.name: cvt.cvterm.dbxref.db_id for cvt in derivesfrom[0].feature_cvterm_collection}

        assert len(terms) == 2, "gff>ontology_term loaded correctly"
        assert '000001' in terms, "gff>ontology_term loaded correctly"
        assert '00002' in terms, "gff>ontology_term loaded correctly"
        assert terms['000001'] == dbs['GO'], "gff>ontology_term loaded correctly"
        assert terms['00002'] == dbs['GO'], "gff>ontology_term loaded correctly"

        # Target location
        assert len(derivesfrom[0].featureloc_collection) == 2, "gff>target loc ok"
        if derivesfrom[0].featureloc_collection[0].fmin == 120:
            checkedloc = 0
        else:
            checkedloc = 1
        assert derivesfrom[0].featureloc_collection[checkedloc].fmin == 120, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].fmax == 320, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].strand == -1, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].rank == 1, "gff>gene located correctly"

    def test_load_gff_pepregex(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], re_protein="foo\\1-bar", re_protein_capture="PAC:([0-9]+)", no_seq_compute=True)

        # Check mrna
        rna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136219") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        # Check relationships
        parents = {x.object_id: x.type_id for x in rna_f.subject_in_relationships}
        assert len(parents) == 1, "mRNA relationships"
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')
        assert rna_f.subject_in_relationships[0].type_id == partofterm, "mRNA relationships"

        derivesfromterm = self.ci.get_cvterm_id('derives_from', 'sequence')
        peps = [x.subject_id for x in rna_f.object_in_relationships if x.type_id == derivesfromterm]
        assert len(peps) == 1, "mRNA relationships, single peptide"

        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=peps[0]) \
            .one()

        # Check pep feature
        assert pep_f.uniquename == "foo18136219-bar", "gff>pep loaded correctly"

    def test_load_gff_pepregex2(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], re_protein="foo\\1-bar", no_seq_compute=True)

        # Check mrna
        rna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136219") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        # Check relationships
        parents = {x.object_id: x.type_id for x in rna_f.subject_in_relationships}
        assert len(parents) == 1, "mRNA relationships"
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')
        assert rna_f.subject_in_relationships[0].type_id == partofterm, "mRNA relationships"

        derivesfromterm = self.ci.get_cvterm_id('derives_from', 'sequence')
        peps = [x.subject_id for x in rna_f.object_in_relationships if x.type_id == derivesfromterm]
        assert len(peps) == 1, "mRNA relationships, single peptide"

        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=peps[0]) \
            .one()

        # Check pep feature
        assert pep_f.uniquename == "fooPAC:18136219-bar", "gff>pep loaded correctly"

    def test_load_gff_twice(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        # Adding twice the same gff should not change anything in the db
        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015632m.g") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        # Check gene feature
        assert gene_f.dbxref_id is None, "gff>gene loaded correctly"
        assert gene_f.organism_id == org['organism_id'], "gff>gene loaded correctly"
        assert gene_f.name == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.uniquename == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.residues is None, "gff>gene loaded correctly"
        assert gene_f.seqlen is None, "gff>gene loaded correctly"
        assert gene_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>gene loaded correctly"
        assert gene_f.type_id == geneterm, "gff>gene loaded correctly"
        assert gene_f.is_analysis is False, "gff>gene loaded correctly"
        assert gene_f.is_obsolete is False, "gff>gene loaded correctly"

        # Check gene loc
        assert len(gene_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert gene_f.featureloc_collection[0].fmin == 4058459, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].fmax == 4062210, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmin_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmax_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].strand == 1, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].phase is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].residue_info is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].locgroup == 0, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].rank == 0, "gff>gene located correctly"

        src_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=gene_f.featureloc_collection[0].srcfeature_id) \
            .one()

        assert src_f.uniquename == "scaffold00001", "gff>gene loaded correctly"
        scaff1_id = src_f.feature_id

        # Check gene aliases
        exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
        syns = {synf.synonym.name: synf.synonym.type_id for synf in gene_f.feature_synonym_collection}

        assert len(syns) == 2, "gff>gene aliases loaded correctly"
        assert 'some-synonym' in syns, "gff>gene aliases loaded correctly"
        assert 'another synonym' in syns, "gff>gene aliases loaded correctly"
        assert syns['some-synonym'] == exactterm, "gff>gene aliases loaded correctly"
        assert syns['another synonym'] == exactterm, "gff>gene aliases loaded correctly"

        # Check gene dbxref
        dbs = self.ci.session.query(self.ci.model.db.db_id, self.ci.model.db.name, self.ci.model.db.description) \
            .filter((self.ci.model.db.name == 'GO') | (self.ci.model.db.name == 'FOOBAR') | (self.ci.model.db.name == 'FOOBARXX') | (self.ci.model.db.name == 'GFF_source'))
        for db in dbs:
            if db.name == "FOOBAR":
                assert db.description == "Added automatically by the GFF loader", "gff>gene dbxrefs db loaded correctly"

        dbs = {db.name: db.db_id for db in dbs}

        assert len(dbs) == 4, "gff>gene dbxrefs db loaded correctly"

        xrefs = {dbx.dbxref.accession: dbx.dbxref.db_id for dbx in gene_f.feature_dbxref_collection}

        assert len(xrefs) == 3, "gff>gene dbxrefs loaded correctly"
        assert '0061611' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert '6528B' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert 'phytozome6' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert xrefs['0061611'] == dbs['GO'], "gff>gene dbxrefs loaded correctly"
        assert xrefs['6528B'] == dbs['FOOBAR'], "gff>gene dbxrefs loaded correctly"
        assert xrefs['phytozome6'] == dbs['GFF_source'], "gff>gene dbxrefs loaded correctly"

        # Check gene featureprop
        expected = [
            'Gap___BLABLA___0',
            'Gap___BLOBLO___1',
            'Note___that\'s fantastic___0',
            'Note___really___1',
            'Poutrelle___test___1',
            'Poutrelle___lapinou___0',
        ]

        assert len(gene_f.featureprop_collection) == 6, "gff>gene loaded correctly"

        for prop in gene_f.featureprop_collection:
            assert prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank) in expected, "gff>gene loaded correctly"
            expected.remove(prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank))

        # Check mrna
        rna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136219") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        rnaterm = self.ci.get_cvterm_id('mRNA', 'sequence')

        # Check mRNA feature
        assert rna_f.dbxref_id is None, "gff>mRNA loaded correctly"
        assert rna_f.organism_id == org['organism_id'], "gff>mRNA loaded correctly"
        assert rna_f.name == "orange1.1g015615m", "gff>mRNA loaded correctly"
        assert rna_f.uniquename == "PAC:18136219", "gff>mRNA loaded correctly"
        assert rna_f.residues is None, "gff>mRNA loaded correctly"
        assert rna_f.seqlen is None, "gff>mRNA loaded correctly"
        assert rna_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>mRNA loaded correctly"
        assert rna_f.type_id == rnaterm, "gff>mRNA loaded correctly"
        assert rna_f.is_analysis is False, "gff>mRNA loaded correctly"
        assert rna_f.is_obsolete is False, "gff>mRNA loaded correctly"

        # Check mRNA loc
        assert len(rna_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert rna_f.featureloc_collection[0].fmin == 4058759, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].fmax == 4062210, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].is_fmin_partial is False, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].is_fmax_partial is False, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].strand == 1, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].phase is None, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].residue_info is None, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].locgroup == 0, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].rank == 0, "gff>mRNA located correctly"

        assert scaff1_id == rna_f.featureloc_collection[0].srcfeature_id, "gff>mRNA loaded correctly"

        # Check mRNA aliases
        exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
        syns = {synf.synonym.name: synf.synonym.type_id for synf in rna_f.feature_synonym_collection}

        assert len(syns) == 2, "gff>mRNA aliases loaded correctly"
        assert 'some-synonym' in syns, "gff>mRNA aliases loaded correctly"
        assert 'another synonym' in syns, "gff>mRNA aliases loaded correctly"
        assert syns['some-synonym'] == exactterm, "gff>mRNA aliases loaded correctly"
        assert syns['another synonym'] == exactterm, "gff>mRNA aliases loaded correctly"

        # Check mRNA dbxref
        xrefs = {dbx.dbxref.accession: dbx.dbxref.db_id for dbx in rna_f.feature_dbxref_collection}

        assert len(xrefs) == 3, "gff>mRNA dbxrefs loaded correctly"
        assert '0061621' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert '6528A' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert 'phytozome6' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['0061621'] == dbs['GO'], "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['6528A'] == dbs['FOOBARXX'], "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['phytozome6'] == dbs['GFF_source'], "gff>mRNA dbxrefs loaded correctly"

        # Check relationships
        parents = {x.object_id: x.type_id for x in rna_f.subject_in_relationships}
        assert len(parents) == 1, "mRNA relationships"
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')
        assert rna_f.subject_in_relationships[0].type_id == partofterm, "mRNA relationships"

        derivesfromterm = self.ci.get_cvterm_id('derives_from', 'sequence')
        peps = [x.subject_id for x in rna_f.object_in_relationships if x.type_id == derivesfromterm]
        assert len(peps) == 1, "mRNA relationships, single peptide"

        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=peps[0]) \
            .one()

        # Check pep feature
        pepterm = self.ci.get_cvterm_id('polypeptide', 'sequence')
        assert pep_f.dbxref_id is None, "gff>pep loaded correctly"
        assert pep_f.organism_id == org['organism_id'], "gff>pep loaded correctly"
        assert pep_f.name == "orange1.1g015615m", "gff>pep loaded correctly"
        assert pep_f.uniquename == "PAC:18136219-protein", "gff>pep loaded correctly"
        assert pep_f.residues is None, "gff>pep loaded correctly"
        assert pep_f.seqlen is None, "gff>pep loaded correctly"
        assert pep_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>pep loaded correctly"
        assert pep_f.type_id == pepterm, "gff>pep loaded correctly"
        assert pep_f.is_analysis is False, "gff>pep loaded correctly"
        assert pep_f.is_obsolete is False, "gff>pep loaded correctly"

        # Check pep loc
        assert len(pep_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].fmin == 4059234, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].fmax == 4061905, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].is_fmin_partial is False, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].is_fmax_partial is False, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].strand == 1, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].phase is None, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].residue_info is None, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].locgroup == 0, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].rank == 0, "gff>pep located correctly"

        assert scaff1_id == pep_f.featureloc_collection[0].srcfeature_id, "gff>pep loaded correctly"

        children = {x.subject_id: x for x in rna_f.object_in_relationships if x.type_id != derivesfromterm}
        assert len(children) == 30, "mRNA relationships, single peptide"  # relations to utr/cds/exons are duplicarted as they don't have an ID or Name attribute in gff => random uniquename in db

        cdsterm = self.ci.get_cvterm_id('CDS', 'sequence')
        exonterm = self.ci.get_cvterm_id('exon', 'sequence')
        utr3term = self.ci.get_cvterm_id('three_prime_UTR', 'sequence')
        utr5term = self.ci.get_cvterm_id('five_prime_UTR', 'sequence')
        for c in children:
            assert children[c].type_id == partofterm, "subsubfeatures"

            if children[c].subject.type_id == utr3term:
                subsub_f = children[c].subject

            assert children[c].subject.type_id in (cdsterm, exonterm, utr3term, utr5term), "subsubfeatures"

        # Check a subsubfeature
        assert subsub_f.dbxref_id is None, "gff>utr loaded correctly"
        assert subsub_f.organism_id == org['organism_id'], "gff>utr loaded correctly"
        assert subsub_f.name.endswith("-three_prime_UTR-scaffold00001:4061905..4062210"), "gff>utr loaded correctly"
        assert subsub_f.uniquename.endswith("-three_prime_UTR-scaffold00001:4061905..4062210"), "gff>utr loaded correctly"
        assert subsub_f.residues is None, "gff>utr loaded correctly"
        assert subsub_f.seqlen is None, "gff>utr loaded correctly"
        assert subsub_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>utr loaded correctly"
        assert subsub_f.type_id == utr3term, "gff>utr loaded correctly"
        assert subsub_f.is_analysis is False, "gff>utr loaded correctly"
        assert subsub_f.is_obsolete is False, "gff>utr loaded correctly"

        assert len(subsub_f.featureloc_collection) == 1, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].fmin == 4061905, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].fmax == 4062210, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].is_fmin_partial is False, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].is_fmax_partial is False, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].strand == 1, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].phase is None, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].residue_info is None, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].locgroup == 0, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].rank == 0, "gff>utr located correctly"

        # Check utr with 2 parents
        confused_child_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename='an_utr_with_two_parents') \
            .all()

        assert len(confused_child_f) == 1, "1 utr with 2 parents"

        confused_rels = confused_child_f[0].subject_in_relationships

        assert len(confused_rels) == 2, "1 utr with 2 parents"

        for r in confused_rels:
            assert (r.object.uniquename == 'PAC:18136239') or (r.object.uniquename == 'PAC:18136238'), "1 utr with 2 parents"

        # Check Derives_from
        derivesfrom = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename='some_special_cds') \
            .all()

        assert len(derivesfrom) == 1, "derives_from"

        derivesfrom_rels = derivesfrom[0].subject_in_relationships

        assert len(derivesfrom_rels) == 2, "derives_from"

        for r in derivesfrom_rels:
            assert (r.object.uniquename == 'PAC:18136217') or (r.object.uniquename == 'PAC:18136225'), "derives_from"

        terms = {cvt.cvterm.name: cvt.cvterm.dbxref.db_id for cvt in derivesfrom[0].feature_cvterm_collection}

        assert len(terms) == 2, "gff>ontology_term loaded correctly"
        assert '000001' in terms, "gff>ontology_term loaded correctly"
        assert '00002' in terms, "gff>ontology_term loaded correctly"
        assert terms['000001'] == dbs['GO'], "gff>ontology_term loaded correctly"
        assert terms['00002'] == dbs['GO'], "gff>ontology_term loaded correctly"

        # Target location
        assert len(derivesfrom[0].featureloc_collection) == 2, "gff>target loc ok"
        if derivesfrom[0].featureloc_collection[0].fmin == 120:
            checkedloc = 0
        else:
            checkedloc = 1
        assert derivesfrom[0].featureloc_collection[checkedloc].fmin == 120, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].fmax == 320, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].strand == -1, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].rank == 1, "gff>target loc ok"

    def test_load_gff_landmarktype(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        # there's a contig loaded by fasta and a supercontig in gff
        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015632m.g") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        # Check gene feature
        assert gene_f.dbxref_id is None, "gff>gene loaded correctly"
        assert gene_f.organism_id == org['organism_id'], "gff>gene loaded correctly"
        assert gene_f.name == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.uniquename == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.residues is None, "gff>gene loaded correctly"
        assert gene_f.seqlen is None, "gff>gene loaded correctly"
        assert gene_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>gene loaded correctly"
        assert gene_f.type_id == geneterm, "gff>gene loaded correctly"
        assert gene_f.is_analysis is False, "gff>gene loaded correctly"
        assert gene_f.is_obsolete is False, "gff>gene loaded correctly"

        # Check gene loc
        assert len(gene_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert gene_f.featureloc_collection[0].fmin == 4058459, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].fmax == 4062210, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmin_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmax_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].strand == 1, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].phase is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].residue_info is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].locgroup == 0, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].rank == 0, "gff>gene located correctly"

    def test_load_gff_nolandmark(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        # Should create the landmark
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], landmark_type="contig", no_seq_compute=True)

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015632m.g") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        # Check gene feature
        assert gene_f.dbxref_id is None, "gff>gene loaded correctly"
        assert gene_f.organism_id == org['organism_id'], "gff>gene loaded correctly"
        assert gene_f.name == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.uniquename == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.residues is None, "gff>gene loaded correctly"
        assert gene_f.seqlen is None, "gff>gene loaded correctly"
        assert gene_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>gene loaded correctly"
        assert gene_f.type_id == geneterm, "gff>gene loaded correctly"
        assert gene_f.is_analysis is False, "gff>gene loaded correctly"
        assert gene_f.is_obsolete is False, "gff>gene loaded correctly"

        # Check gene loc
        assert len(gene_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert gene_f.featureloc_collection[0].fmin == 4058459, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].fmax == 4062210, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmin_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmax_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].strand == 1, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].phase is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].residue_info is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].locgroup == 0, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].rank == 0, "gff>gene located correctly"

    @raises(Exception)
    def test_load_gff_nolandmark_fail(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        # Should create the landmark
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015632m.g") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        # Check gene feature
        assert gene_f.dbxref_id is None, "gff>gene loaded correctly"
        assert gene_f.organism_id == org['organism_id'], "gff>gene loaded correctly"
        assert gene_f.name == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.uniquename == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.residues is None, "gff>gene loaded correctly"
        assert gene_f.seqlen is None, "gff>gene loaded correctly"
        assert gene_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>gene loaded correctly"
        assert gene_f.type_id == geneterm, "gff>gene loaded correctly"
        assert gene_f.is_analysis is False, "gff>gene loaded correctly"
        assert gene_f.is_obsolete is False, "gff>gene loaded correctly"

        # Check gene loc
        assert len(gene_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert gene_f.featureloc_collection[0].fmin == 4058459, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].fmax == 4062210, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmin_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmax_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].strand == 1, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].phase is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].residue_info is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].locgroup == 0, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].rank == 0, "gff>gene located correctly"

    def test_load_gff_match(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')
        an_match = self._create_fake_an('matches')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)
        self.ci.feature.load_gff(gff="./test-data/matches.gff", analysis_id=an_match['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        # Check match
        match_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238-protein_XP_012228303.1_match_0001") \
            .one()

        matchterm = self.ci.get_cvterm_id('match', 'sequence')
        matchpartterm = self.ci.get_cvterm_id('match_part', 'sequence')

        assert match_f.dbxref_id is None, "gff>match loaded correctly"
        assert match_f.organism_id == org['organism_id'], "gff>match loaded correctly"
        assert match_f.name == "PAC:18136238-protein_XP_012228303.1_match_0001", "gff>match loaded correctly"
        assert match_f.uniquename == "PAC:18136238-protein_XP_012228303.1_match_0001", "gff>match loaded correctly"
        assert match_f.residues is None, "gff>match loaded correctly"
        assert match_f.seqlen is None, "gff>match loaded correctly"
        assert match_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>match loaded correctly"
        assert match_f.type_id == matchterm, "gff>match loaded correctly"
        assert match_f.is_analysis is False, "gff>match loaded correctly"
        assert match_f.is_obsolete is False, "gff>match loaded correctly"

        assert len(match_f.featureloc_collection) == 1, "gff>match loaded correctly"
        assert match_f.featureloc_collection[0].srcfeature.uniquename == 'PAC:18136238-protein', "gff>match loaded correctly"
        assert match_f.featureloc_collection[0].fmin == 50, "gff>match loaded correctly"
        assert match_f.featureloc_collection[0].fmax == 325, "gff>match loaded correctly"
        assert match_f.featureloc_collection[0].strand is None, "gff>match loaded correctly"
        assert match_f.featureloc_collection[0].phase == 0, "gff>match loaded correctly"

        # Check relationships
        assert len(match_f.object_in_relationships) == 1, "match_part relationship"
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')
        assert match_f.object_in_relationships[0].type_id == partofterm, "match_part relationship"

        # Check gene featureprop
        expected = [
            'e-value___1e-27___0',
            'hit_description___PREDICTED: uncharacterized protein LOC105675603 [Linepithema humile]___0',
            'hit_name___XP_012228303.1___0'
        ]

        assert len(match_f.featureprop_collection) == 3, "gff>match loaded correctly"

        for prop in match_f.featureprop_collection:
            assert prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank) in expected, "gff>match loaded correctly"
            expected.remove(prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank))

        match_part = match_f.object_in_relationships[0].subject
        assert match_part.dbxref_id is None, "gff>match loaded correctly"
        assert match_part.organism_id == org['organism_id'], "gff>match loaded correctly"
        assert match_part.name == "PAC:18136238-protein_XP_012228303.1_match_0001_1", "gff>match loaded correctly"
        assert match_part.uniquename == "PAC:18136238-protein_XP_012228303.1_match_0001_1", "gff>match loaded correctly"
        assert match_part.residues is None, "gff>match loaded correctly"
        assert match_part.seqlen is None, "gff>match loaded correctly"
        assert match_part.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>match loaded correctly"
        assert match_part.type_id == matchpartterm, "match_part"
        assert match_part.is_analysis is False, "gff>match loaded correctly"
        assert match_part.is_obsolete is False, "gff>match loaded correctly"

        assert len(match_part.featureloc_collection) == 1, "gff>match loaded correctly"

        assert match_part.featureloc_collection[0].srcfeature.uniquename == 'PAC:18136238-protein', "gff>match loaded correctly"
        assert match_part.featureloc_collection[0].fmin == 50, "gff>match loaded correctly"
        assert match_part.featureloc_collection[0].fmax == 325, "gff>match loaded correctly"
        assert match_part.featureloc_collection[0].strand is None, "gff>match loaded correctly"
        assert match_part.featureloc_collection[0].phase == 1, "gff>match loaded correctly"

        # Check gene featureprop
        expected = [
            'target___XP_012228303.1+21+302___0'
        ]

        assert len(match_part.featureprop_collection) == 1, "gff>match loaded correctly"

        for prop in match_part.featureprop_collection:
            assert prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank) in expected, "gff>match loaded correctly"
            expected.remove(prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank))

        # Check analysisfeature
        assert match_f.analysisfeature_collection[0].analysis_id == an_match['analysis_id'], "gff>match loaded correctly"
        assert match_f.analysisfeature_collection[0].significance == 303, "gff>match loaded correctly"
        assert match_part.analysisfeature_collection[0].analysis_id == an_match['analysis_id'], "gff>match loaded correctly"
        assert match_part.analysisfeature_collection[0].significance == 303, "gff>match loaded correctly"

    def test_load_gff_withpepfasta(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'])
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], fasta="./test-data/prots.fa")

        # Check pep
        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238-protein") \
            .one()

        assert pep_f.residues == "SGTRGVDFSVFDC", "gff>fasta seq loaded correctly"
        assert pep_f.seqlen == 13, "gff>fasta seq loaded correctly"
        assert pep_f.md5checksum == "744bbb7c3f619a479ea90b4e9f627bd1", "gff>fasta seq loaded correctly"

    def test_load_gff_withlandmarkfasta(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], fasta="./test-data/genome.fa", landmark_type="supercontig")

        # Check landmark
        lm_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="scaffold00001") \
            .one()

        assert lm_f.residues.startswith("TTTTGTATTCTATGTCCTCTGATCTTT"), "gff>fasta seq loaded correctly"
        assert lm_f.seqlen == 5927163, "gff>fasta seq loaded correctly"
        assert lm_f.md5checksum == "80db0e5ccdc07e200c035d23c5951271", "gff>fasta seq loaded correctly"

        # Check mrna
        mrna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238") \
            .one()

        assert mrna_f.residues.startswith("AAAGGAATTGAGTTTCATTAAGAATTTAAATAAAACAATGTCATAATCCGGGTATTTGGAATATT"), "gff>fasta seq loaded correctly"
        assert mrna_f.seqlen == 1212, "gff>fasta seq loaded correctly"
        assert mrna_f.md5checksum == "ad0d8a5031b63bacfe23296c80072550", "gff>fasta seq loaded correctly"

        # Check pep
        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238-protein") \
            .one()

        assert pep_f.residues.startswith("KGIEFH*EFK*NNVIIRVFGIFKLQPGLVVMQRPTR*QNDNLALVLGFRSFVHSFSS*AKANWNLTKCNAYTSSEPEQHSSYKXXXXXXXXXXXXXXXXXXXXXX"), "gff>fasta seq loaded correctly"
        assert pep_f.seqlen == 404, "gff>fasta seq loaded correctly"
        assert pep_f.md5checksum == "fbf522da5203405c620eb708afd3cc9f", "gff>fasta seq loaded correctly"

    def test_load_gff_withlandmarkonly(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        # Here the gff will create a supercontig, and other features will be mapped on it.
        # No fasta => no computed seq
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], landmark_type="supercontig")

        # Check landmark
        lm_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="scaffold00001") \
            .one()

        assert lm_f.residues is None, "gff>fasta seq loaded correctly"
        assert lm_f.seqlen is None, "gff>fasta seq loaded correctly"
        assert lm_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>fasta seq loaded correctly"

        # Check mrna
        mrna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238") \
            .one()

        assert mrna_f.residues is None, "gff>fasta seq loaded correctly"
        assert mrna_f.seqlen is None, "gff>fasta seq loaded correctly"
        assert mrna_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>fasta seq loaded correctly"

        # Check pep
        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238-protein") \
            .one()

        assert pep_f.residues is None, "gff>fasta seq loaded correctly"
        assert pep_f.seqlen is None, "gff>fasta seq loaded correctly"
        assert pep_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>fasta seq loaded correctly"

    @raises(Exception)
    def test_load_gff_withoutlandmark(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        # Here the gff will create a supercontig, and the loader features will try to map on a contig => fail
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'])

    def test_load_gff_withwronglandmarkonly(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        # Here the gff will create a supercontig and a contig, and other features will be mapped on contig.
        # No fasta => no computed seq
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], landmark_type="contig")

        contigterm = self.ci.get_cvterm_id('contig', 'sequence')

        # Check landmark
        lm_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="scaffold00001") \
            .filter_by(type_id=contigterm) \
            .one()

        assert lm_f.residues is None, "gff>fasta seq loaded correctly"
        assert lm_f.seqlen is None, "gff>fasta seq loaded correctly"
        assert lm_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>fasta seq loaded correctly"

        # Check mrna
        mrna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238") \
            .one()

        assert mrna_f.residues is None, "gff>fasta seq loaded correctly"
        assert mrna_f.seqlen is None, "gff>fasta seq loaded correctly"
        assert mrna_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>fasta seq loaded correctly"

        # Check pep
        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238-protein") \
            .one()

        assert pep_f.residues is None, "gff>fasta seq loaded correctly"
        assert pep_f.seqlen is None, "gff>fasta seq loaded correctly"
        assert pep_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>fasta seq loaded correctly"

    def test_load_gff_relranks(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], landmark_type="contig")

        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')

        # Check mrna
        mrna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136238") \
            .one()

        rels = mrna_f.object_in_relationships

        locsorted_rels = []
        for rel in rels:
            if rel.type_id == partofterm:
                locsorted_rels.append((rel.subject.featureloc_collection[0].fmin, rel.rank))

        locsorted_rels = sorted(locsorted_rels, key=lambda x: x[0])

        sorted_rels = sorted(locsorted_rels, key=lambda x: x[1])

        assert locsorted_rels == sorted_rels, "children sorted correctly"

    def test_load_gff_addonly(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], fasta="./test-data/genome.fa", landmark_type="supercontig", no_seq_compute=True, add_only=True)

        gene_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="orange1.1g015632m.g") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        # Check gene feature
        assert gene_f.dbxref_id is None, "gff>gene loaded correctly"
        assert gene_f.organism_id == org['organism_id'], "gff>gene loaded correctly"
        assert gene_f.name == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.uniquename == "orange1.1g015632m.g", "gff>gene loaded correctly"
        assert gene_f.residues is None, "gff>gene loaded correctly"
        assert gene_f.seqlen is None, "gff>gene loaded correctly"
        assert gene_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>gene loaded correctly"
        assert gene_f.type_id == geneterm, "gff>gene loaded correctly"
        assert gene_f.is_analysis is False, "gff>gene loaded correctly"
        assert gene_f.is_obsolete is False, "gff>gene loaded correctly"

        # Check gene loc
        assert len(gene_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert gene_f.featureloc_collection[0].fmin == 4058459, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].fmax == 4062210, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmin_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].is_fmax_partial is False, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].strand == 1, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].phase is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].residue_info is None, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].locgroup == 0, "gff>gene located correctly"
        assert gene_f.featureloc_collection[0].rank == 0, "gff>gene located correctly"

        src_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=gene_f.featureloc_collection[0].srcfeature_id) \
            .one()

        assert src_f.uniquename == "scaffold00001", "gff>gene loaded correctly"
        scaff1_id = src_f.feature_id

        # Check gene aliases
        exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
        syns = {synf.synonym.name: synf.synonym.type_id for synf in gene_f.feature_synonym_collection}

        assert len(syns) == 2, "gff>gene aliases loaded correctly"
        assert 'some-synonym' in syns, "gff>gene aliases loaded correctly"
        assert 'another synonym' in syns, "gff>gene aliases loaded correctly"
        assert syns['some-synonym'] == exactterm, "gff>gene aliases loaded correctly"
        assert syns['another synonym'] == exactterm, "gff>gene aliases loaded correctly"

        # Check gene dbxref
        dbs = self.ci.session.query(self.ci.model.db.db_id, self.ci.model.db.name, self.ci.model.db.description) \
            .filter((self.ci.model.db.name == 'GO') | (self.ci.model.db.name == 'FOOBAR') | (self.ci.model.db.name == 'FOOBARXX') | (self.ci.model.db.name == 'GFF_source'))
        for db in dbs:
            if db.name == "FOOBAR":
                assert db.description == "Added automatically by the GFF loader", "gff>gene dbxrefs db loaded correctly"

        dbs = {db.name: db.db_id for db in dbs}

        assert len(dbs) == 4, "gff>gene dbxrefs db loaded correctly"

        xrefs = {dbx.dbxref.accession: dbx.dbxref.db_id for dbx in gene_f.feature_dbxref_collection}

        assert len(xrefs) == 3, "gff>gene dbxrefs loaded correctly"
        assert '0061611' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert '6528B' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert 'phytozome6' in xrefs, "gff>gene dbxrefs loaded correctly"
        assert xrefs['0061611'] == dbs['GO'], "gff>gene dbxrefs loaded correctly"
        assert xrefs['6528B'] == dbs['FOOBAR'], "gff>gene dbxrefs loaded correctly"
        assert xrefs['phytozome6'] == dbs['GFF_source'], "gff>gene dbxrefs loaded correctly"

        # Check gene featureprop
        expected = [
            'Gap___BLABLA___0',
            'Gap___BLOBLO___1',
            'Note___that\'s fantastic___0',
            'Note___really___1',
            'Poutrelle___test___1',
            'Poutrelle___lapinou___0',
        ]

        assert len(gene_f.featureprop_collection) == 6, "gff>gene loaded correctly"

        for prop in gene_f.featureprop_collection:
            assert prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank) in expected, "gff>gene loaded correctly"
            expected.remove(prop.cvterm.name + '___' + prop.value + '___' + str(prop.rank))

        # Check mrna
        rna_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename="PAC:18136219") \
            .join(self.ci.model.featureloc, self.ci.model.featureloc.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.feature_synonym, self.ci.model.feature_synonym.feature_id == self.ci.model.feature.feature_id) \
            .join(self.ci.model.synonym, self.ci.model.feature_synonym.synonym_id == self.ci.model.synonym.synonym_id) \
            .one()

        rnaterm = self.ci.get_cvterm_id('mRNA', 'sequence')

        # Check mRNA feature
        assert rna_f.dbxref_id is None, "gff>mRNA loaded correctly"
        assert rna_f.organism_id == org['organism_id'], "gff>mRNA loaded correctly"
        assert rna_f.name == "orange1.1g015615m", "gff>mRNA loaded correctly"
        assert rna_f.uniquename == "PAC:18136219", "gff>mRNA loaded correctly"
        assert rna_f.residues is None, "gff>mRNA loaded correctly"
        assert rna_f.seqlen is None, "gff>mRNA loaded correctly"
        assert rna_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>mRNA loaded correctly"
        assert rna_f.type_id == rnaterm, "gff>mRNA loaded correctly"
        assert rna_f.is_analysis is False, "gff>mRNA loaded correctly"
        assert rna_f.is_obsolete is False, "gff>mRNA loaded correctly"

        # Check mRNA loc
        assert len(rna_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert rna_f.featureloc_collection[0].fmin == 4058759, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].fmax == 4062210, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].is_fmin_partial is False, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].is_fmax_partial is False, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].strand == 1, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].phase is None, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].residue_info is None, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].locgroup == 0, "gff>mRNA located correctly"
        assert rna_f.featureloc_collection[0].rank == 0, "gff>mRNA located correctly"

        assert scaff1_id == rna_f.featureloc_collection[0].srcfeature_id, "gff>mRNA loaded correctly"

        # Check mRNA aliases
        exactterm = self.ci.get_cvterm_id('exact', 'synonym_type')
        syns = {synf.synonym.name: synf.synonym.type_id for synf in rna_f.feature_synonym_collection}

        assert len(syns) == 2, "gff>mRNA aliases loaded correctly"
        assert 'some-synonym' in syns, "gff>mRNA aliases loaded correctly"
        assert 'another synonym' in syns, "gff>mRNA aliases loaded correctly"
        assert syns['some-synonym'] == exactterm, "gff>mRNA aliases loaded correctly"
        assert syns['another synonym'] == exactterm, "gff>mRNA aliases loaded correctly"

        # Check mRNA dbxref
        xrefs = {dbx.dbxref.accession: dbx.dbxref.db_id for dbx in rna_f.feature_dbxref_collection}

        assert len(xrefs) == 3, "gff>mRNA dbxrefs loaded correctly"
        assert '0061621' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert '6528A' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert 'phytozome6' in xrefs, "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['0061621'] == dbs['GO'], "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['6528A'] == dbs['FOOBARXX'], "gff>mRNA dbxrefs loaded correctly"
        assert xrefs['phytozome6'] == dbs['GFF_source'], "gff>mRNA dbxrefs loaded correctly"

        # Check relationships
        parents = {x.object_id: x.type_id for x in rna_f.subject_in_relationships}
        assert len(parents) == 1, "mRNA relationships"
        partofterm = self.ci.get_cvterm_id('part_of', 'sequence')
        assert rna_f.subject_in_relationships[0].type_id == partofterm, "mRNA relationships"

        derivesfromterm = self.ci.get_cvterm_id('derives_from', 'sequence')
        peps = [x.subject_id for x in rna_f.object_in_relationships if x.type_id == derivesfromterm]
        assert len(peps) == 1, "mRNA relationships, single peptide"

        pep_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(feature_id=peps[0]) \
            .one()

        # Check pep feature
        pepterm = self.ci.get_cvterm_id('polypeptide', 'sequence')
        assert pep_f.dbxref_id is None, "gff>pep loaded correctly"
        assert pep_f.organism_id == org['organism_id'], "gff>pep loaded correctly"
        assert pep_f.name == "orange1.1g015615m", "gff>pep loaded correctly"
        assert pep_f.uniquename == "PAC:18136219-protein", "gff>pep loaded correctly"
        assert pep_f.residues is None, "gff>pep loaded correctly"
        assert pep_f.seqlen is None, "gff>pep loaded correctly"
        assert pep_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>pep loaded correctly"
        assert pep_f.type_id == pepterm, "gff>pep loaded correctly"
        assert pep_f.is_analysis is False, "gff>pep loaded correctly"
        assert pep_f.is_obsolete is False, "gff>pep loaded correctly"

        # Check pep loc
        assert len(pep_f.featureloc_collection) == 1, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].fmin == 4059234, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].fmax == 4061905, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].is_fmin_partial is False, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].is_fmax_partial is False, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].strand == 1, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].phase is None, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].residue_info is None, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].locgroup == 0, "gff>pep located correctly"
        assert pep_f.featureloc_collection[0].rank == 0, "gff>pep located correctly"

        assert scaff1_id == pep_f.featureloc_collection[0].srcfeature_id, "gff>pep loaded correctly"

        children = {x.subject_id: x for x in rna_f.object_in_relationships if x.type_id != derivesfromterm}
        assert len(children) == 15, "mRNA relationships, single peptide"

        cdsterm = self.ci.get_cvterm_id('CDS', 'sequence')
        exonterm = self.ci.get_cvterm_id('exon', 'sequence')
        utr3term = self.ci.get_cvterm_id('three_prime_UTR', 'sequence')
        utr5term = self.ci.get_cvterm_id('five_prime_UTR', 'sequence')
        for c in children:
            assert children[c].type_id == partofterm, "subsubfeatures"

            if children[c].subject.type_id == utr3term:
                subsub_f = children[c].subject

            assert children[c].subject.type_id in (cdsterm, exonterm, utr3term, utr5term), "subsubfeatures"

        # Check a subsubfeature
        assert subsub_f.dbxref_id is None, "gff>utr loaded correctly"
        assert subsub_f.organism_id == org['organism_id'], "gff>utr loaded correctly"
        assert subsub_f.name.endswith("-three_prime_UTR-scaffold00001:4061905..4062210"), "gff>utr loaded correctly"
        assert subsub_f.uniquename.endswith("-three_prime_UTR-scaffold00001:4061905..4062210"), "gff>utr loaded correctly"
        assert subsub_f.residues is None, "gff>utr loaded correctly"
        assert subsub_f.seqlen is None, "gff>utr loaded correctly"
        assert subsub_f.md5checksum == "d41d8cd98f00b204e9800998ecf8427e", "gff>utr loaded correctly"
        assert subsub_f.type_id == utr3term, "gff>utr loaded correctly"
        assert subsub_f.is_analysis is False, "gff>utr loaded correctly"
        assert subsub_f.is_obsolete is False, "gff>utr loaded correctly"

        assert len(subsub_f.featureloc_collection) == 1, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].fmin == 4061905, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].fmax == 4062210, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].is_fmin_partial is False, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].is_fmax_partial is False, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].strand == 1, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].phase is None, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].residue_info is None, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].locgroup == 0, "gff>utr located correctly"
        assert subsub_f.featureloc_collection[0].rank == 0, "gff>utr located correctly"

        # Check utr with 2 parents
        confused_child_f = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename='an_utr_with_two_parents') \
            .all()

        assert len(confused_child_f) == 1, "1 utr with 2 parents"

        confused_rels = confused_child_f[0].subject_in_relationships

        assert len(confused_rels) == 2, "1 utr with 2 parents"

        for r in confused_rels:
            assert (r.object.uniquename == 'PAC:18136239') or (r.object.uniquename == 'PAC:18136238'), "1 utr with 2 parents"

        # Check Derives_from
        derivesfrom = self.ci.session.query(self.ci.model.feature) \
            .filter_by(uniquename='some_special_cds') \
            .all()

        assert len(derivesfrom) == 1, "derives_from"

        derivesfrom_rels = derivesfrom[0].subject_in_relationships

        assert len(derivesfrom_rels) == 2, "derives_from"

        for r in derivesfrom_rels:
            assert (r.object.uniquename == 'PAC:18136217') or (r.object.uniquename == 'PAC:18136225'), "derives_from"

        terms = {cvt.cvterm.name: cvt.cvterm.dbxref.db_id for cvt in derivesfrom[0].feature_cvterm_collection}

        assert len(terms) == 2, "gff>ontology_term loaded correctly"
        assert '000001' in terms, "gff>ontology_term loaded correctly"
        assert '00002' in terms, "gff>ontology_term loaded correctly"
        assert terms['000001'] == dbs['GO'], "gff>ontology_term loaded correctly"
        assert terms['00002'] == dbs['GO'], "gff>ontology_term loaded correctly"

        # Target location
        assert len(derivesfrom[0].featureloc_collection) == 2, "gff>target loc ok"
        if derivesfrom[0].featureloc_collection[0].fmin == 120:
            checkedloc = 0
        else:
            checkedloc = 1
        assert derivesfrom[0].featureloc_collection[checkedloc].fmin == 120, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].fmax == 320, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].strand == -1, "gff>target loc ok"
        assert derivesfrom[0].featureloc_collection[checkedloc].rank == 1, "gff>gene located correctly"

    def test_load_gff_twice_addonly(self):
        org = self._create_fake_org()
        an_gff = self._create_fake_an('gff')

        # Adding twice the same gff with --add_only should raise some exception
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], fasta="./test-data/genome.fa", landmark_type="supercontig", no_seq_compute=True)

        try:
            self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True, add_only=True)
        except Exception:
            self.ci.session.rollback()
            assert True

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
