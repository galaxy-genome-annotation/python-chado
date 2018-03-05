import sys
import unittest

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from . import ci


class ExportTest(unittest.TestCase):

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

    def test_export_fasta(self):
        org = self._create_fake_org()
        an = self._create_fake_an()

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')

        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out

            self.ci.export.export_fasta(org['organism_id'])
            output = out.getvalue().strip()
            assert output.startswith(">scaffold00001"), "fasta exported"
            assert "TTTTGTATTCTATGTCCTCTGATCTTTATACTTCTTCATTTT" in output, "fasta exported"
        finally:
            sys.stdout = saved_stdout

    def test_export_gff(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out

            self.ci.export.export_gff3(org['organism_id'])
            output = out.getvalue().strip()
            assert output.startswith("##gff-version 3\n##sequence-region Testorg 1 1"), "gff exported"
            assert "ID=PAC:18136238;PACid=18136238;Parent=orange1.1g015633m.g" in output, "gff exported"
        finally:
            sys.stdout = saved_stdout

    def test_export_gbk(self):
        org = self._create_fake_org()
        an = self._create_fake_an()
        an_gff = self._create_fake_an('gff')

        self.ci.feature.load_fasta(fasta="./test-data/genome.fa", analysis_id=an['analysis_id'], organism_id=org['organism_id'], sequence_type='supercontig')
        self.ci.feature.load_gff(gff="./test-data/annot.gff", analysis_id=an_gff['analysis_id'], organism_id=org['organism_id'], no_seq_compute=True)

        saved_stdout = sys.stdout
        try:
            out = StringIO()
            sys.stdout = out

            self.ci.export.export_gbk(org['organism_id'])
            output = out.getvalue().strip()
            assert output.startswith("LOCUS       Testorg              5927163 bp"), "gbk exported"
            assert "three_prime_UTR complement(121..320)" in output, "gbk exported"
        finally:
            sys.stdout = saved_stdout

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
