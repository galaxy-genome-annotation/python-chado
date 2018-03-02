import unittest

from nose.tools import raises

from . import ci


class CvtermTest(unittest.TestCase):

    def test_get_id_by_name(self):

        term = self.ci.get_cvterm_id(name='contig', cv='sequence')

        assert term > 0, "Got contig cvterm id"

    @raises(Exception)
    def test_get_id_by_synonym_fail(self):

        print(self.ci.get_cvterm_id(name='lives inside of', cv='relationship'))

    def test_get_id_by_synonym(self):

        term = self.ci.get_cvterm_id(name='lives inside of', cv='relationship', allow_synonyms=True)

        assert term == 415, "Got cvterm id by synonym"

    def test_get_id_by_name_synmultiple(self):

        term = self.ci.get_cvterm_id(name='three_prime_UTR', cv='sequence', allow_synonyms=True)

        assert term > 0, "Got contig cvterm id multiple"

    def test_get_name(self):

        term = self.ci.get_cvterm_name(415)

        assert term == 'endoparasite_of', "Got cvterm name"

    @raises(Exception)
    def test_get_name_fail(self):

        self.ci.get_cvterm_name(4150000000)

    def create_term(self):
        self.ci.create_cvterm("my_term", 'my_cv', 'my_db', term_definition="my_term is cool", cv_definition="my_cv is awesome", db_definition="look at my_db")

        term = self.ci.session.query(self.ci.model.cvterm) \
            .filter_by(name='my_term') \
            .one()

        assert term.definition == "my_term is cool", 'cvterm creation'
        assert term.cv.name == "my_cv", 'cvterm creation'
        assert term.cv.definition == "my_cv is awesome", 'cvterm creation'
        assert term.dbxref.accession == "my_term", 'cvterm creation'
        assert term.dbxref.db.name == "my_db", 'cvterm creation'
        assert term.dbxref.db.definition == "look at my_db", 'cvterm creation'

    def setUp(self):
        self.ci = ci
