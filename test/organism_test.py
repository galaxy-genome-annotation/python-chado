import unittest

from . import ci


class OrganismTest(unittest.TestCase):

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

    def setUp(self):
        self.ci = ci
        self.ci.organism.delete_organisms()

        self.ci.session.commit()

    def tearDown(self):
        self.ci.organism.delete_organisms()

        self.ci.session.commit()
