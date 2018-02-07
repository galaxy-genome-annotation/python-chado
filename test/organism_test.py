import unittest
from chado import *

class OrganismTest(unittest.TestCase):

    @staticmethod
    def test_get_orgs():

        ci = ChadoInstance(dbuser="postgres", dbpass="postgres", dbname="postgres")

        orgs = [x['common_name'] for x in ci.organism.get_organisms()]

        assert 'human' in orgs, "human organism is loaded"
        assert 'yeast' in orgs, "yeast organism is loaded"


    @staticmethod
    def test_add_orgs():

        ci = ChadoInstance(dbuser="postgres", dbpass="postgres", dbname="postgres")

        genus = "Testus"
        common = "Testorg"
        abbr = "Ttesta"
        species = "testa"
        comment = "A test org"

        ci.organism.add_organism(genus=genus, common=common, abbr=abbr, species=species, comment=comment)

        org = ci.organism.get_organisms(abbr=abbr)[0]

        assert org["genus"] == genus, "org properly created"
        assert org["common_name"] == common, "org properly created"
        assert org["abbreviation"] == abbr, "org properly created"
        assert org["species"] == species, "org properly created"
        assert org["comment"] == comment, "org properly created"
