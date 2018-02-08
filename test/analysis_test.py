import unittest
from chado import *

class AnalysisTest(unittest.TestCase):

    @staticmethod
    def test_add_analysis():

        ci = ChadoInstance(dbuser="postgres", dbpass="postgres", dbname="postgres")

        name = "analysis x"
        program = "Magic"
        programversion = "1.0"
        algorithm = "mind"
        sourcename = "src"
        sourceversion = "2.1beta"
        sourceuri = "http://example.org/"
        date_executed = "2018-02-03"

        ci.analysis.add_analysis(name=name, program=program, programversion=programversion, algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri, date_executed=date_executed)

        ana = ci.analysis.get_analyses(name=name)[0]

        assert ana["name"] == name, "analysis properly created"
        assert ana["program"] == program, "analysis properly created"
        assert ana["programversion"] == programversion, "analysis properly created"
        assert ana["algorithm"] == algorithm, "analysis properly created"
        assert ana["sourcename"] == sourcename, "analysis properly created"
        assert ana["sourceversion"] == sourceversion, "analysis properly created"
        assert ana["sourceuri"] == sourceuri, "analysis properly created"
        assert ana["timeexecuted"] == '2018-02-03 00:00:00', "analysis properly created"
