import unittest

from . import ci


class AnalysisTest(unittest.TestCase):

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

        ana = self.ci.analysis.add_analysis(name=name, program=program, programversion=programversion, algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri, description=description, date_executed=date_executed)

        assert ana["analysis_id"] > 0, "ana properly created"

        ana = self.ci.analysis.get_analyses(name=name)[0]

        assert ana["name"] == name, "analysis properly created"
        assert ana["program"] == program, "analysis properly created"
        assert ana["programversion"] == programversion, "analysis properly created"
        assert ana["algorithm"] == algorithm, "analysis properly created"
        assert ana["sourcename"] == sourcename, "analysis properly created"
        assert ana["sourceversion"] == sourceversion, "analysis properly created"
        assert ana["sourceuri"] == sourceuri, "analysis properly created"
        assert ana["description"] == description, "analysis properly created"
        assert ana["timeexecuted"] == '2018-02-03 00:00:00', "analysis properly created"

        self.ci.analysis.delete_analyses(name=name)

        ans = self.ci.analysis.get_analyses(name=name)

        assert len(ans) == 0, "analysis properly deleted"

    def setUp(self):
        self.ci = ci
        self.ci.analysis.delete_analyses()

        self.ci.session.commit()

    def tearDown(self):
        self.ci.analysis.delete_analyses()

        self.ci.session.commit()
