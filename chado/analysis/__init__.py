"""
Contains possible interactions with the Chado Analysis Module
"""
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from chado.client import Client
from chado.models import *
from datetime import datetime


class AnalysisClient(Client):
    """
    Access to the chado analysis table
    """

    def add_analysis(self, name, program, programversion, algorithm, sourcename, sourceversion, sourceuri, date_executed=None):
        """
        Create an analysis

        :type name: str
        :param name: analysis name

        :type program: str
        :param program: analysis program

        :type programversion: str
        :param programversion: analysis programversion

        :type algorithn: str
        :param algorithn: analysis algorithn

        :type sourcename: str
        :param sourcename: analysis sourcename

        :type sourcevesion: str
        :param sourcevesion: analysis sourcevesion

        :type sourceuri: str
        :param sourceuri: analysis sourceuri

        :type date_executed: str
        :param date_executed: analysis date_executed (yyyy-mm-dd)

        :rtype: dict
        :return: Analysis information
        """
        # check if the analysis exists
        res = self.session.query(Analysis).filter_by(name=name, program=program, programversion=programversion, sourcename=sourcename)

        if res.count() > 0:
            raise Exception("Found a preexisting analysis with the same attributes in the database")

        date = datetime.today()
        if date_executed:
            date = datetime.strptime(date_executed, '%Y-%m-%d')

        newa = Analysis()
        newa.name = name
        newa.program = program
        newa.programversion = programversion
        newa.algorithm = algorithm
        newa.sourcename = sourcename
        newa.sourceversion = sourceversion
        newa.sourceuri = sourceuri
        newa.timeexecuted = date
        self.session.add(newa)
        self.session.commit()

        return {
            'analysis_id': newa.analysis_id,
            'name': newa.name,
            'program': newa.program,
            'programversion': newa.programversion,
            'algorithm': newa.algorithm,
            'sourcename': newa.sourcename,
            'sourceversion': newa.sourceversion,
            'sourceuri': newa.sourceuri,
            'timeexecuted': newa.timeexecuted.isoformat(),
        }
