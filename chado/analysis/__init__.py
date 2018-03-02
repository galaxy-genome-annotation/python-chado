"""
Contains possible interactions with the Chado Analysis Module
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

from datetime import datetime

from chado.client import Client

from future import standard_library

standard_library.install_aliases()


class AnalysisClient(Client):
    """
    Access to the chado analysis table
    """

    def add_analysis(self, name, program, programversion, sourcename, algorithm=None, sourceversion=None, sourceuri=None, description=None, date_executed=None):
        """
        Create an analysis

        :type name: str
        :param name: analysis name

        :type program: str
        :param program: analysis program

        :type programversion: str
        :param programversion: analysis programversion

        :type algorithm: str
        :param algorithm: analysis algorithm

        :type sourcename: str
        :param sourcename: analysis sourcename

        :type sourceversion: str
        :param sourceversion: analysis sourceversion

        :type sourceuri: str
        :param sourceuri: analysis sourceuri

        :type description: str
        :param description: analysis description

        :type date_executed: str
        :param date_executed: analysis date_executed (yyyy-mm-dd)

        :rtype: dict
        :return: Analysis information
        """
        # check if the analysis exists
        res = self.session.query(self.model.analysis).filter_by(name=name, program=program, programversion=programversion, sourcename=sourcename)

        if res.count() > 0:
            raise Exception("Found a preexisting analysis with the same attributes in the database")

        date = datetime.today()
        if date_executed:
            date = datetime.strptime(date_executed, '%Y-%m-%d')

        newa = self.model.analysis()
        newa.name = name
        newa.program = program
        newa.programversion = programversion
        newa.algorithm = algorithm
        newa.sourcename = sourcename
        newa.sourceversion = sourceversion
        newa.sourceuri = sourceuri
        newa.description = description
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
            'description': newa.description,
            'timeexecuted': newa.timeexecuted.isoformat(),
        }

    def get_analyses(self, analysis_id=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None):
        """
        Get all or some analyses

        :type analysis_id: int
        :param analysis_id: analysis_id filter

        :type name: str
        :param name: analysis name filter

        :type program: str
        :param program: analysis program filter

        :type programversion: str
        :param programversion: analysis programversion filter

        :type algorithm: str
        :param algorithm: analysis algorithm filter

        :type sourcename: str
        :param sourcename: analysis sourcename filter

        :type sourceversion: str
        :param sourceversion: analysis sourceversion filter

        :type sourceuri: str
        :param sourceuri: analysis sourceuri filter

        :type description: str
        :param description: analysis description

        :rtype: list of dict
        :return: Analysis information
        """

        # check if the organism exists
        res = self.session.query(self.model.analysis)
        if analysis_id:
            res = res.filter_by(analysis_id=analysis_id)
        if name:
            res = res.filter_by(name=name)
        if program:
            res = res.filter_by(program=program)
        if programversion:
            res = res.filter_by(programversion=programversion)
        if algorithm:
            res = res.filter_by(algorithm=algorithm)
        if sourcename:
            res = res.filter_by(sourcename=sourcename)
        if sourceversion:
            res = res.filter_by(sourceversion=sourceversion)
        if sourceuri:
            res = res.filter_by(sourceuri=sourceuri)
        if description:
            res = res.filter_by(description=description)

        data = []
        for ana in res:
            data.append({
                'analysis_id': ana.analysis_id,
                'name': ana.name,
                'program': ana.program,
                'programversion': ana.programversion,
                'algorithm': ana.algorithm,
                'sourcename': ana.sourcename,
                'sourceversion': ana.sourceversion,
                'sourceuri': ana.sourceuri,
                'description': ana.description,
                'timeexecuted': str(ana.timeexecuted),
            })
        return data

    def delete_analyses(self, analysis_id=None, name=None, program=None, programversion=None, algorithm=None, sourcename=None, sourceversion=None, sourceuri=None, description=None):
        """
        Delete analysis

        :type analysis_id: int
        :param analysis_id: analysis_id filter

        :type name: str
        :param name: analysis name filter

        :type program: str
        :param program: analysis program filter

        :type programversion: str
        :param programversion: analysis programversion filter

        :type algorithm: str
        :param algorithm: analysis algorithm filter

        :type sourcename: str
        :param sourcename: analysis sourcename filter

        :type sourceversion: str
        :param sourceversion: analysis sourceversion filter

        :type sourceuri: str
        :param sourceuri: analysis sourceuri filter

        :type description: str
        :param description: analysis description

        :rtype: None
        :return: None
        """

        res = self.session.query(self.model.analysis)
        if analysis_id:
            res = res.filter_by(analysis_id=analysis_id)
        if name:
            res = res.filter_by(name=name)
        if program:
            res = res.filter_by(program=program)
        if programversion:
            res = res.filter_by(programversion=programversion)
        if algorithm:
            res = res.filter_by(algorithm=algorithm)
        if sourcename:
            res = res.filter_by(sourcename=sourcename)
        if sourceversion:
            res = res.filter_by(sourceversion=sourceversion)
        if sourceuri:
            res = res.filter_by(sourceuri=sourceuri)
        if description:
            res = res.filter_by(description=description)

        res = res.delete()

        self.session.commit()
        return res
