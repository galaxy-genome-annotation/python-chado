"""
Contains possible interactions with the Chado Phylogeny Module
http://gmod.org/wiki/Chado_Phylogeny_Module
As implemented in https://github.com/legumeinfo/tripal_phylotree
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import os.path
from Bio import Phylo
from chado.client import Client
from future import standard_library

standard_library.install_aliases()


class PhylogenyClient(Client):
    """
    Access to the chado phylogeny content
    """

    def load_tree(self, newick, analysis_id, name=None, xref_db='null', xref_accession=None, match_on_name=False, prefix=""):
        """
        Load a phylogenetic tree (Newick format) into Chado db

        :type newick: str
        :param newick: Newick file to load (or directory containing multiple newick files to load)

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type name: str
        :param name: The name given to the phylotree entry in the database (default=<filename>)

        :type xref_db: str
        :param xref_db: The name of the db to link dbxrefs for the trees (default: "null")

        :type xref_accession: str
        :param xref_accession: The accession to use for dbxrefs for the trees (assumed same as name unless otherwise specified)

        :type match_on_name: bool
        :param match_on_name: Match polypeptide features usnig their name instead of their uniquename

        :type prefix: str
        :param prefix: Comma-separated list of prefix to be removed from identifiers (e.g species prefixes when using loading OrthoFinder output)

        :rtype: None
        :return: None
        """

        if not os.path.exists(newick):
            raise Exception("Could not read input file/dir '{}'".format(newick))

        if not os.path.isdir(newick):
            return self._load_single_tree(newick, analysis_id, name, xref_db, xref_accession, match_on_name)

        out = []
        for nf in os.listdir(newick):
            name = os.path.splitext(os.path.basename(nf))[0]
            if name.endswith('_tree'):  # OrthoFinder file
                name = name[:-5]
            print("Loading newick '{}' from file '{}'".format(name, nf))
            out.append(self._load_single_tree(newick + nf, analysis_id, name, xref_db, xref_accession, match_on_name))

        return out

    def _load_single_tree(self, newick, analysis_id, name=None, xref_db='null', xref_accession=None, match_on_name=False, prefix=""):
        """
        Load a phylogenetic tree (Newick format) into Chado db

        :type newick: str
        :param newick: Newick file to load

        :type analysis_id: int
        :param analysis_id: Analysis ID

        :type name: str
        :param name: The name given to the phylotree entry in the database (default=<filename>)

        :type xref_db: str
        :param xref_db: The name of the db to link dbxrefs for the trees (default: "null")

        :type xref_accession: str
        :param xref_accession: The accession to use for dbxrefs for the trees (assumed same as name unless otherwise specified)

        :type match_on_name: bool
        :param match_on_name: Match polypeptide features usnig their name instead of their uniquename

        :type prefix: str
        :param prefix: Comma-separated list of prefix to be removed from identifiers (e.g species prefixes when using loading OrthoFinder output)

        :rtype: None
        :return: None
        """

        cvterms = self.add_cvterms()
        for cvt in cvterms:
            cvterms[cvt] = self.session.query(self.model.cvterm).filter_by(name=cvterms[cvt]['name'], cv_id=cvterms[cvt]['cv_id']).one()

        res = self.session.query(self.model.db).filter_by(name=xref_db)
        if res.count() == 0:
            raise Exception("Could not find xref_db '{}'".format(xref_db))
        db = res.one()

        res = self.session.query(self.model.analysis).filter_by(analysis_id=analysis_id)
        if res.count() == 0:
            raise Exception("Could not find analysis with id {}".format(analysis_id))
        analysis = res.one()

        if not os.path.exists(newick):
            raise Exception("Could not read input file '{}'".format(newick))

        tree_file = open(newick, 'r')

        if not name:
            name = os.path.splitext(os.path.basename(newick))[0]
            if name.endswith('_tree'):  # OrthoFinder file
                name = name[:-5]

        # Add phylotree and its dbxref
        res = self.session.query(self.model.dbxref).filter_by(accession=name, db=db)
        if res.count() > 0:
            print("Found existing '{}' dbxref, using it".format(name))
            dbxref = res.one()
        else:
            dbxref = self.model.dbxref()
            dbxref.accession = name
            dbxref.db = db

            self.session.add(dbxref)

        res = self.session.query(self.model.phylotree).filter_by(name=name)
        if res.count() > 0:
            print("Found existing '{}' phylotree, using it".format(name))
            db_tree = res.one()
        else:
            db_tree = self.model.phylotree()
            db_tree.name = name
            db_tree.dbxref = dbxref
            db_tree.comment = tree_file.read()
            db_tree.analysis = analysis

            self.session.add(db_tree)

            tree_file.seek(0)

        # Read input newick file
        trees = Phylo.parse(tree_file, 'newick')

        # Retrieve leaf features
        res = self.session.query(self.model.feature, self.model.cvterm).filter(self.model.cvterm.name == 'polypeptide').all()
        if match_on_name:
            peps = {p.feature.uniquename: p.feature for p in res}
        else:
            peps = {p.feature.name: p.feature for p in res}

        prefixes = prefix.split(',')

        indexes = []
        for tree in trees:

            self._create_nodes(tree.root, 0, db_tree, cvterms, peps, indexes, prefixes=prefixes)

        tree_file.close()

        self.session.commit()

        return {
            'phylotree_id': db_tree.phylotree_id,
            'dbxref_id': db_tree.dbxref.dbxref_id,
            'name': db_tree.name,
            'type_id': db_tree.type_id,
            'analysis_id': db_tree.analysis.analysis_id,
            'comment': db_tree.comment,
        }

    def _create_nodes(self, clade, depth, db_tree, cvterms, peps, indexes, parent_node=None, prefixes=[]):
        """
        Recursive function which will insert into db all nodes of a phylo tree
        """

        lindex = 1
        while lindex in indexes:
            lindex += 1
        rindex = lindex + 2 * (len(list(clade.find_clades())) - 1) + 1
        indexes.append(lindex)
        indexes.append(rindex)

        if depth == 0:
            term = cvterms['phylo_root']
        elif clade.is_terminal():
            term = cvterms['phylo_leaf']
        else:
            term = cvterms['phylo_interior']

        node = self.model.phylonode()
        node.phylotree = db_tree
        node.phylonode = parent_node
        node.left_idx = lindex
        node.right_idx = rindex
        node.cvterm = term
        if clade.name:
            cname = clade.name
            # Remove prefix from id (typically added by orthofinder)
            for p in prefixes:
                if cname.startswith(p):
                    cname = cname[len(p) - 1:]
                elif cname.startswith(p + '_'):
                    cname = cname[len(p):]
            if cname not in peps:
                raise Exception("Could not find polypeptide '{}', rolling back".format(cname))
            node.feature = peps[cname]
        node.label = clade.name
        node.distance = clade.branch_length
        self.session.add(node)

        # print('--' * depth + 'lindex: '+str(lindex)+' rindex: '+str(rindex)+' branch_len:'+str(clade.branch_length)+' name:'+str(clade.name)+' term:'+str(term.name)+' parent_li:'+str(node.phylonode.left_idx if node.phylonode else None))

        clades = clade.clades
        for c in clades:
            self._create_nodes(c, depth + 1, db_tree, cvterms, peps, indexes, node, prefixes=[])

    def add_cvterms(self):
        """
        Make sure required cvterms are loaded

        :rtype: list of dict
        :return: created cvterms
        """

        # check if the db exists
        res = self.session.query(self.model.db).filter_by(name='tripal')

        if res.count() > 0:
            db = res.one()
        else:
            db = self.model.db()
            db.name = 'tripal'
            db.definition = 'Used as a database placeholder for tripal defined objects such as tripal cvterms'

            self.session.add(db)

        # check if the cv exists
        res = self.session.query(self.model.cv).filter_by(name='tripal_phylotree')

        if res.count() > 0:
            cv = res.one()
        else:
            cv = self.model.cv()
            cv.name = 'tripal_phylotree'
            cv.definition = 'Terms used by the Tripal phylotree module for phylogenetic and taxonomic trees.'

            self.session.add(cv)

        terms = {
            'phylo_interior': 'An interior node in a phylogenetic tree.',
            'phylo_root': 'The root node of a phylogenetic tree.',
            'phylo_leaf': 'A leaf node in a phylogenetic tree.',
        }

        out = {}

        for term in terms:
            res = self.session.query(self.model.dbxref).filter_by(accession=term, db=db)
            if res.count() > 0:
                dbxref = res.one()
            else:
                dbxref = self.model.dbxref()
                dbxref.accession = term
                dbxref.db = db

                self.session.add(dbxref)

            res = self.session.query(self.model.cvterm).filter_by(name=term, cv=cv)
            if res.count() > 0:
                cvterm = res.one()
            else:
                cvterm = self.model.cvterm()
                cvterm.name = term
                cvterm.cv = cv
                cvterm.definition = terms[term]
                cvterm.dbxref = dbxref

                self.session.add(cvterm)

            out[cvterm.name] = {
                'name': cvterm.name,
                'cv_id': cvterm.cv.cv_id,
                'definition': cvterm.definition,
                'dbxref_id': cvterm.dbxref.dbxref_id
            }

        return out
