"""
Contains possible interactions with the Chado Phylogeny Module
http://gmod.org/wiki/Chado_Phylogeny_Module
As implemented in https://github.com/legumeinfo/tripal_phylotree
and https://github.com/legumeinfo/lis_context_viewer/
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os.path
import warnings

from Bio import Phylo

import chado
from chado.client import Client

from future import standard_library

from sqlalchemy import Column, ForeignKey, Index, Integer, String, Table, UniqueConstraint
from sqlalchemy import exc as sa_exc
from sqlalchemy.orm import aliased

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
        :param match_on_name: Match polypeptide features using their name instead of their uniquename

        :type prefix: str
        :param prefix: Comma-separated list of prefix to be removed from identifiers (e.g species prefixes when using loading OrthoFinder output)

        :rtype: dict
        :return: Number of inserted trees
        """

        if not os.path.exists(newick):
            raise Exception("Could not read input file/dir '{}'".format(newick))

        if not os.path.isdir(newick):
            return self._load_single_tree(newick, analysis_id, name, xref_db, xref_accession, match_on_name, prefix)

        # prefetch this for performances
        peps = self._fetch_peptides(match_on_name)

        out = []
        files = os.listdir(newick)
        cur = 1
        total = len(files)
        for nf in sorted(files):
            name = os.path.splitext(os.path.basename(nf))[0]
            if name.endswith('_tree'):  # OrthoFinder file
                name = name[:-5]
            print("Loading newick '{}' from file '{}' {}/{}".format(name, os.path.join(newick, nf), cur, total))
            out.append(self._load_single_tree(os.path.join(newick, nf), analysis_id, name, xref_db, xref_accession, match_on_name, prefix, peps, commit=False))
            cur += 1

        self.session.commit()

        return {'inserted': cur}

    def _load_single_tree(self, newick, analysis_id, name=None, xref_db='null', xref_accession=None, match_on_name=False, prefix="", peps=None, commit=True):
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
        :param match_on_name: Match polypeptide features using their name instead of their uniquename

        :type prefix: str
        :param prefix: Comma-separated list of prefix to be removed from identifiers (e.g species prefixes when using loading OrthoFinder output)

        :rtype: dict
        :return: Tree description
        """

        self.add_cvterms()

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
            dbxref = res.one()
        else:
            dbxref = self.model.dbxref()
            dbxref.accession = name
            dbxref.db = db

            self.session.add(dbxref)

        res = self.session.query(self.model.phylotree).filter_by(name=name)
        if res.count() > 0:
            print("Found existing '{}' phylotree, updating it".format(name))
            db_tree = res.one()
            self.session.query(self.model.phylonode).filter_by(phylotree_id=db_tree.phylotree_id).delete()
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
        if not peps:
            peps = self._fetch_peptides(match_on_name)

        prefixes = prefix.split(',')

        indexes = []
        for tree in trees:

            self._create_nodes(tree.root, 0, db_tree, peps, indexes, prefixes=prefixes)

        tree_file.close()

        if commit:
            self.session.commit()

        return {
            'phylotree_id': db_tree.phylotree_id,
            'dbxref_id': db_tree.dbxref.dbxref_id,
            'name': db_tree.name,
            'type_id': db_tree.type_id,
            'analysis_id': db_tree.analysis.analysis_id,
            'comment': db_tree.comment,
        }

    def _fetch_peptides(self, match_on_name):

        ppterm = self.ci.get_cvterm_id('polypeptide', 'sequence')
        res = self.session.query(self.model.feature.feature_id, self.model.feature.name, self.model.feature.uniquename) \
            .filter_by(type_id=ppterm) \
            .all()

        if match_on_name:
            peps = {p.name: p.feature_id for p in res}
        else:
            peps = {p.uniquename: p.feature_id for p in res}

        return peps

    def _create_nodes(self, clade, depth, db_tree, peps, indexes, parent_node=None, prefixes=[]):
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
            term = self.ci.get_cvterm_id('phylo_root', 'tripal_phylotree')
        elif clade.is_terminal():
            term = self.ci.get_cvterm_id('phylo_leaf', 'tripal_phylotree')
        else:
            term = self.ci.get_cvterm_id('phylo_interior', 'tripal_phylotree')

        node = self.model.phylonode()
        node.phylotree = db_tree
        node.phylonode = parent_node
        node.left_idx = lindex
        node.right_idx = rindex
        node.cvterm_id = term
        if clade.name:
            cname = clade.name
            # Remove prefix from id (typically added by orthofinder)
            for p in prefixes:
                if cname.startswith(p + '_'):
                    cname = cname[len(p) + 1:]
                    break
                elif cname.startswith(p):
                    cname = cname[len(p):]
                    break
            if cname not in peps:
                raise Exception("Could not find polypeptide '{}', rolling back".format(cname))
            node.feature_id = peps[cname]
        node.label = clade.name
        node.distance = clade.branch_length
        self.session.add(node)

        # print('--' * depth + 'lindex: '+str(lindex)+' rindex: '+str(rindex)+' branch_len:'+str(clade.branch_length)+' name:'+str(clade.name)+' term:'+str(term.name)+' parent_li:'+str(node.phylonode.left_idx if node.phylonode else None))

        clades = clade.clades
        for c in clades:
            self._create_nodes(c, depth + 1, db_tree, peps, indexes, node, prefixes=prefixes)

    def add_cvterms(self):
        """
        Make sure required cvterms are loaded

        :rtype: list of dict
        :return: created cvterms
        """

        terms = {
            'phylo_interior': 'An interior node in a phylogenetic tree.',
            'phylo_root': 'The root node of a phylogenetic tree.',
            'phylo_leaf': 'A leaf node in a phylogenetic tree.',
        }

        out = {}

        for term in terms:
            try:
                cvterm_id = self.ci.get_cvterm_id(term, 'tripal_phylotree')
            except chado.RecordNotFoundError:
                cvterm_id = self.ci.create_cvterm(term=term, term_definition=terms[term], cv='tripal_phylotree', db='tripal', cv_definition="Terms used by the Tripal phylotree module for phylogenetic and taxonomic trees.", db_definition="Used as a database placeholder for tripal defined objects such as tripal cvterms")

            out[term] = {
                'cvterm_id': cvterm_id,
                'name': term,
                'definition': terms[term],
            }

        # Some cvterm for LIS-GCV
        try:
            cvterm_id = self.ci.get_cvterm_id('gene family', 'GCV_properties')
        except chado.RecordNotFoundError:
            cvterm_id = self.ci.create_cvterm(term='gene family', term_definition='A group of genes presumed to be related by common ancestry', cv='GCV_properties', db='null', cv_definition="Used by https://github.com/legumeinfo/lis_context_viewer/", db_definition="A fake database for local items")

        out['gene family'] = {
            'cvterm_id': cvterm_id,
            'name': 'gene family',
            'definition': 'A group of genes presumed to be related by common ancestry',
        }

        return out

    def gene_order(self, nuke=False):
        """
        Orders all the genes in the database by their order on their respective chromosomes in the gene_order table (for use in https://github.com/legumeinfo/lis_context_viewer/).

        :type nuke: bool
        :param nuke: Removes all previous gene ordering data

        :rtype: None
        :return: None
        """

        if not hasattr(self.model, 'gene_order'):
            # Create gene_order table if it doesn't exist yet
            gene_order_table = Table(
                'gene_order', self.metadata,
                Column('gene_order_id', Integer, primary_key=True),
                Column('chromosome_id', Integer, ForeignKey(self.model.feature.feature_id)),
                Column('gene_id', Integer, ForeignKey(self.model.feature.feature_id)),
                Column('number', Integer, nullable=False),
                UniqueConstraint('chromosome_id', 'number', name='gene_order_c1'),
                UniqueConstraint('gene_id', name='gene_order_gene_id_key'),
                schema=self.ci.dbschema
            )
            gene_order_table.create(self.engine)

            # Reload the db schema
            with warnings.catch_warnings():
                # https://stackoverflow.com/a/5225951
                warnings.simplefilter("ignore", category=sa_exc.SAWarning)
                self.ci._reflect_tables()
                self.model = self.ci.model

        if nuke:
            # Remove old ordering
            self.session.query(self.model.gene_order).delete()

        # Insert ordering info
        chroterm = self.ci.get_cvterm_id('chromosome', 'sequence')
        contigterm = self.ci.get_cvterm_id('contig', 'sequence')
        scontigterm = self.ci.get_cvterm_id('supercontig', 'sequence')
        geneterm = self.ci.get_cvterm_id('gene', 'sequence')

        chromosomes = self.session.query(self.model.feature.feature_id) \
            .filter((self.model.feature.type_id == chroterm) | (self.model.feature.type_id == contigterm) | (self.model.feature.type_id == scontigterm)) \
            .all()

        existing = self.session.query(self.model.gene_order) \
            .all()

        existing = {'{}_{}'.format(ex.chromosome_id, ex.gene_id): ex.number for ex in existing}

        for chro in chromosomes:
            genes = self.session.query(self.model.feature.feature_id) \
                .join(self.model.featureloc, self.model.feature.feature_id == self.model.featureloc.feature_id) \
                .filter(self.model.feature.type_id == geneterm, self.model.featureloc.srcfeature_id == chro.feature_id) \
                .order_by(self.model.featureloc.fmin.asc()) \
                .all()

            pos = 0
            for g in genes:
                pos += 1
                dup_check_id = '{}_{}'.format(chro.feature_id, g.feature_id)
                if dup_check_id in existing:
                    if existing[dup_check_id] == pos:
                        continue  # We already have this order record, skip it
                    else:
                        raise Exception("Found an existing gene_order record with different value. Rolling back, use the --nuke option to replace all existing values.")

                order = self.model.gene_order()
                order.chromosome_id = chro.feature_id
                order.gene_id = g.feature_id
                order.number = pos
                self.session.add(order)

        self.session.commit()

    def gene_families(self, family_name='', nuke=False):
        """
        Adds an entry in the featureprop table in a chado database for each each family a gene belongs to (for use in https://github.com/legumeinfo/lis_context_viewer/).

        :type family_name: str
        :param family_name: Restrict to families beginning with given prefix

        :type nuke: bool
        :param nuke: Removes all previous gene families data

        :rtype: None
        :return: None
        """

        if not hasattr(self.model, 'gene_family_assignment'):
            # Create gene_order table if it doesn't exist yet
            gfa_table = Table(
                'gene_family_assignment', self.metadata,
                Column('gene_family_assignment_id', Integer, primary_key=True),
                Column('gene_id', Integer, ForeignKey(self.model.feature.feature_id)),
                Column('family_label', String, nullable=False),
                Index('gene_family_assignment_idx1', 'family_label'),
                schema=self.ci.dbschema
            )
            gfa_table.create(self.engine)

            # Reload the db schema
            with warnings.catch_warnings():
                # https://stackoverflow.com/a/5225951
                warnings.simplefilter("ignore", category=sa_exc.SAWarning)
                self.ci._reflect_tables()
                self.model = self.ci.model

        self.add_cvterms()

        gfterm = self.ci.get_cvterm_id('gene family', 'GCV_properties')

        if nuke:
            # Remove old ordering
            self.session.query(self.model.gene_family_assignment).delete()
            self.session.query(self.model.featureprop).filter_by(type_id=gfterm).delete()

        # Insert ordering info
        trees = self.session.query(self.model.phylotree) \
            .filter(self.model.phylotree.name != 'NCBI taxonomy tree')

        if family_name:
            trees.filter(self.model.phylotree.name.like(family_name + '%'))

        trees = trees.all()

        mrna_alias = aliased(self.model.feature)
        gene_alias = aliased(self.model.feature)
        rel_alias = aliased(self.model.feature_relationship)

        assignements = {}

        for tree in trees:
            genes = self.session.query(self.model.phylotree.name, gene_alias.feature_id) \
                .join(self.model.phylonode, self.model.phylonode.phylotree_id == self.model.phylotree.phylotree_id) \
                .join(self.model.feature_relationship, self.model.phylonode.feature_id == self.model.feature_relationship.subject_id) \
                .join(mrna_alias, self.model.feature_relationship.object_id == mrna_alias.feature_id) \
                .join(rel_alias, mrna_alias.feature_id == rel_alias.subject_id) \
                .join(gene_alias, rel_alias.object_id == gene_alias.feature_id) \
                .filter(self.model.phylonode.phylotree_id == tree.phylotree_id) \
                .all()

            for gene in genes:
                if gene.feature_id not in assignements:
                    assignements[gene.feature_id] = []

                assignements[gene.feature_id].append(gene.name)

        # Cache all existing gene_family_assignment rows
        existing_gfa = self.session.query(self.model.gene_family_assignment) \
            .all()
        existing_gfa = [(ex.gene_id, ex.family_label) for ex in existing_gfa]

        # Cache all existing featureprop rows
        existing_fp = self.session.query(self.model.featureprop) \
            .filter_by(type_id=gfterm) \
            .all()
        existing_fp = {(ex.feature_id, ex.value): ex.rank for ex in existing_fp}

        for gene in assignements:

            rank = 0
            tree_list = ','.join(assignements[gene])
            if (gene, tree_list) in existing_fp:
                rank = existing_fp[(gene, tree_list)] + 1

            fprop = self.model.featureprop()
            fprop.feature_id = gene
            fprop.type_id = gfterm
            fprop.value = tree_list
            fprop.rank = rank
            self.session.add(fprop)

            for tree in assignements[gene]:
                if (gene, tree) not in existing_gfa:
                    gfa = self.model.gene_family_assignment()
                    gfa.gene_id = gene
                    gfa.family_label = tree
                    self.session.add(gfa)

        self.session.commit()
