import click
from chakin.commands.phylogeny.add_cvterms import cli as add_cvterms
from chakin.commands.phylogeny.gene_families import cli as gene_families
from chakin.commands.phylogeny.gene_order import cli as gene_order
from chakin.commands.phylogeny.load_tree import cli as load_tree


@click.group()
def cli():
    """
    Access to the chado phylogeny content
    """
    pass


cli.add_command(add_cvterms)
cli.add_command(gene_families)
cli.add_command(gene_order)
cli.add_command(load_tree)
