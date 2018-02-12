import click
from chakin.commands.phylogeny.add_cvterms import cli as func0
from chakin.commands.phylogeny.gene_families import cli as func1
from chakin.commands.phylogeny.gene_order import cli as func2
from chakin.commands.phylogeny.load_tree import cli as func3


@click.group()
def cli():
    """
    Access to the chado phylogeny content
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
cli.add_command(func3)
