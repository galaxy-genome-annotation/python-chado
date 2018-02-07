import click
from chakin.commands.feature.load_fasta import cli as func0
from chakin.commands.feature.load_gff import cli as func1


@click.group()
def cli():
    """
    Access to the chado features
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
