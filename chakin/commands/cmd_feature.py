import click
from chakin.commands.feature.delete_features import cli as func0
from chakin.commands.feature.get_features import cli as func1
from chakin.commands.feature.load_fasta import cli as func2
from chakin.commands.feature.load_gff import cli as func3


@click.group()
def cli():
    """
    Access to the chado features
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
cli.add_command(func3)
