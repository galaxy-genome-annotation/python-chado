import click
from chakin.commands.feature.delete_features import cli as func0
from chakin.commands.feature.get_feature_analyses import cli as func1
from chakin.commands.feature.get_feature_cvterms import cli as func2
from chakin.commands.feature.get_features import cli as func3
from chakin.commands.feature.load_fasta import cli as func4
from chakin.commands.feature.load_featureprops import cli as func5
from chakin.commands.feature.load_gff import cli as func6
from chakin.commands.feature.load_go import cli as func7


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
cli.add_command(func4)
cli.add_command(func5)
cli.add_command(func6)
cli.add_command(func7)
