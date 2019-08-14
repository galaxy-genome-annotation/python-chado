import click
from chakin.commands.feature.delete_features import cli as delete_features
from chakin.commands.feature.get_feature_analyses import cli as get_feature_analyses
from chakin.commands.feature.get_feature_cvterms import cli as get_feature_cvterms
from chakin.commands.feature.get_features import cli as get_features
from chakin.commands.feature.load_fasta import cli as load_fasta
from chakin.commands.feature.load_featureprops import cli as load_featureprops
from chakin.commands.feature.load_gff import cli as load_gff
from chakin.commands.feature.load_go import cli as load_go


@click.group()
def cli():
    """
    Access to the chado features
    """
    pass


cli.add_command(delete_features)
cli.add_command(get_feature_analyses)
cli.add_command(get_feature_cvterms)
cli.add_command(get_features)
cli.add_command(load_fasta)
cli.add_command(load_featureprops)
cli.add_command(load_gff)
cli.add_command(load_go)
