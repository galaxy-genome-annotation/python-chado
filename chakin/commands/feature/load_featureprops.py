import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('load_featureprops')
@click.argument("tab_file", type=str)
@click.argument("analysis_id", type=int)
@click.argument("organism_id", type=int)
@click.argument("prop_type", type=str)
@click.option(
    "--feature_type",
    help="Type of the target features in sequence ontology (will speed up loading if specified)",
    type=str
)
@click.option(
    "--match_on_name",
    help="Match features using their name instead of their uniquename",
    is_flag=True
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, tab_file, analysis_id, organism_id, prop_type, feature_type="", match_on_name=False):
    """Load feature properties from a tabular file (Column1: feature name or uniquename, Column2: property value)

Output:

    Number of inserted featureprop
    """
    return ctx.gi.feature.load_featureprops(tab_file, analysis_id, organism_id, prop_type, feature_type=feature_type, match_on_name=match_on_name)
