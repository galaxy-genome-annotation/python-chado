import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, list_output


@click.command('delete_features')
@click.option(
    "--organism_id",
    help="organism_id filter",
    type=str
)
@click.option(
    "--analysis_id",
    help="analysis_id filter",
    type=str
)
@click.option(
    "--name",
    help="name filter",
    type=str
)
@click.option(
    "--uniquename",
    help="uniquename filter",
    type=str
)
@pass_context
@custom_exception
@list_output
def cli(ctx, organism_id="", analysis_id="", name="", uniquename=""):
    """Get all or some features

Output:

    Features information
    """
    return ctx.gi.feature.delete_features(organism_id=organism_id, analysis_id=analysis_id, name=name, uniquename=uniquename)
