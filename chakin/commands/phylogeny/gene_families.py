import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('gene_families')
@click.option(
    "--family_name",
    help="Restrict to families beginning with given prefix",
    type=str
)
@click.option(
    "--nuke",
    help="Removes all previous gene families data",
    is_flag=True
)
@pass_context
@custom_exception
@None_output
def cli(ctx, family_name="", nuke=False):
    """Adds an entry in the featureprop table in a chado database for each each family a gene belongs to (for use in https://github.com/legumeinfo/lis_context_viewer/).

Output:

    None
    """
    return ctx.gi.phylogeny.gene_families(family_name=family_name, nuke=nuke)
