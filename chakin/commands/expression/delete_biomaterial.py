import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, str_output


@click.command('delete_biomaterial')
@click.option(
    "--names",
    help="JSON list of biomaterial names to delete.",
    default="[]",
    show_default=True,
    type=str
)
@click.option(
    "--ids",
    help="JSON list of biomaterial ids to delete.",
    default="[]",
    show_default=True,
    type=str
)
@click.option(
    "--organism_id",
    help="Delete all biomaterial associated with this organism id.",
    type=str
)
@click.option(
    "--analysis_id",
    help="Delete all biomaterial associated with this analysis id.",
    type=str
)
@pass_context
@custom_exception
@str_output
def cli(ctx, names="[]", ids="[]", organism_id="", analysis_id=""):
    """

Output:

    I have no idea
    """
    return ctx.gi.expression.delete_biomaterial(names=names, ids=ids, organism_id=organism_id, analysis_id=analysis_id)
