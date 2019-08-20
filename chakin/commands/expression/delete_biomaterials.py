import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, str_output


@click.command('delete_biomaterials')
@click.option(
    "--names",
    help="JSON list of biomaterial names to delete.",
    type=str
)
@click.option(
    "--ids",
    help="JSON list of biomaterial ids to delete.",
    type=str
)
@click.option(
    "--organism_id",
    help="Delete all biomaterial associated with this organism id.",
    type=int
)
@click.option(
    "--analysis_id",
    help="Delete all biomaterial associated with this analysis id.",
    type=int
)
@pass_context
@custom_exception
@str_output
def cli(ctx, names="", ids="", organism_id="", analysis_id=""):
    """Will delete biomaterials based on selector. Only one selector will be used.

Output:

    Number of deleted biomaterials
    """
    return ctx.gi.expression.delete_biomaterials(names=names, ids=ids, organism_id=organism_id, analysis_id=analysis_id)
