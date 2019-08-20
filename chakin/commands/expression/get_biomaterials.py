import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, list_output


@click.command('get_biomaterials')
@click.option(
    "--provider_id",
    help="Limit query to the selected provider",
    type=int
)
@click.option(
    "--biomaterial_id",
    help="Limit query to the selected biomaterial id",
    type=int
)
@click.option(
    "--organism_id",
    help="Limit query to the selected organism",
    type=int
)
@click.option(
    "--biomaterial_name",
    help="Limit query to the selected biomaterial name",
    type=str
)
@click.option(
    "--analysis_id",
    help="Limit query to the selected analysis_id",
    type=int
)
@pass_context
@custom_exception
@list_output
def cli(ctx, provider_id="", biomaterial_id="", organism_id="", biomaterial_name="", analysis_id=""):
    """List biomaterials in the database

Output:

    List of biomaterials
    """
    return ctx.gi.expression.get_biomaterials(provider_id=provider_id, biomaterial_id=biomaterial_id, organism_id=organism_id, biomaterial_name=biomaterial_name, analysis_id=analysis_id)
