import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, list_output


@click.command('get_biomaterials')
@click.option(
    "--provider_id",
    help="Limit query to the selected provider",
    type=str
)
@click.option(
    "--biomaterial_id",
    help="Limit query to the selected biomaterial id",
    type=str
)
@click.option(
    "--organism_id",
    help="Limit query to the selected organism",
    type=str
)
@click.option(
    "--dbxref_id",
    help="Limit query to the selected ref",
    type=str
)
@click.option(
    "--biomaterial_name",
    help="Limit query to the selected biomaterial name",
    type=str
)
@pass_context
@custom_exception
@list_output
def cli(ctx, provider_id="", biomaterial_id="", organism_id="", dbxref_id="", biomaterial_name=""):
    """List biomaterials in the database

Output:

    List of biomaterials
    """
    return ctx.gi.expression.get_biomaterials(provider_id=provider_id, biomaterial_id=biomaterial_id, organism_id=organism_id, dbxref_id=dbxref_id, biomaterial_name=biomaterial_name)
