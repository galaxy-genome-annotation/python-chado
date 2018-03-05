import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, list_output


@click.command('get_organisms')
@click.option(
    "--organism_id",
    help="organism_id filter",
    type=int
)
@click.option(
    "--genus",
    help="genus filter",
    type=str
)
@click.option(
    "--species",
    help="species filter",
    type=str
)
@click.option(
    "--common",
    help="common filter",
    type=str
)
@click.option(
    "--abbr",
    help="abbr filter",
    type=str
)
@click.option(
    "--comment",
    help="comment filter",
    type=str
)
@pass_context
@custom_exception
@list_output
def cli(ctx, organism_id="", genus="", species="", common="", abbr="", comment=""):
    """Get all or some organisms

Output:

    Organisms information
    """
    return ctx.gi.organism.get_organisms(organism_id=organism_id, genus=genus, species=species, common=common, abbr=abbr, comment=comment)
