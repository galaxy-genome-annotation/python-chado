import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('delete_organisms')
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
@None_output
def cli(ctx, organism_id="", genus="", species="", common="", abbr="", comment=""):
    """Delete all organisms

Output:

    None
    """
    return ctx.gi.organism.delete_organisms(organism_id=organism_id, genus=genus, species=species, common=common, abbr=abbr, comment=comment)
