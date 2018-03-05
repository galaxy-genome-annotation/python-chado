import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('add_organism')
@click.argument("genus", type=str)
@click.argument("common", type=str)
@click.argument("abbr", type=str)
@click.option(
    "--species",
    help="The species of the organism",
    type=str
)
@click.option(
    "--comment",
    help="A comment / description",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, genus, common, abbr, species="", comment=""):
    """Add a new organism to the Chado database

Output:

    Organism information
    """
    return ctx.gi.organism.add_organism(genus, common, abbr, species=species, comment=comment)
