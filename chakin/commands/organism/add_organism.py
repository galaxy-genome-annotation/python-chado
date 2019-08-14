import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output


@click.command('add_organism')
@click.argument("genus", type=str)
@click.argument("species", type=str)
@click.argument("common", type=str)
@click.argument("abbr", type=str)
@click.option(
    "--comment",
    help="A comment / description",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, genus, species, common, abbr, comment=""):
    """Add a new organism to the Chado database

Output:

    Organism information
    """
    return ctx.gi.organism.add_organism(genus, species, common, abbr, comment=comment)
