import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, dict_output, _arg_split

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
@chado_exception
@dict_output
def cli(ctx, genus, common, abbr, species="", comment=""):
    """Add a new organism to the Chado database

Output:

     Organism information
        
    """
    return ctx.gi.organism.add_organism(genus, common, abbr, species=species, comment=comment)
