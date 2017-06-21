import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, list_output, _arg_split

@click.command('get_organisms')

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
@chado_exception
@list_output
def cli(ctx, genus="", species="", common="", abbr="", comment=""):
    """Get all organisms

Output:

     Organisms information
        
    """
    return ctx.gi.organism.get_organisms(genus=genus, species=species, common=common, abbr=abbr, comment=comment)
