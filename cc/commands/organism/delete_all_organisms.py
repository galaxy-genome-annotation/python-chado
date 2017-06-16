import click
from cc.cli import pass_context, json_loads
from cc.decorators import chado_exception, None_output, _arg_split

@click.command('delete_all_organisms')

@click.option(
    "--confirm",
    help="Confirm that you really do want to delete ALL of the organisms.",
    is_flag=True
)

@pass_context
@chado_exception
@None_output
def cli(ctx, confirm=False):
    """Get all organisms

Output:

     None
        
    """
    return ctx.gi.organism.delete_all_organisms(confirm=confirm)
