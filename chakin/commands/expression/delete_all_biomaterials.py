import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, None_output


@click.command('delete_all_biomaterials')
@click.option(
    "--confirm",
    help="Confirm that you really do want to delete ALL of the biomaterials.",
    is_flag=True
)
@pass_context
@custom_exception
@None_output
def cli(ctx, confirm=False):
    """Delete all biomaterials

Output:

    None
    """
    return ctx.gi.expression.delete_all_biomaterials(confirm=confirm)
