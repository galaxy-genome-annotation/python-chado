import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, str_output


@click.command('add_expression')
@click.argument("organism_id", type=str)
@click.argument("analysis_id", type=str)
@click.argument("file_path", type=str)
@click.option(
    "--separator",
    help="Separating character in the matrix file (ex : ','). Default character is tab.",
    default="	",
    show_default=True,
    type=str
)
@click.option(
    "--unit",
    help="The units associated with the loaded values (ie, FPKM, RPKM, raw counts)",
    type=str
)
@pass_context
@custom_exception
@str_output
def cli(ctx, organism_id, analysis_id, file_path, separator="	", unit=""):
    """Add an expression matrix file to the database

Output:

    Number of expression data loaded
    """
    return ctx.gi.expression.add_expression(organism_id, analysis_id, file_path, separator=separator, unit=unit)
