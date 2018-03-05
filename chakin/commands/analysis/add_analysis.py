import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('add_analysis')
@click.argument("name", type=str)
@click.argument("program", type=str)
@click.argument("programversion", type=str)
@click.argument("sourcename", type=str)
@click.option(
    "--algorithm",
    help="analysis algorithm",
    type=str
)
@click.option(
    "--sourceversion",
    help="analysis sourceversion",
    type=str
)
@click.option(
    "--sourceuri",
    help="analysis sourceuri",
    type=str
)
@click.option(
    "--description",
    help="analysis description",
    type=str
)
@click.option(
    "--date_executed",
    help="analysis date_executed (yyyy-mm-dd)",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, name, program, programversion, sourcename, algorithm="", sourceversion="", sourceuri="", description="", date_executed=""):
    """Create an analysis

Output:

    Analysis information
    """
    return ctx.gi.analysis.add_analysis(name, program, programversion, sourcename, algorithm=algorithm, sourceversion=sourceversion, sourceuri=sourceuri, description=description, date_executed=date_executed)
