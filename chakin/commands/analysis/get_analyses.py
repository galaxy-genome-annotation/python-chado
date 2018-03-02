import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, list_output


@click.command('get_analyses')
@click.option(
    "--analysis_id",
    help="analysis_id filter",
    type=int
)
@click.option(
    "--name",
    help="analysis name filter",
    type=str
)
@click.option(
    "--program",
    help="analysis program filter",
    type=str
)
@click.option(
    "--programversion",
    help="analysis programversion filter",
    type=str
)
@click.option(
    "--algorithm",
    help="analysis algorithm filter",
    type=str
)
@click.option(
    "--sourcename",
    help="analysis sourcename filter",
    type=str
)
@click.option(
    "--sourceversion",
    help="analysis sourceversion filter",
    type=str
)
@click.option(
    "--sourceuri",
    help="analysis sourceuri filter",
    type=str
)
@click.option(
    "--description",
    help="analysis description",
    type=str
)
@pass_context
@custom_exception
@list_output
def cli(ctx, analysis_id="", name="", program="", programversion="", algorithm="", sourcename="", sourceversion="", sourceuri="", description=""):
    """Get all or some analyses

Output:

    Analysis information
    """
    return ctx.gi.analysis.get_analyses(analysis_id=analysis_id, name=name, program=program, programversion=programversion, algorithm=algorithm, sourcename=sourcename, sourceversion=sourceversion, sourceuri=sourceuri, description=description)
