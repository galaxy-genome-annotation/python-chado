import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, None_output


@click.command('launch_docker_image')
@click.option(
    "--background",
    help="Launch the image in the background",
    is_flag=True
)
@click.option(
    "--no_yeast",
    help="Disable loading of example yeast data",
    is_flag=True
)
@pass_context
@custom_exception
@None_output
def cli(ctx, background=False, no_yeast=False):
    """Launch a chado docker image.

Output:

    None
    """
    return ctx.gi.util.launch_docker_image(background=background, no_yeast=no_yeast)
