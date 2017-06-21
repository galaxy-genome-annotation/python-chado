import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import chado_exception, None_output, _arg_split

@click.command('launch_docker_image')

@click.option(
    "--background",
    help="Launch the image in the background",
    is_flag=True
)

@pass_context
@chado_exception
@None_output
def cli(ctx, background=False):
    """Launch a chado docker image.

Output:

     None
        
    """
    return ctx.gi.util.launch_docker_image(background=background)
