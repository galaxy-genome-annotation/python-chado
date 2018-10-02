import click
from chakin.commands.expression.add_biomaterial import cli as func0
from chakin.commands.expression.get_biomaterials import cli as func1


@click.group()
def cli():
    """
    Interact with expressions
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
