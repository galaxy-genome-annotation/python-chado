import click
from chakin.commands.expression.add_biomaterial import cli as func0
from chakin.commands.expression.add_expression import cli as func1
from chakin.commands.expression.delete_all_biomaterials import cli as func2
from chakin.commands.expression.delete_biomaterial import cli as func3
from chakin.commands.expression.delete_biomaterials import cli as func4
from chakin.commands.expression.get_biomaterials import cli as func5


@click.group()
def cli():
    """
    Interact with expressions
    """
    pass


cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
cli.add_command(func3)
cli.add_command(func4)
cli.add_command(func5)
