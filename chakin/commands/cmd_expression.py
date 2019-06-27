import click
from chakin.commands.expression.add_biomaterial import cli as add_biomaterial
from chakin.commands.expression.add_expression import cli as add_expression
from chakin.commands.expression.delete_all_biomaterials import cli as delete_all_biomaterials
from chakin.commands.expression.delete_biomaterial import cli as delete_biomaterial
from chakin.commands.expression.delete_biomaterials import cli as delete_biomaterials
from chakin.commands.expression.get_biomaterials import cli as get_biomaterials


@click.group()
def cli():
    """
    Interact with expressions
    """
    pass


cli.add_command(add_biomaterial)
cli.add_command(add_expression)
cli.add_command(delete_all_biomaterials)
cli.add_command(delete_biomaterial)
cli.add_command(delete_biomaterials)
cli.add_command(get_biomaterials)
