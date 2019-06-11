import click
from chakin.commands.load.blast import cli as func0
from chakin.commands.load.go import cli as func1
from chakin.commands.load.load_blast import cli as func2


@click.group()
def cli():
    pass


cli.add_command(func0)
cli.add_command(func1)
cli.add_command(func2)
