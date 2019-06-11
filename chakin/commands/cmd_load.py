import click
from chakin.commands.load.blast import cli as func0
from chakin.commands.load.go import cli as func1


@click.group()
def cli():
    pass


cli.add_command(func0)
cli.add_command(func1)
