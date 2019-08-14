import click
from chakin.commands.load.blast import cli as blast
from chakin.commands.load.go import cli as go
from chakin.commands.load.interpro import cli as interpro


@click.group()
def cli():
    pass


cli.add_command(blast)
cli.add_command(go)
cli.add_command(interpro)
