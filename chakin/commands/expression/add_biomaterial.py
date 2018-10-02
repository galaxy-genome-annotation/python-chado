import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, dict_output


@click.command('add_biomaterial')
@click.argument("organism_id", type=str)
@click.argument("biomaterial_name", type=str)
@click.option(
    "--organism_name",
    help="Orgnaism name (can be different from the organism ID provided)",
    type=str
)
@click.option(
    "--provider",
    help="Biomaterial provider name",
    type=str
)
@click.option(
    "--accession",
    help="Biomaterial accession",
    type=str
)
@click.option(
    "--sra_accession",
    help="SRA acession",
    type=str
)
@click.option(
    "--bioproject_accession",
    help="Bioproject accession",
    type=str
)
@click.option(
    "--attributes",
    help="Custom attributes (In JSON dict form)",
    type=str
)
@pass_context
@custom_exception
@dict_output
def cli(ctx, organism_id, biomaterial_name, organism_name="", provider="", accession="", sra_accession="", bioproject_accession="", attributes={}):
    """Add a new biomaterial to the database

Output:

    Job information
    """
    return ctx.gi.expression.add_biomaterial(organism_id, biomaterial_name, organism_name=organism_name, provider=provider, accession=accession, sra_accession=sra_accession, bioproject_accession=bioproject_accession, attributes=attributes)
