import click
from chakin.cli import pass_context
from chakin.decorators import custom_exception, str_output


@click.command('add_biomaterial')
@click.argument("biomaterial_name", type=str)
@click.argument("organism_id", type=str)
@click.option(
    "--description",
    help="Description of the biomaterial",
    type=str
)
@click.option(
    "--analysis_id",
    help="Analysis ID",
    type=str
)
@click.option(
    "--provider",
    help="Biomaterial provider name",
    type=str
)
@click.option(
    "--biosample_accession",
    help="Biosample accession",
    type=str
)
@click.option(
    "--sra_accession",
    help="SRA acession",
    type=str
)
@click.option(
    "--attributes",
    help="Custom attributes (In JSON dict form)",
    type=str
)
@pass_context
@custom_exception
@str_output
def cli(ctx, biomaterial_name, organism_id, description="", analysis_id="", provider="", biosample_accession="", sra_accession="", attributes={}):
    """Add a new biomaterial to the database

Output:

    Biomaterial id
    """
    return ctx.gi.expression.add_biomaterial(biomaterial_name, organism_id, description=description, analysis_id=analysis_id, provider=provider, biosample_accession=biosample_accession, sra_accession=sra_accession, attributes=attributes)
