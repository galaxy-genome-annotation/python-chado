import click
from chakin.cli import pass_context, json_loads
from chakin.decorators import custom_exception, dict_output


@click.command('add_biomaterial')
@click.argument("biomaterial_name", type=str)
@click.argument("organism_id", type=int)
@click.option(
    "--description",
    help="Description of the biomaterial",
    type=str
)
@click.option(
    "--biomaterial_provider",
    help="Biomaterial provider name",
    type=str
)
@click.option(
    "--biosample_accession",
    help="Biosample accession number",
    type=str
)
@click.option(
    "--sra_accession",
    help="SRA accession number",
    type=str
)
@click.option(
    "--bioproject_accession",
    help="Bioproject accession number",
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
def cli(ctx, biomaterial_name, organism_id, description="", biomaterial_provider="", biosample_accession="", sra_accession="", bioproject_accession="", attributes={}):
    """Add a new biomaterial to the database

Output:

    Biomaterial details
    """
    return ctx.gi.expression.add_biomaterial(biomaterial_name, organism_id, description=description, biomaterial_provider=biomaterial_provider, biosample_accession=biosample_accession, sra_accession=sra_accession, bioproject_accession=bioproject_accession, attributes=attributes)
