"""Command-line interfaces for OpenASChem."""

import click

from ._version import __version__
from .art import logo
from .parameters import DEFAULT_RT_LEFT_LIMIT, DEFAULT_RT_RIGHT_LIMIT, validate_rt_limits

# Configure CLI context settings
CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"], max_content_width=120)


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo(f"OpenASChem Version {__version__}")
    click.echo(logo)
    ctx.exit()


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--version",
    is_flag=True,
    callback=print_version,
    expose_value=False,
    is_eager=True,
    help="Print version and exit.",
)
def openaschem():
    """OpenASChem: A data analysis software suite for Affinity Selection Mass Spectrometry."""
    pass


@openaschem.command(name="openidentify", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--file_dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory that contains only ASMS raw data files. Centroided .mzXML or .mzML files are supported.",
)
@click.option(
    "--chemical_pool_library",
    required=True,
    type=click.Path(exists=True),
    help="File name of the chemical pool library. .tsv, .csv and .xlsx formats are supported.",
)
@click.option(
    "--standard_prefix",
    default="STD100nM",
    required=True,
    type=str,
    help="File name prefix of the standard mass spectrometry runs.",
)
@click.option(
    "--r_home",
    default=None,
    hidden=True,
    help="Deprecated compatibility option; ignored by the Python alignment backend.",
)
@click.option(
    "--out", required=True, type=click.Path(exists=False), help="Directory for output files."
)
@click.option(
    "--mz_unit",
    default="ppm",
    show_default=True,
    type=click.Choice(["Da", "ppm"]),
    help="m/z unit for m/z tolerance settings.",
)
@click.option("--mz_tol", default="3", show_default=True, type=float, help="m/z tolerance.")
@click.option(
    "--rt_tol",
    default="0.45",
    show_default=True,
    type=float,
    help="Retention time tolerance in second.",
)
@click.option(
    "--acquisition_mode",
    default="positive",
    show_default=True,
    type=click.Choice(["positive", "negative"]),
    help="Mass spectrometry data acquisition mode.",
)
@click.option(
    "--quantification_method",
    default="intensity",
    show_default=True,
    type=click.Choice(["intensity", "area"]),
    help="Quantification method.",
)
@click.option(
    "--background_intensity",
    default="5000",
    show_default=True,
    type=float,
    help="Minimum intensity of a valid signal.",
)
@click.option(
    "--min_intensity_standard",
    default="100000",
    show_default=True,
    type=float,
    help="Minimum intensity of a valid signal in standard runs.",
)
@click.option(
    "--min_fold_change",
    default="5",
    show_default=True,
    type=float,
    help="Minimum fold change threshold for output.",
)
@click.option(
    "--min_fold_change_single",
    default="1",
    show_default=True,
    type=float,
    help="Minimum fold change threshold for pair-wise comparisons.",
)
@click.option(
    "--p_value_cutoff",
    default="0.05",
    show_default=True,
    type=float,
    help="Maximum p value threshold for output.",
)
@click.option(
    "--p_value_cutoff_single",
    default="0.01",
    show_default=True,
    type=float,
    help="Maximum p value threshold for pair-wise t-tests.",
)
@click.option("--n_threads", default=16, show_default=True, type=int, help="Number of threads.")
@click.option(
    "--peak_selection",
    type=click.Choice(["identity", "top1"]),
    default="identity",
    show_default=True,
    help="Standard-only identity scoring, or legacy highest-peak selection.",
)
@click.option(
    "--reference_relative_height",
    type=click.FloatRange(min=0.0, max=1.0, min_open=True),
    default=0.5,
    show_default=True,
    help="Minimum candidate height relative to the strongest accepted standard candidate.",
)
@click.option(
    "--reference_ambiguity_margin",
    type=click.FloatRange(min=0.0),
    default=0.15,
    show_default=True,
    help="Defer distinct competing references whose score difference is below this value.",
)
@click.option(
    "--identity_parent_library",
    type=click.Path(exists=True, dir_okay=False),
    default=None,
    help="Optional full CSV/TSV/XLSX compound library for parent-isotope checks. Defaults to input library.",
)
@click.option(
    "--rt_left_limit",
    default=DEFAULT_RT_LEFT_LIMIT,
    type=click.FloatRange(min=0),
    show_default=True,
    help="Lower RT limit in seconds; mean standard peak start must be above this value.",
)
@click.option(
    "--rt_right_limit",
    default=DEFAULT_RT_RIGHT_LIMIT,
    type=click.FloatRange(min=0),
    show_default=True,
    help="Upper RT limit in seconds; mean standard peak end must be below this value.",
)
def identify(
    file_dir,
    chemical_pool_library,
    standard_prefix,
    r_home,
    out,
    mz_unit,
    mz_tol,
    rt_tol,
    acquisition_mode,
    quantification_method,
    background_intensity,
    min_intensity_standard,
    min_fold_change,
    min_fold_change_single,
    p_value_cutoff,
    p_value_cutoff_single,
    n_threads,
    peak_selection,
    reference_relative_height,
    reference_ambiguity_margin,
    identity_parent_library,
    rt_left_limit,
    rt_right_limit,
):
    """OpenASChem binding compounds identification."""
    try:
        validate_rt_limits(rt_left_limit, rt_right_limit)
    except ValueError as error:
        raise click.BadParameter(
            str(error), param_hint="--rt_left_limit / --rt_right_limit"
        ) from error
    from .open_identify import open_identify

    open_identify(
        file_dir,
        chemical_pool_library,
        standard_prefix,
        out,
        mz_unit,
        mz_tol,
        rt_tol,
        acquisition_mode,
        quantification_method,
        background_intensity,
        min_intensity_standard,
        min_fold_change,
        min_fold_change_single,
        p_value_cutoff,
        p_value_cutoff_single,
        n_threads,
        peak_selection=peak_selection,
        reference_relative_height=reference_relative_height,
        reference_ambiguity_margin=reference_ambiguity_margin,
        identity_parent_library=identity_parent_library,
        rt_left_limit=rt_left_limit,
        rt_right_limit=rt_right_limit,
    )


@openaschem.command(name="openfilling", context_settings=CONTEXT_SETTINGS)
@click.option(
    "--file_dir",
    required=True,
    type=click.Path(exists=True),
    help="Directory that contains only ASMS raw data files. Centroided .mzXML or .mzML files are supported.",
)
@click.option(
    "--out",
    required=True,
    type=click.Path(exists=False),
    help="Directory for outputing file filling logs.",
)
@click.option(
    "--n_replicates",
    default="3",
    show_default=True,
    type=int,
    help="Number of replicates for each sample.",
)
def fill(file_dir, out, n_replicates):
    """OpenASChem file filling."""
    from .open_filling import open_filling

    open_filling(file_dir, out, n_replicates)


if __name__ == "__main__":
    openaschem()
