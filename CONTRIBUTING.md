# Contributing

Install with `python -m pip install -e ".[gui,dev]"`.
Use a feature branch and keep changes focused. Before submitting a pull request,
run `ruff check .`, `ruff format --check .` and `python -m pytest`.
CI also builds the source distribution and wheel on Windows and Linux.

Use four-space indentation, double-quoted strings, snake_case Python functions
and a 100-character formatting target. Ruff configuration is in `pyproject.toml`.
Keep numerical behavior changes separate from cleanup and include meaningful
regression cases for extraction, alignment, quantification or reference changes.
CLI option names and output schemas are compatibility interfaces.

Bug reports should include application/Python versions, command or GUI settings,
expected behavior and relevant logs. Prefer small synthetic or shareable examples.
Do not commit raw experimental data, generated reports, private local paths,
caches or credentials.

Before a release, update the changelog, run tests and the smoke example, check
installed CLI/GUI entry points and report templates, then `python -m build`.
Validate result identities, reference IDs and numerical quantities against the
previous release when changing the analysis code. Release archives and validation
manifests are generated artifacts, not source files maintained by hand.
