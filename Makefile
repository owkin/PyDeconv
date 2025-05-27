.PHONY: checks, docs-serve, docs-build, fmt, config, type, tests

checks:
	uvx pre-commit run --all-files

docs-serve:
	uv run mkdocs serve

docs-build:
	uv run mkdocs build

fmt:
	uv run ruff format pydeconv tests

lint:
	uv run ruff check --fix pydeconv tests

type:
	uv run mypy pydeconv --install-types --non-interactive --show-traceback

tests:
	uv run pytest --cov=pydeconv --cov-report=term-missing tests/ -s -vv
