# Contributor Guide

Thank you for your interest in improving this project!

We welcome contributions from anyone from Owkin in the form of suggestions, bug reports, pull requests, and feedback. This document will guide you through the process of contributing to this project.

Here is a list of important resources for contributing:

- [Source Code]
- [Issue Tracker]
- [Documentation]

## How to report a bug

Report bugs on the [Issue Tracker].

When filing an issue, make sure to answer these questions:

- Which Python version are you using?
- Which version of this project are you using?
- What did you do?
- What did you expect to see?
- What did you see instead?

The best way to get your bug fixed is to provide a test case,
and/or steps to reproduce the issue.

## How to request a feature

Request features on the [Issue Tracker].

## How to set up your development environment

You need [uv]. (Yes that's it, you don't need to install anything else)

### Set up the environment

Retrieve the project from the repository and navigate to the root directory of the project.

```console
git clone git@github.com:owkin/PyDeconv.git
cd PyDeconv
```

To install all required dependencies and the pinned version of Python, run:

```console
uv sync
```

It will install the python dependencies needed for the project.

## Quality checks

If you have installed the pre-commit hooks, they will be run automatically before each commit.

You can also run them manually with:

```console
# All checks used in the pre-commit hooks
make checks

# Format the code
make fmt

# Run the linter
make lint
```

## Test suite

Unit tests are located in the _tests_ directory,
and are written using the [pytest] testing framework.

To run the test suite, run:

```console
make tests
```

## Documentation

The documentation is written using [mkdocs-material].

To serve the documentation locally, run:

```console
make docs-serve
```

To build the documentation, run:

```console
make docs-build
```

## Start working on a new feature

**It is recommended to open an issue before starting work on anything.**
This will allow a chance to talk it over with the owners and validate your approach.

> We follow the [GitHub flow] branching model.

To start working on a new feature, create a new branch from the `main` branch:

```console
$ git checkout main
$ git pull
$ git checkout -b my-new-feature
```

You can also create a new branch directly from the GitHub interface.

## How to submit changes

Once you are done working on your feature, push your branch to the remote repository and open a [Pull Request].

We will review your changes and merge them into the `main` branch if everything looks good.

Your pull request should include:
- a clear description of the changes you made
- an update of the documentation (if needed)
- a link to the issue you are fixing (if any)
- tests for the new code (if applicable)

Feel free to submit early drafts of your pull request if you want to get feedback on your work in progress. We can always help iterate on your changes.

[conda]: https://docs.conda.io/en/latest/
[documentation]: https://owkin.github.io/PyDeconv/
[github flow]: https://docs.github.com/en/get-started/using-github/github-flow
[issue tracker]: https://github.com/owkin/PyDeconv/issues
[uv]: https://docs.astral.sh/uv/getting-started/installation/
[pull request]: https://github.com/owkin/PyDeconv/pulls
[pytest]: https://pytest.readthedocs.io/
[source code]: https://github.com/owkin/PyDeconv
