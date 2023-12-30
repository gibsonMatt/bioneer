import click
from click.testing import CliRunner

from bioneer.cli import ask


def test_ask_command_no_input():
    runner = CliRunner()
    result = runner.invoke(ask, input="")
    assert result.exit_code == 2
    # assert result.exception is AssertionError


def test_ask_command_long_input():
    runner = CliRunner()
    result = runner.invoke(ask, input="a" * 301)
    assert result.exit_code == 2
