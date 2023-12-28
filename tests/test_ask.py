import click
from click.testing import CliRunner
from bioneer.cli import ask


def test_ask_command_no_input():
    runner = CliRunner()
    result = runner.invoke(ask, input="")
    assert result.exception is AssertionError
