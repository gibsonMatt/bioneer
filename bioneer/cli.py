import click
from dataclasses import dataclass

from bioneer.auth import AuthHandle
from bioneer.query import Query
from bioneer.vectorstore import VectorStoreHandle

import logging
import sys


@click.group()
@click.option(
    "-v", "--verbose", "verbose", default=False, type=bool, help="Verbose", is_flag=True
)
def cli(verbose: bool):
    if verbose:
        logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
    else:
        logging.basicConfig(stream=sys.stdout, level=logging.ERROR)


@cli.command()
@click.argument("query")
@click.option(
    "-d",
    "--degree",
    "degree",
    default=3,
    type=int,
    help="Number of similar prompts to include",
)
@click.option(
    "-f",
    "--force",
    "force",
    default=False,
    is_flag=True,
    type=bool,
    help="Force re-create Chroma db",
)
def ask(query: str, degree: int, force: bool):
    """
    ask bioneer\n

    Parameters:\n
        query (str): The query to ask the LLM\n
        degree (int): Number of similar prompts to include in dynamic few-shot injection\n
        force (bool): Force re-create Chroma db\n

    Returns:
    """

    logger = logging.getLogger(__name__)

    # verify that api key is provided
    auth = AuthHandle()
    auth.configure()

    # pick model
    available_models = auth.get_available_models()
    if "gpt-4-1106-preview" in available_models:
        auth.model = "gpt-4-1106-preview"
    else:
        auth.model = "gpt-4"
    logger.debug(f"Using model {auth.model}")

    # configure the vector store for getting related prompt examples.
    # will create new local Chroma db at path `VECTORSTORE` using prompts
    # at `PROMPT_EXAMPLES_PATH` or if data base exists will load it in directly.
    # Creates db with embedding model "all-MiniLM-L6-v2".
    my_vector_store = VectorStoreHandle(auth, force, "default", degree)

    # Make query object
    my_query = Query(query=query, auth=auth)

    # Get related examples from the vector store Chroma db
    prompt_response_examples = my_vector_store.get_examples(my_query)

    # Run the bioneer query
    my_query.run(prompt_response_examples)

    # Validate
    # my_query.validate_response()

    # Pretty-print response
    click.echo(click.style(my_query.response.content, fg="green", bold=True))
