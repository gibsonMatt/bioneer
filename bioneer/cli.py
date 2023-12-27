import click
from dataclasses import dataclass

from bioneer.auth import AuthHandle
from bioneer.query import Query
from bioneer.vectorstore import VectorStoreHandle

from bioneer.logging import setup_logger


@click.group()
def cli():
    logger = setup_logger("bioneer")


@cli.command()
@click.argument("query")
@click.option(
    "-d",
    "--degree",
    "degree",
    default=3,
    type=int,
    help="Number of similar prompts to include in dynamic few-shot injection",
)
@click.option(
    "-f", "--force", "force", default=False, type=bool, help="Force re-create Chroma db"
)
def ask(query: str, degree: int, force: bool):
    """
    Main function for querying the model using dynamic few-shot injection.

    Parameters:
        query (str): The query to ask the LLM
        degree (int): Number of similar prompts to include in dynamic few-shot injection
        force (bool): Force re-create Chroma db

    Returns:
    """

    # verify that api key is provided
    auth = AuthHandle()
    auth.configure()

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

    # Pretty-print response
    click.echo(click.style(my_query.response.content, fg="green", bold=True))
