from dataclasses import dataclass
from langchain.embeddings import OpenAIEmbeddings
from langchain.prompts import SemanticSimilarityExampleSelector
from langchain.vectorstores import Chroma
from bioneer.auth import AuthHandle
from bioneer.query import Query
from dotenv import load_dotenv
from langchain.document_loaders import TextLoader
from langchain.embeddings.sentence_transformer import SentenceTransformerEmbeddings
import logging
import os
import json
import os
import requests


@dataclass
class VectorStoreHandle:
    auth: AuthHandle()
    force: bool
    type: str = "default"
    degree: int = 3
    initialized: bool = True
    examples_url: str = "https://raw.githubusercontent.com/gibsonMatt/bioneer/main/data/data_bcftools_view.txt"

    def __post_init__(self):
        logger = logging.getLogger(__name__)

        # re-initializing will update the examples and recreate the database

        if self.force:
            self.initialized = False

        load_dotenv()
        if self.initialized == False:
            path = os.getenv("PROMPT_EXAMPLES_PATH")
            if path == None:
                file_path = os.path.expanduser("~/.bioneer/bioneer_data.txt")
                directory = os.path.dirname(file_path)

                if not os.path.exists(directory):
                    os.makedirs(directory)

                if not os.path.exists(file_path):
                    url = self.examples_url
                    response = requests.get(url)
                    with open(file_path, "wb") as file:
                        file.write(response.content)
            else:
                logger.debug(f"Using prompt examples from {path}")

            path = file_path

            persistent = os.getenv("VECTORSTORE")

            if persistent == None:
                persistent = os.path.expanduser("~/.bioneer/vectorstore")
                logger.debug(f"Using default vectorstore {persistent}")
            else:
                logger.debug(f"Using vectorstore {persistent}")

            embedding_function = SentenceTransformerEmbeddings(
                model_name="all-MiniLM-L6-v2"
            )
            logger.debug(f"Using embedding model {embedding_function.model_name}")

            if self.type == "default":
                with open(path, "r") as file:
                    data = json.load(file)
                to_vectorize = [" ".join(d.values()) for d in data]
                vectorstore = Chroma.from_texts(
                    to_vectorize,
                    embedding_function,
                    metadatas=data,
                    persist_directory=persistent,
                )
                self.initialized = True

        elif self.initialized == True:
            vectorstore = Chroma(
                persist_directory=persistent, embedding_function=embedding_function
            )

        # define the example_selector function, set as attribute
        self.example_selector = SemanticSimilarityExampleSelector(
            vectorstore=vectorstore, k=self.degree
        )

    def get_examples(self, query: Query):
        return self.example_selector.select_examples({"query": query.query})
