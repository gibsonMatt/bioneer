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


@dataclass
class VectorStoreHandle:
    auth: AuthHandle()
    force: bool
    type: str = "default"
    degree: int = 3
    initialized: bool = True

    def __post_init__(self):
        logger = logging.getLogger(__name__)

        load_dotenv()
        path = os.getenv("PROMPT_EXAMPLES_PATH")
        logger.debug(f"Using prompt examples from {path}")

        persistent = os.getenv("VECTORSTORE")
        logger.debug(f"Using vectorstore {persistent}")

        embedding_function = SentenceTransformerEmbeddings(
            model_name="all-MiniLM-L6-v2"
        )
        logger.debug(f"Using embedding model {embedding_function.model_name}")

        if os.listdir(persistent) == []:
            self.initialized = False
        else:
            self.initialized = True

        if self.force:
            self.initialized = False

        if self.initialized == False:
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

        self.example_selector = SemanticSimilarityExampleSelector(
            vectorstore=vectorstore, k=self.degree
        )

    def get_examples(self, query: Query):
        return self.example_selector.select_examples({"query": query.query})
