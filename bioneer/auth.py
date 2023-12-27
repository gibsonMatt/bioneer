# auth.py

from dataclasses import dataclass
from dotenv import load_dotenv
import os
from openai import OpenAI
from openai import AuthenticationError
from bioneer.logging import setup_logger
from logging import Logger


@dataclass
class AuthHandle:
    # TO DO: handle selecting model. Default to GPT4, but select 3.5 if not available.

    logger: Logger = setup_logger("auth")

    def configure(self):
        self.logger.debug("Loading env variables")
        load_dotenv()
        self.is_configured()

    def is_configured(self):
        key = os.getenv("OPENAI_API_KEY")
        if key != None:
            self.api_key = key
            self.logger.debug("Authenticated")
        else:
            self.logger.error("Invalid or missing api key")
