from dataclasses import dataclass
from uuid import uuid4, UUID
from bioneer.auth import AuthHandle
from langchain.prompts import (
    ChatPromptTemplate,
    FewShotChatMessagePromptTemplate,
)
from langchain.chat_models import ChatOpenAI


@dataclass
class Query():
    query: str
    auth: AuthHandle

    id: UUID = uuid4()

    def __post_init__(self):
        assert len(self.query) < 300 and len(self.query) > 5, "Query either too small or too large"

    def run(self, examples):
        # assemble prompt

        formatting = """Return only the command in Markdown format, e.g.:

        ```bash
        ...
        ```

        Do not include any additional explanation.
        """

        example_prompt = ChatPromptTemplate.from_messages(
            [
                ("human", "{query}"),
                ("ai", "{response}"),
            ]
        )

        few_shot_prompt = FewShotChatMessagePromptTemplate(
            example_prompt=example_prompt,
            examples=examples,
        )

        final_prompt = ChatPromptTemplate.from_messages(
            [
                ("system", "You are a helpful AI bot that writes bcftools command line commands to solve user's requestd operation on their vcf or bcf file"),
                ("system", formatting),
                few_shot_prompt,
                ("human", "{query}")
            ]
        )

        chain = final_prompt | ChatOpenAI(temperature=0.0, openai_api_key = self.auth.api_key)

        self.response = chain.invoke({"query": self.query})

        

    def validate_response():
        pass

    def format_response():
        pass

    


    
