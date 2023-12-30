# Bioneer

A bioinformatics LLM companion tool. Ask questions in your terminal (e.g., "subset 5 vcfs to 'chr1' and merge"), and immediately get a `bcftools` command. 


bioneer uses a training set of data curated from the bcftools vx.x. manual. Related examples of valid queries and repsonses are included with the user prompt. Related examples (stored locally in a Chroma database) are chosen based on semantic similarity using embeddings generated with `MiniLM-L6-v2`. 


Bioneer streamlines and optimizes queries by using dynamic few-shot prompt engineering to produce high quality results. The only output are "valid" bcftools commands.


## Usage
```
bioneer ask "remove all indels from a bcf file, subset to contig '1'"
```

### Result
```
bcftools view -V indels -r 1 your_file.bcf
```

### Advanced usage

More complex operations are also possible

```
bioneer ask "exclude all indels, subset to contig named 'contig1', remove SNPs if more than 20 samples have missing genotypes, remove all multi-allelic sites"
```

#### Result
```
bcftools view -V indels,mnps,other -r contig1 -i 'F_MISSING <= 0.2' -m2 -M2 your_vcf.vcf
```

## Configure

`bioneer` will automatically configure itself, with the exception of your `OPENAI_APIKEY`. Set this variable yourself, for example by running

```
export OPENAI_APIKEY="you-api-key"
```

Vectorstores and example datasets will be stored at `~/.bioneer/`. To change where `bioneer` looks for example data or where the Chroma database is stored, you can optionally set the following environmental variables:

```
VECTORSTORE="path"
PROMPT_EXAMPLES_PATH="path"
```

## Options
```
$ poetry run bioneer ask --help
Usage: bioneer ask [OPTIONS] QUERY

  ask bioneer

  Parameters:

      query (str): The query to ask the LLM

      degree (int): Number of similar prompts to include in dynamic few-shot
      injection

      force (bool): Force re-create Chroma db

  Returns:

Options:
  -d, --degree INTEGER  Number of similar prompts to include
  -f, --force           Force re-create Chroma db
  --help                Show this message and exit.
```


## Installation and Development

`bioneer` uses Poetry

### Install dependencies with dev requirements

```
poetry install --with dev # installs with isort, pytest, and black
```

### Execute in Poetry environment
```
poetry run bioneer ask [query]
```

### Hooks

pytest, black, and isort are executed as pre-commit hooks


### Tests

Tests as defined in [tests/](file)

## Update examples

Append examples to [data/data_bcftools_view.txt](file)

