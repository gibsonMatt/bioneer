# Bioneer

A bioinformatics LLM companion tool. Ask questions in your terminal (e.g., "subset 5 vcfs to 'chr1' and merge"), and immediately get a `bcftools` command. 


Bioneer uses a training set of data curated from the bcftools vx.x. manual. Related examples of valid queries and repsonses are included with the user prompt. Related examples (stored locally in a Chroma database) are chosen based on semantic similarity using embeddings generated with `MiniLM-L6-v2`. 


Bioneer streamlines and optimizes queries by using dynamic few-shot prompt engineering to produce high quality results. The only output are "valid" bcftools commands.


## Usage
```
bioneer ask "filter vcf to remove all indels, subset to 'chr1', and then normalize all multiallelics to their own records"
```
```
bcftools view -i 'TYPE="snp"' -r chr1 -H your_vcf.vcf | bcftools norm -m- your_vcf.vcf
```

### Options
```
$ poetry run bioneer ask --help
Usage: bioneer ask [OPTIONS] QUERY

  Main function for querying the model using dynamic few-shot injection.

  Parameters:     
  
    query (str): The query to ask the LLM     
  
    degree (int): Number of similar prompts to include   
  
    force (bool): Force re-create Chroma db

Options:
  -d, --degree INTEGER  Number of similar prompts to include in dynamic few-
                        shot injection
  -f, --force BOOLEAN   Force re-create Chroma db
  --help                Show this message and exit.

```


## Development

```
poetry install --with dev
```

```
poetry run bioneer ask [query]
```

or

```
poetry shell
bioneer ask [query]
```

## Update examples

Append examples to this file [./data/data_bcftools_view.txt](file).

