# Bioneer

Your bioinformatics LLM companion tool

As a bioinformatician I constantly find myself looking up the same common commands every day. Common searches include:

> how to filter a VCF by minor allele frequency?   
> how to remove multiallelics from a VCF?   
> how to merge VCF files   
etc

Bioneer makes it easy to ask these questions in your terminal, and immediately get a `bcftools` command. 

Bioneer streamlines and optimizes queries to an LLM by using dynamic few-shot prompt engineering to produce high quality results every time. The only output is a single bcftools one-liner to accomplish your task. 



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

  Parameters:     query (str): The query to ask the LLM     degree (int):
  Number of similar prompts to include in dynamic few-shot injection     force
  (bool): Force re-create Chroma db

  Returns:

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