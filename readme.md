# RNA-seq workshop @ ASHG 2019

This is the main readme for the the RNA-seq hands-on workshop @ ASHG 2019

## First try with Snakemake and aberrant-expression-analysis pipeline
```
# Dry run for the pipeline with snakemake/wbuild
snakemake \
    --configfile src/wbuild/expression_analysis.snake \
    --snakefile ../aberrant-expression-pipeline/Snakefile \
    -n
```

