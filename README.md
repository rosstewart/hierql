# hierql

## Repository layout

Main experiments: `experiments.py`

Transpiler: `hierql_transpiler.py`

LLM prompt: `llm_prompt.txt`
LLM raw table outputs: `llm_outputs.txt`
LLM output justifications: `llm_justifications.txt`

## Example usage:

To run on real data (you must download `gnomad.exomes.v4.1.1.sites.chrY.vcf.bgz`, `mondo.json`, and `gencc-submissions.tsv` from [https://gnomad.broadinstitute.org/downloads#v4], [https://mondo.monarchinitiative.org/pages/download/], and [https://thegencc.org/download], respectively):
```
python experiments.py --dataset real --repeats 5 --preview-rows 5
```

To run on synthetic toy data:
```
python experiments.py --dataset sample
```
