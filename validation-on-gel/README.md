
Describe steps:

`main.py` runs the whole thing
1. 
2. gives you fasta and db format
2b. you need to get structures from farfar2
    ```
    python main.py `
  --farfar2 `
  --farfar2-executable "path/to/rna_denovo.default.linuxgccrelease" `
  --farfar2-nstruct 10
    ```
3.
4. 

```
python main.py `
  --model-command "your-modeler {fasta} {dotbracket} {output}" `
  --farfar2 `
  --farfar2-executable "path/to/rna_denovo.default.linuxgccrelease" `
  --farfar2-nstruct 10 `
  --relax-command "your-oxrna-command {input} {output} {steps} {seed}" `
  --steps 1000
```



FARFAR2 uses the step-2 FASTA and dot-bracket files as its sequence and
secondary-structure inputs. It writes Rosetta silent output and dumped PDB
models under `farfar2/<construct>/`. Install Rosetta separately and provide
the path to its `rna_denovo` executable with `--farfar2-executable`.