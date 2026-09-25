
Describe steps:

`main.py` runs the whole thing
1. 
2. gives you fasta and db format for RNA-Composer
2b. you need to get structures from RNA-Composer yourself
3.
4. 

```
python main.py `
  --model-command "your-modeler {fasta} {dotbracket} {output}" `
  --relax-command "your-oxrna-command {input} {output} {steps} {seed}" `
  --steps 1000
```