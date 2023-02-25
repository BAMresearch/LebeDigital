# Snakemake
Implementation of LEBEDIGITAL snakemake with various real scripts: ver#1 including units extensivley


IMPORTANT!!!
Constant parameters for each subprocess are defined in the "Inputs/Constant_parameters.json" file

The "Lebedigital_SUPERDICTIONARY.json" stores EVERY input parameter and result for each sub-process

The temporal files are individual files containing the output of each rule, their content is also stored in the SUPERDICTIONARY. Their main purposes for snakemake to form the DAG and run it, by checking which rules' output individual file is fed to other rule as input. If requiered, they can be saved as well.

each rule loads its own environment based on the /Envs folder files. Is not mandatory to do so in every rule.

because the homogeneization rule inputs parameters that come from other rules still not implemented, I just wrote them manually for now in the SUPERDICTIONARY, soon I will adapt these two rules to the big mock algorithm done prevoiusly.

the FEM_model rule also inputs parameters from other subprocesses; for now I wrote them in the constant parameters inputs.

for the moment, im just running the program with a forced run: "snakemake --cores 1  -F --use-conda" ,other types of running are not tested yet

Due to the fact that some of the scripts are currently in developement, some values have been stored directly (for the moment) into the SUPERDICTIONARY. It is not recommended to erase its content before each run. The same applies to some files in the "Results" folder.

Finally, detailed examples of Snakemake were done on passed versions of the algorithm scripts. They can be accessed as passed commits in the GitHub 52 branch. 

## Running
Open a terminal in the library where the snakefile is contained (in this case LebeDigital/usecases/demonstrator/workflow) and activate snakemake:

```sh
conda activate snakemake
```

The workflow can be run by targeting the Final_Output_Files rule with
```sh
snakemake --cores 1 Final_Output_Files --use-conda
```
or by doing a forced run
```sh
snakemake --cores 1 -F --use-conda
```
if there is interest just to target one specific file:
```sh
snakemake --cores 1 Results/Out_Demolding_Time.txt --use-conda
```

The workflow DAG plot can be obtained by running the snakemake command:
```sh
snakemake --forceall --dag | dot -Tpdf > dag.pdf
```
