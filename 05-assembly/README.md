# 05. Genome assembly

Two assemblers were used on the filtered (non-human) reads from [step 4](../04-human-read-removal/).

## Option A: Flye (metagenome mode)

To assemble the long reads, Flye was run with the `--meta` flag, using enough threads to speed up the process. The input is the filtered file output from the human read removal step, which contains only microbial reads. The SLURM script is [`run_flye.sh`](run_flye.sh):

```sh
sbatch --nodes=1 --ntasks=1 --cpus-per-task=14 --mem=100G --time=12:00:00 --partition=interactive run_flye.sh
```

Visualize the `.gfa` file from the assembly using Bandage to easily identify circularized contigs.

## Option B: Myloasm

Myloasm is a recently developed assembler designed specifically for long-read sequencing data. In the authors' benchmarks of real Oxford Nanopore metagenomes, Myloasm recovered approximately three times as many complete circular contigs as the next-best assembler tested.

Shaw, J., Marin, M.G. & Li, H. High-resolution metagenome assembly for modern long reads with myloasm. _Nat Biotechnol_ (2026). https://doi.org/10.1038/s41587-026-03053-z

The `assembly_primary.fa` file is the main output contig file produced by Myloasm and is the file that should usually be used for downstream analysis. See [step 7](../07-closed-contigs/) for selecting circular contigs from either assembler.
