# 02. Super high accuracy (sup) basecalling using Dorado

You can do this either locally using MinKNOW or using Dorado on the ARC. The latter is faster.

## Installing Dorado

Dorado can be downloaded [here](https://github.com/nanoporetech/dorado) onto a local device and copied via `scp` to its desired directory. For more information on installing Dorado see the Bacterial-Methylation-Workflow (https://github.com/jonbakerlab/Baker-Lab-Scripts/commit/ee94e5f2da2773b992c020b29905b9baf3f4c267).

Alongside `dorado`, the sup model also needs to be downloaded:

```sh
./dorado download --model dna_r10.4.1_e8.2_400bps_sup@v5.0.0
```

## Running the sup model with Dorado on the ARC

With all of the `.pod5` files, Dorado can be used to basecall with the sup model by providing the appropriate directory. Adapter sequences are trimmed by default once demultiplexing is completed. The script [`run_dorado.sh`](run_dorado.sh) outputs quality-filtered (>Q9), demultiplexed `.fastq` files. Set `working_dir` before submitting:

```sh
sbatch run_dorado.sh
```

## Notes

If using MinKNOW, make sure to select the appropriate basecalling model and check the options for quality filtering (e.g. >Q9), demultiplexing and barcode trimming. The basecalling model can be found in the first line of the `.fastq` files (check with `head`), in `basecall_model_version_id`, where `sup` stands for 'super high accuracy'. Other acronyms include `fast` (fast basecalling model) and `hac` (high accuracy).
