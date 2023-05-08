# Pipeline setup

For a detailed description of the steps carried out during read processing, genotyping and annotation, see [README.md](README.md).

| Article                | PALEOMIX configuration                          | PALEOMIX version                                                                                           | Conda environment                     | AnnoVEP version                                                                                           |
| ---------------------- | ----------------------------------------------- | ---------------------------------------------------------------------------------------------------------- | ------------------------------------- | --------------------------------------------------------------------------------------------------------- |
| Gul *et al*. 2023      | [genotyping.yaml](gul2023/genotyping.yaml)      | [4bbc759](https://github.com/Hansen-Group/NGSPipeline/releases/download/snapshots/paleomix-4bbc759.tar.gz) | [conda.yaml](gul2023/conda.yaml)      | [f418dcf](https://github.com/Hansen-Group/NGSPipeline/releases/download/snapshots/annovep-f418dcf.tar.gz) |
| Johansen *et al.* 2023 | [genotyping.yaml](johansen2023/genotyping.yaml) | [394a2bf](https://github.com/Hansen-Group/NGSPipeline/releases/download/snapshots/paleomix-394a2bf.tar.gz) | [conda.yaml](johansen2023/conda.yaml) | [f418dcf](https://github.com/Hansen-Group/NGSPipeline/releases/download/snapshots/annovep-f418dcf.tar.gz) |
| Peña *et al*. 2023     | [genotyping.yaml](pena2023/genotyping.yaml)     | [304d4ef](https://github.com/Hansen-Group/NGSPipeline/releases/download/snapshots/paleomix-304d4ef.tar.gz) | [conda.yaml](pena2023/conda.yaml)     | [f418dcf](https://github.com/Hansen-Group/NGSPipeline/releases/download/snapshots/annovep-f418dcf.tar.gz) |

## The genotyping pipeline

The genotyping pipeline is based on [PALEOMIX](https://github.com/MikkelSchubert/paleomix). The table above contains snapshots of the versions used for a given study, as well as conda environment corresponding to the environment in which the pipelines were run.

### Installing the genotyping pipeline

1. Install [Conda](https://docs.conda.io/projects/conda/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

2. Download the `conda.yaml` file, the `genotyping.yaml` file, and the PALEOMIX snapshot corresponding to the study of interest. All are listed in the table above.

3. Create a Conda environment using the downloaded `conda.yaml` file:

    ```bash
    conda env create --name paleomix --file conda.yaml
    ```

4. Install the PALEOMIX snapshot in the Conda environment:

    ```bash
    tar xvzf paleomix-*.tar.gz
    cd ./paleomix/
    conda run --name paleomix python3 setup.py install
    ```

5. Download the [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) and place the files a folder named `resources` in the same folder as the pipeline configuration file (`genotyping.yaml`).

### Running the genotyping pipeline

The pipeline can be run either using the `conda run` command as shown above:

```bash
conda run --name paleomix paleomix ngs run genotyping.yaml
```

Or by first activating the conda environment:

```bash
conda activate paleomix
paleomix ngs run genotyping.yaml
```

### Genotyping results

By defaults, resulting genotypes will be placed in `${filename}.output` where `${filename}` corresponds to the configuration filename without the `.yaml` extension. E.g. for the configuration file included in this repository, output will be placed in `genotyping.output`:

* The final genotypes (VCF) for all samples can be found in `genotyping.output/genotypes/`.
* Per sample alignment files (BAMs) can be found `genotyping.output/alignments/`.
* Per sample haplotypes (g.VCFs) can be found in `genotyping.output/haplotypes/`.
* Statistics and reports for all stages of the pipeline are located in `genotyping.output/statistics/`.
* Temporary files can be found in `genotyping.output/cache/`. This folder can safely be deleted once the pipeline has been run to completion.

## The annotation pipeline

The annotation pipeline is based on [AnnoVEP](https://github.com/cbmrphenomics/annovep). The annotation pipeline is intended to be run using podman, but may used with docker docker as well.

### Installing the annotation pipeline

1. Download the AnnoVEP snapshot corresponding to the study of interest (see the table above).

2. Unpack and build image using `podman`:

    ```bash
    tar xvzf annovep-*.tar.gz
    cd ./annovep/
    make
    ```

    If using `docker`, run `make` as follows:

    ```bash
    make MANAGER=docker
    ```

3. Download custom annotations used by the pipeline:

    ```bash
    ./bin/annovep setup
    ```

    If using `docker`, instead run

    ```bash
    export ANNOVEP_RUNNER=docker
    ./bin/annovep setup
    ```

    By default this will place the annotation files in `~/annovep`, but this behavior may be changed by setting the `ANNOVEP_CACHE` environmental variable. Note that the `setup` step downloads approximately 350 GB of data and requires approximately 150 GB of space when done.

### Running the annotation pipeline

The pipeline takes as input a single VCF file and generated several output files using a user-supplied output prefix:

```bash
./bin/annovep pipeline input.vcf.gz output
ls output.tsv*
output.tsv output.tsv.columns
```

If no output prefix is supplied, the resulting files will be named using the base name of the input file (e.g. `input` for `input.vcf.gz`).
