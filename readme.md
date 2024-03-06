
# Setup

### 1. Install miniconda 


### 2. checkout BioTracer

```bash 
git clone git@github.com:calculs-ca/BioTracer.git
```

### 3. checkout DryPipe

```bash 
git clone git@github.com:calculs-ca/DryPipe.git
```

### 4. create conda env

```bash 
conda env create -f <Location of BioTracer repo clone>/conda-env.yml
```

The BioTracerEnv conda environment is now created


### 5. customize env.sh file


customize the file BioTracer/env.sh

```bash 
#!/usr/bin/bash

BIO_TRACER_HOME=...set appropriately...
DRYPIPE_HOME=...set appropriately...


conda activate BioTracerEnv

export DRYPIPE_PIPELINE_INSTANCE_DIR=$1
export PYTHONPATH=$BIO_TRACER_HOME:$DRYPIPE_HOME

echo "DRYPIPE_PIPELINE_INSTANCE_DIR=$DRYPIPE_PIPELINE_INSTANCE_DIR"
```

# Usage

### Initialize shell:

To run a given pipeline, create a dir for the pipeline instance <PIPELINE_INSTANCE_DIR>, then set the environment with: 

```bash 
source BioTracer/env.sh <PIPELINE_INSTANCE_DIR>
```

### Run

```bash 

python3 -m dry_pipe.cli run --generator=bio_tracer:pipeline
```

### Prepare

will create tasks, without executing

```bash 

python3 -m dry_pipe.cli prepare --generator=bio_tracer:pipeline

```

### Execute single task, in foreground

```bash 

python3 -m dry_pipe.cli task $PIPELINE_INSTANCE_DIR/.drypipe/fastp.RK1723_20231106_R00 --wait
```
