
# Setup

### 1. Install miniconda 


### 2. checkout BioTracer

```bash 
git clone git@github.com:calculs-ca/BioTracer.git
```

### 3. checkout DryPipe

```bash 
git clone git@github.com:calculs-ca/DryPipe.git

cd DryPipe
git fetch origin wip-dsl-overhaul
git checkout wip-dsl-overhaul
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

### Run the pipeline

```bash 

python3 -m dry_pipe.cli run --generator=bio_tracer:pipeline
```

When tasks fail, they can be run individualy. 

# Debugging

### Prepare Tasks

will create tasks, without executing

```bash 

python3 -m dry_pipe.cli prepare --generator=bio_tracer:pipeline

```

### Execute single task, --wait --tail

+ option "--wait" launches the task in the foreground 
+ option "--tail" tails the task logs:  out.log (stdout and stderr), and pipeline.log
+ option "--tail-all" tails the task logs:  out.log (stdout and stderr), and pipeline.log to STDOUT


### Relaunch failed tasks in array

```bash 


python3 -m dry_pipe.cli restart-failed-array-tasks --task-key=array_parent
```

# Run Example:

```bash 

mkdir <the_example_instance_dir>

cp BioTracer/example-config.yml <the_example_instance_dir>          

```
