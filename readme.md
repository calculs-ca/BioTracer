
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


Copy and customize the file BioTracer/env.sh

```bash 

cp <Location of BioTracer repo clone>/env.sh <Location of BioTracer repo clone>/myenv.sh

# edit <Location of BioTracer repo clone>/myenv.sh

```

# Usage

### Initialize shell:

To run a given pipeline, create a dir for the pipeline instance <PIPELINE_INSTANCE_DIR>. **You need to provide full path 
to instance directory including NFS mount path for access within compute nodes.** Set the environment with: 

```bash 
source <Location of BioTracer repo clone>/myenv.sh <PIPELINE_INSTANCE_DIR_FULLPATH>
```

# Run with test data to test installation:

```bash 

cd ${BIO_TRACER_HOME}
unzip bio-tracer-example.zip
cp example-config.yml <PIPELINE_INSTANCE_DIR>/config.yaml
# edit the config.yaml accordingly

# run to test install


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


