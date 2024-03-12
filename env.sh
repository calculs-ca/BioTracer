

# Adjust the following vars
CONDA_HOME=/nfs3_ib/nfs-ip34/home/maxl/miniconda3
BIO_TRACER_HOME=/nfs3_ib/nfs-ip34/home/maxl/BioTracer
DRYPIPE_HOME=/nfs3_ib/nfs-ip34/home/maxl/BioTracer/DryPipe


conda activate BioTracerEnv

export DRYPIPE_PIPELINE_INSTANCE_DIR=$1
export PYTHONPATH=$BIO_TRACER_HOME:$DRYPIPE_HOME

echo "DRYPIPE_PIPELINE_INSTANCE_DIR=$DRYPIPE_PIPELINE_INSTANCE_DIR"

