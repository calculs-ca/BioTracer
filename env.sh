#!/usr/bin/bash

# Adjust the following vars
export CONDA_MYHOME=/nfs3_ib/nfs-ip34/home/maxl/miniconda3
export BIO_TRACER_HOME=/nfs3_ib/nfs-ip34/home/maxl/BioTracer
export DRYPIPE_HOME=/nfs3_ib/nfs-ip34/home/maxl/BioTracer/DryPipe


conda activate BioTracerEnv

export DRYPIPE_PIPELINE_INSTANCE_DIR=$1
export PYTHONPATH=$BIO_TRACER_HOME:$DRYPIPE_HOME

echo "DRYPIPE_PIPELINE_INSTANCE_DIR=$DRYPIPE_PIPELINE_INSTANCE_DIR"

