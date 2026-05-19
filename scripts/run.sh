#!/bin/bash

JOB_NAME=$1

sbatch -p gpu -n 1 -J "${JOB_NAME}" --wrap="./Rectangular2DIsing"
