#!/bin/bash

for i in *.job
do
        sbatch $i
done
