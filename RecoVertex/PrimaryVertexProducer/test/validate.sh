#!/bin/bash


eval `scramv1 runtime -sh`

cmsRun vertexTest.py
cmsRun vertexTest.py gpu=False

harvestTrackValidationPlots.py test_dqm_gpu.root -o gpu.root
harvestTrackValidationPlots.py test_dqm_cpu.root -o cpu.root

rm -r plots/

makeTrackValidationPlots.py gpu.root cpu.root --png 

cp -r plots/ /nfshome0/gpizzati/

