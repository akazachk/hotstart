#!/usr/bin/env bash

PROJ_DIR=".."
EXECUTABLE=${PROJ_DIR}/Debug/main

# Command below _works_
${EXECUTABLE} ./neos-5078479-escaut_presolved.mps.gz 0

# Command below gives memory access error
#${EXECUTABLE} ./neos-5078479-escaut_presolved.mps.gz -1

# Command below gives incorrect determination of infeasibility
#${EXECUTABLE} ./neos-5078479-escaut_presolved.mps.gz 2
