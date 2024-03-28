#!/bin/bash


# Inputs
GDNAM=SC_2019
BRAVEShome=/media/leohoinaski/HDD/BRAVES_database
mcipPath=/media/leohoinaski/HDD/SC_2019
YEAR=2019

# ---------------------Running BRAVES_database-----------------------
# Check shRunnerBRAVES_database.py for more input configurations
./shRunnerBRAVES_database.py ${BRAVEShome} ${mcipPath} ${GDNAM} ${YEAR}
