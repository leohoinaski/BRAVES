-----------------------------BRAVESdatabase-----------------------------

This repository contains all functions and data required to create your
own vehicular emissions files from BRAVES outputs. You can setup the grid
definition, output names, time window, and other features. You can choose
to create netCDF files for a single chemical specie or full pack ready to 
use in CMAQ.

If you want to run BRAVESdatabase you just need to run the script 
BRAVESdatabase_main.py. You can adjust this scrip to create your database. 

This repository also contains:

FOLDERS:
	ChemicalSpec folder: folder with chemical speciation profiles, 
			molecular weight and specie names.

	TemporalAloc: folder with temporal profiles for temporal 
			disaggregation.

FUNCTIONS:
	BRAVES_temporalDisag_v1.py : function used to create temporal 
		disaggregated emissions.
	
	roadDensity_v1.py: function used to calculate the road density
		from openstreetmaps shapefiles.

	mergeRoadEmiss_v1.py: function used to merge spatial disaggre-
		gated files from roadEmiss_v1.py. 
		
	roadEmiss_v1.py: generates the disaggregated vehicular emissions
		using BRAVES and roadDensity_v1.py outputs. 

	netCDFcreator_v1.py: creates the netCDF files

	BRAVES_ChemicalSpec_v1.py: generates the chemical speciated from 
		mergeRoadEmiss_v1.py data

	BRAVES2netCDF_v1.py: controls the creation of netCDF files

OUTPUTS:

	The BRAVESdatabase_main.py will create the Outputs folder in your 
	computer with intermediate and final outputs. 

REQUIRED DATA:

	You can download the openstreetmaps files and BRAVES outputs data 
	at: https://arquivos.ufsc.br/d/c04a184981a94dab87e9/ and 
	https://arquivos.ufsc.br/d/7f126a12c9f04b228031/

	The folders "Shapefiles" and BRAVESoutputs must be in the same 
	directory as the BRAVESdatabase_main.py. You should keep this 
	names for this folders. 

	Before you run the BRAVESdatabase_main.py script you should install
	the requires python packages in requirements.txt. You can do this by
	pip install -r requirements.txt


	See https://stackoverflow.com/questions/64038673/could-not-build-wheels-for-which-use-pep-517-and-cannot-be-installed-directly
	sudo apt-get install libproj-dev proj-data proj-bin  
	sudo apt-get install libgeos-dev  
	sudo apt-get install libgdal-dev

