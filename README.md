# OBIS

This is a set of scripts that is run in three stages to generate the full set of data available from OBIS:
1. A list of all available taxon IDs is generated with generate_species_list.py
2. The species occurrence data are downloaded in json format.
3. The data are binned to the 1x1 degree World Ocean Atlas grid with 33 depth levels and are matched with monthly T/S/O2 data from WOA for all data points that include depth data. Output files are created in mat format.

Notes:
The pyobis package is used for the OBIS API calls. (It needs to be revised for v3 of the OBIS API.)

The scripts contain some hard-coded paths that need to be adapted for a local environment.
