# wfc3-telemetry

Searches for telemetry correlations with all HST telemetry paramters.

The following outputs are created (note: 2 and 3 are only created if the user inputs a 
table containing both the MJDs and the parameter of interest tied to those MJDs - see 
the 'Use' section below):

1. telemetry_table.txt : table showing the values of all HST telemetry parameters at the
                         given MJDs
2. correlations.txt : the ranked telemetry correlations with the user-defined parameter
3. correlation_matrix.png : correlation matrix of the top telemetry correlations from
                            correlations.txt
    
Use
---
    This script can be run via the command line as such:
        
        python wfc3_telemetry.py --t <t>
        
    --t [Required]: The name of the table containing the MJDs of interest. This table should
                    take the form:
                    
                    mjd param
                    55555.1 23
                    55556.5 24.5
                    55560 -17.22
                    
                    If only the mjd column is included, the only output that is created is
                    telemetry_table.txt.
