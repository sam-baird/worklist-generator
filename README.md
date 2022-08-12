# Worklist Generator

Select COVID-19 samples for sequencing based on RT-qPCR results. For use at the Colorado State Public Health Laboratory Genomic Surveillance Unit.

## Usage

Ensure that the Numpy and Pandas Python packages are installed.

Copy the PCR result CSV files (they must start with 'nCoV...') to the same directory as the `worklist.py` script. Only include files for batches that are not on the previous BatchSmash file. Copy the previous BatchSmash to the current directory. The last COVSEQ ID in the BatchSmash file should correspond to the first COVSEQ ID to be worklisted.

Run the script using `python worklist.py`. Any samples with a Ct value less than 30 (but not '0') for the N gene and ORF1ab targets will be selected for sequencing. This value can be changed by modifying the `CT_CUTOFF` constant in the script.

An `Ouput` directory will be created, which will contain the new BatchSmash, UNSAT, and Freedom Evo Tecan files.
