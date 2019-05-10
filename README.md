## Deletions pipeline

### Running the pipeline

Update `config.sh`, `included_samples_info.txt`. Input file `processing.txt` should contain list of all samples with the processing stage (starts at 0, ends at 12).

Steps:

1. Steps 1-7 of the pipeline: `pipeline-stage1.sh`
2. Mapping quality: check, update if needed and run `MAPQ_slurm_submit_sep.sh` for new samples
3. Number of SNPs: check, update if needed and run `hadoop_query.py` at hdp-master2 for all samples
4. Steps 8-12 of the pipeline: `pipeline-stage2.sh`

### Information about pipeline and output

See [README.txt](../blob/master/README.txt) for description and details.
