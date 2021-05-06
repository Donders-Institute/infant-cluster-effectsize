
These scripts and the data in BIDS format are part of Meyer, M., Lamers, D., Kayhan,
E., Hunnius, S., & Oostenveld, R. (2021) Fostering reproducibility in developmental
EEG research by using BIDS, cluster-based permutation tests and reporting
effectsizes (in preparation)

The infant EEG dataset is originally described in Kayhan, E., Meyer, M., O'Reilly,
J. X., Hunnius, S., & Bekkering, H. (2019). Nine-month-old infants update their
predictive models of a changing environment. Developmental cognitive neuroscience,
38, 100680.)

The complete EEG dataset is available from https://doi.org/10.34973/gvr3-6g88

The MATLAB code should be executed as follows

- do_setpath
- do_convert_data_to_BIDS (only to be done once)
- do_complete_analysis, this will call
  - do_singlesubject_analysis for each subject
  - do_group_analysis
- do_convert_results_to_BIDS
