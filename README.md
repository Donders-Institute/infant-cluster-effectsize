These scripts and the data in BIDS format are part of Meyer, M., Lamers, D., Kayhan,
E., Hunnius, S., & Oostenveld, R. (2021). Enhancing reproducibility in developmental
EEG research: BIDS, cluster-based permutation tests, and effect sizes, Developmental 
Cognitive Neuroscience, 52, 101036, https://doi.org/10.1016/j.dcn.2021.101036.

The infant EEG dataset is originally described in Kayhan, E., Meyer, M., O'Reilly,
J. X., Hunnius, S., & Bekkering, H. (2019). Nine-month-old infants update their
predictive models of a changing environment. Developmental cognitive neuroscience,
38, 100680, https://doi.org/10.1016/j.dcn.2019.100680.

The raw EEG dataset is available at https://doi.org/10.34973/gvr3-6g88.

The processed EEG data (in MATLAB and BIDS format) is available at https://doi.org/10.34973/g4we-5v66.

The MATLAB code should be executed as follows

- do_setpath
- do_convert_data_to_BIDS (Only to be done once; It refers to function 'trialfun_BeeG'; The code in this script is referred to as Script Section 1 in the manuscript)
- do_prepare_neighbours   (Only to be done once; The code in this script is not meant to be executed as-is, but is provided for reference)
- do_complete_analysis, this will call
  - do_singlesubject_analysis for each subject (The code in this script is referred to as Script Section 2 in the manuscript)
  - do_group_analysis (The code in this script is referred to as Script Section 3 in the manuscript)
- do_convert_results_to_BIDS

Note that the original source data referred to in the do_convert_data_to_BIDS file 
