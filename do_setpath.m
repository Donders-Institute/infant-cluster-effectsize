%% Set the path for these scripts and for the data, depending on who is executing the code
%
% These scripts and the data in BIDS format are part of Meyer, M., Lamers, D., Kayhan,
% E., Hunnius, S., & Oostenveld, R. (2021). Enhancing reproducibility in developmental
% EEG research: BIDS, cluster-based permutation tests, and effect sizes (in preparation).
%
% The infant EEG dataset is originally described in Kayhan, E., Meyer, M., O'Reilly,
% J. X., Hunnius, S., & Bekkering, H. (2019). Nine-month-old infants update their
% predictive models of a changing environment. Developmental cognitive neuroscience,
% 38, 100680. https://doi.org/10.1016/j.dcn.2019.100680

user = getenv('USER');
if isempty(user)
  user = getenv('UserName');
end

switch user
  case 'Didi'
    scripts     = 'C:\Users\Didi\Documents\GitHub\BeeG analysis';
    sourcedata  = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset';
    bidsroot    = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\BIDS';
    results     = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\results';
    derivatives = 'C:\Users\Didi\Documents\GitHub\Donders Datasets\BeeG dataset\derivatives';
  case 'roboos'
    scripts     = '/Volumes/SamsungT3/data/di.dcc.DSC_2020.00134_473/scripts';
    sourcedata  = '/Volumes/SamsungT3/data/di.dcc.DSC_2020.00134_473/sourcedata';
    bidsroot    = '/Volumes/SamsungT3/data/di.dcc.DSC_2020.00134_473/bidsdata';
    results     = '/Volumes/SamsungT3/data/di.dcc.DSC_2020.00134_473/results';
    derivatives = '/Volumes/SamsungT3/data/di.dcc.DSC_2020.00134_473/bidsresults';
  case 'Marlene'
    scripts     = 'F:\Documents\Donders\S_ClusterBasedMethods\Analysis\BeeG-analysis';
    sourcedata  = 'F:\Documents\Donders\S_ClusterBasedMethods\Analysis\Sourcedata';
    bidsroot    = 'F:\Documents\Donders\S_ClusterBasedMethods\Analysis\BIDS version 1';
    results     = 'F:\Documents\Donders\S_ClusterBasedMethods\Analysis\BeeG dataset results';
    derivatives = 'F:\Documents\Donders\S_ClusterBasedMethods\Analysis\bidsresults';
  otherwise
    error('The directories for the input and output data could not be found');
end

addpath(scripts)
