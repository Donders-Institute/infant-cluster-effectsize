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
        scripts     = '/Volumes/Samsung T3/data/BeeG-analysis';
        sourcedata  = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/sourcedata';
        bidsroot    = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/bids';
        results     = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/results';
        derivatives = '/Volumes/Samsung T3/data/di.dcc.DSC_2020.00134_473/derivatives';
    otherwise
        error('directories for data could not be found');
end

addpath(scripts)
