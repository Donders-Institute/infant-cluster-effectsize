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
    case 'U567154'
        scripts     = '\\cnas.ru.nl\u567154\Documents\Donders\S_ClusterBasedMethods\Analysis\BeeG-analysis';
        sourcedata  = '\\CNAS.RU.NL\U567154\Documents\Donders\S_ClusterBasedMethods\Analysis\Sourcedata';
        bidsroot    = '\\CNAS.RU.NL\U567154\Documents\Donders\S_ClusterBasedMethods\Analysis\BIDS version 1';
        results     = '\\cnas.ru.nl\u567154\Documents\Donders\S_ClusterBasedMethods\Analysis\BeeG dataset results';
       %derivatives = '';
    otherwise
        error('directories for data could not be found');
end

addpath(scripts)
