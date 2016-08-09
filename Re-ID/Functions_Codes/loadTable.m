%
% Load the average pairwise distance
% Written by Paul for Hamid
% 6 Feb 2015
%

function distAvg = loadTable(dataset, trialNo)

distAvg = [];
if nargin ~= 2, printUsage(); return; end

if ~ (isequal(dataset, '3dpes') || isequal(dataset, 'cuhk01') || ...
      isequal(dataset, 'cuhk03') || isequal(dataset, 'ilids') || ...
      isequal(dataset, 'viper'))
      printUsage(); return;
end    

if trialNo < 1 || trialNo > 10, printUsage(); return; end

load(sprintf('%s/LFDA_%02d.mat', dataset, trialNo));

distAll = cat(3,dist.single{:});
distAvg = mean(distAll,3);
% cmc = evaluate_pwdist(distAvg');
% fprintf('Rank-1 Accuracy: %f\n', cmc(1));


function printUsage
fprintf('Usage: table = loadTable(DATASET, trialNo)\n');
fprintf('       where DATASET = {3dpes, cuhk01, cuhk03, ilids, viper}\n');
fprintf('             trialNo is between 1 to 10\n\n');
