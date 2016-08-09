function distAvg = loadTable2(dataset, trialNo)

distAvg = [];
if nargin ~= 2, printUsage(); return; end

if ~ (isequal(dataset, 'RAiD') || isequal(dataset, 'WARD')) 
      printUsage(); return;
end    

if trialNo < 1 || trialNo > 5, printUsage(); return; end

load(sprintf('%s/FT_%02d.mat', dataset, trialNo));
score_mt= max(score_mt,0.0000001);
distAvg = -log(score_mt);
% cmc = evaluate_pwdist(distAvg');
% fprintf('Rank-1 Accuracy: %f\n', cmc(1));


function printUsage
fprintf('Usage: table = loadTable(DATASET, trialNo)\n');
fprintf('       where DATASET = {RAiD, WARD}\n');
fprintf('             trialNo is between 1 to 5\n\n');