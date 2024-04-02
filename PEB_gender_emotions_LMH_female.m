clear all;
close all;
%%

GCM_dir = 'F:\Dua\NED_BME\FYDP\Data\Anger_female'; % main folder hcp 1079
 
cd(GCM_dir)

GCM_hcp1079_female = cellstr((spm_select('FPListRec', fullfile(GCM_dir,'Medium_903'), '^DCM.*\.mat$')));% full path
Anger_female_moderate = importdata('F:\Dua\NED_BME\FYDP\HCP\Text files_903\Scores\Modangscore_female.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First provide estimated DCM as GCM (Group DCM) cell array. Individual DCMs can be
% estimated by using spm_dcm_fit.m

M = struct();
M.alpha = 1;
M.beta  = 16;
M.hE    = 0;
M.hC    = 1/16;
% M.Q     = 'all';
M.Q = 'single';

N1 = length(GCM_hcp1079_female);

%Mean

Anger_female_moderate = Anger_female_moderate - mean(Anger_female_moderate);

%% Design Matrix
% Specify design matrix for N subjects. It should start with a constant column

% Choose field 
field = {'A'};

% For AngAffect (female)
X1 = [ones(N1,1) Anger_female_moderate];
M.X = X1;
PEB = spm_dcm_peb(GCM_hcp1079_female,M,field);
BMA = spm_dcm_peb_bmc(PEB);

save('PEB_hcp903_female_AngerAffect_moderate.mat','BMA','PEB', 'GCM_hcp1079_female', 'X1');

spm_dcm_peb_review(BMA)

