clear all;
close all;
%%

GCM_dir = 'F:\Dua\NED_BME\FYDP\Data\Anger_male'; % main folder hcp 1079
 
cd(GCM_dir)


GCM_hcp1079_male = cellstr((spm_select('FPListRec', fullfile(GCM_dir,'Medium'), '^DCM.*\.mat$')));
Anger_male_mod = importdata('F:\Dua\NED_BME\FYDP\HCP\Text files_903\Scores\Modangscore_male.txt');

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

N2 = length(GCM_hcp1079_male);

%Mean

Anger_male_mod = Anger_male_mod - mean(Anger_male_mod);
%% Design Matrix
% Specify design matrix for N subjects. It should start with a constant column

% Choose field 
field = {'A'};

% For AngAffect (male)
X1 = [ones(N2,1) Anger_male_mod];
M.X = X1;
PEB = spm_dcm_peb(GCM_hcp1079_male,M,field);
BMA = spm_dcm_peb_bmc(PEB);


save('PEB_hcp903_male_AngerAffect_moderate.mat','BMA','PEB', 'GCM_hcp1079_male', 'X1');

spm_dcm_peb_review(BMA)

