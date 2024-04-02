% This batch script analyses the resting state fMRI dataset
% available from the SPM website using DCM:
%   http://www.fil.ion.ucl.ac.uk/spm/data/spDCM/
% as described in the SPM manual:
%   http://www.fil.ion.ucl.ac.uk/spm/doc/spm12_manual.pdf#Chap:DCM_rsfmri
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Adeel Razi
% $Id$

clear;
close all;

fs              = '/';       % platform-specific file separator
dir_base        = '';
dir_functional  = 'func'; % base directory of functional scans (Nifti)
dir_struct      = 'anat';

GLM          = 1;
specify_DCM  = 0;
estimate_DCM = 0;

RT = 0.72; % TR or the reptition time

cd(dir_base)
name_subj = dir(); %list all files and folders
dirFlags = [name_subj.isdir]; % logical array for directory (folder) only
name_subj = name_subj(dirFlags); % Extract only directories


TIME = tic;

for s0 = 3 : length(name_subj) 
   
    disp(['Analysing Subject : ', name_subj(s0).name]);
   
    subj_dir = [dir_base fs name_subj(s0).name fs 'MNINonLinear' fs 'Results' fs 'rfMRI_REST1_LR'];

    f = spm_select('FPList', subj_dir, '^rfMRI.*\.nii$');
    mvmnt = sprintf('Movement_Regressors_firstsix.txt');

    clear matlabbatch

  
    if GLM
        glmdir = fullfile(subj_dir,'glm');
        if ~exist(glmdir,'file'), mkdir(glmdir); end
    

        % First GLM specification
        %-----------------------------------------------------------------
        matlabbatch{1}.spm.stats.fmri_spec.dir = {glmdir};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = RT;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(f);
        %     matlabbatch{glm1}.spm.stats.fmri_spec.global = 'Scaling';

        % First GLM estimation
        %-----------------------------------------------------------------
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(glmdir,'SPM.mat')};

        % Extraction of time series from WM and CSF
        %-----------------------------------------------------------------
        matlabbatch{3}.spm.util.voi.spmmat = {fullfile(glmdir,'SPM.mat')};
        matlabbatch{3}.spm.util.voi.adjust = NaN;
        matlabbatch{3}.spm.util.voi.session = 1;
        matlabbatch{3}.spm.util.voi.name = 'CSF';
        matlabbatch{3}.spm.util.voi.roi{1}.sphere.centre = [0 -40 -5]; % 0 -40 -5   
        matlabbatch{3}.spm.util.voi.roi{1}.sphere.radius = 6;
        matlabbatch{3}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
        matlabbatch{3}.spm.util.voi.roi{2}.mask.image = {fullfile(glmdir,'mask.nii')};
        matlabbatch{3}.spm.util.voi.expression = 'i1 & i2';

        matlabbatch{4} = matlabbatch{3};
        matlabbatch{4}.spm.util.voi.name = 'WM';
        matlabbatch{4}.spm.util.voi.roi{1}.sphere.centre = [0 -24 -33]; 
        
        % Second GLM specification
        %-----------------------------------------------------------------
        
        glmdir = fullfile(subj_dir,'glm_corrected');
        if ~exist(glmdir,'file'), mkdir(glmdir); end
        
        matlabbatch{5}.spm.stats.fmri_spec.dir = {glmdir};
        matlabbatch{5}.spm.stats.fmri_spec.timing.units = 'scans';
        matlabbatch{5}.spm.stats.fmri_spec.timing.RT = 0.72;
        matlabbatch{5}.spm.stats.fmri_spec.sess.scans = cellstr(f);
        matlabbatch{5}.spm.stats.fmri_spec.sess.multi_reg = {
            fullfile(subj_dir,mvmnt),...
            fullfile(subj_dir,'glm','VOI_CSF_1.mat'),...
            fullfile(subj_dir,'glm','VOI_WM_1.mat'),...
            }';
        %     matlabbatch{4}.spm.stats.fmri_spec.global = 'Scaling';

        % Second GLM estimation
        %-----------------------------------------------------------------
        matlabbatch{6}.spm.stats.fmri_est.spmmat = {fullfile(glmdir,'SPM.mat')};

        % Extraction of time series from specific 11 nodes
        %-----------------------------------------------------------------
        matlabbatch{7}.spm.util.voi.spmmat = {fullfile(glmdir,'SPM.mat')};
        matlabbatch{7}.spm.util.voi.adjust = NaN;
        matlabbatch{7}.spm.util.voi.session = 1;
        matlabbatch{7}.spm.util.voi.name = 'PCC';
        matlabbatch{7}.spm.util.voi.roi{1}.sphere.centre = [-4 -52 24]; 
        matlabbatch{7}.spm.util.voi.roi{1}.sphere.radius = 8;
        matlabbatch{7}.spm.util.voi.roi{2}.mask.image = {fullfile(glmdir,'mask.nii')};
        matlabbatch{7}.spm.util.voi.expression = 'i1 & i2';

        matlabbatch{8} = matlabbatch{7};
        matlabbatch{8}.spm.util.voi.name = 'mPFC';
        matlabbatch{8}.spm.util.voi.roi{1}.sphere.centre = [-2 54 18];
        matlabbatch{8}.spm.util.voi.expression = 'i1 & i2';
        
        matlabbatch{9} = matlabbatch{7};
        matlabbatch{9}.spm.util.voi.name = 'lHP';
        matlabbatch{9}.spm.util.voi.roi{1}.sphere.centre = [-28 -18 -16]; 
        matlabbatch{9}.spm.util.voi.roi{3}.mask.image = {fullfile(dir_base,'stanfordfROI_mask_Lhippo.nii')};
        matlabbatch{9}.spm.util.voi.expression = 'i1 & i2 & i3';

        matlabbatch{10} = matlabbatch{7};
        matlabbatch{10}.spm.util.voi.name = 'rHP';
        matlabbatch{10}.spm.util.voi.roi{1}.sphere.centre = [28 -18 -16]; 
        matlabbatch{10}.spm.util.voi.roi{3}.mask.image = {fullfile(dir_base,'stanfordfROI_mask_Rhippo.nii')};
        matlabbatch{10}.spm.util.voi.expression = 'i1 & i2 & i3';
        
        matlabbatch{11} = matlabbatch{7};
        matlabbatch{11}.spm.util.voi.name = 'lAMG';
        matlabbatch{11}.spm.util.voi.roi{1}.sphere.centre = [-24 -2 -22]; 
        matlabbatch{11}.spm.util.voi.roi{3}.mask.image = {fullfile(dir_base,'wfu_mask_LAmygdala.nii')};
        matlabbatch{11}.spm.util.voi.expression = 'i1 & i2 & i3';
        
        matlabbatch{12} = matlabbatch{7};
        matlabbatch{12}.spm.util.voi.name = 'rAMG';
        matlabbatch{12}.spm.util.voi.roi{1}.sphere.centre = [24 -2 -22]; 
        matlabbatch{12}.spm.util.voi.roi{3}.mask.image = {fullfile(dir_base,'wfu_mask_RAmygdala.nii')};
        matlabbatch{12}.spm.util.voi.expression = 'i1 & i2 & i3';
        
        matlabbatch{13} = matlabbatch{7};
        matlabbatch{13}.spm.util.voi.name = 'dACC';
        matlabbatch{13}.spm.util.voi.roi{1}.sphere.centre = [4 32 22]; 
        matlabbatch{13}.spm.util.voi.expression = 'i1 & i2';
        
        matlabbatch{14} = matlabbatch{7};
        matlabbatch{14}.spm.util.voi.name = 'lAI'; 
        matlabbatch{14}.spm.util.voi.roi{1}.sphere.centre = [-34 22 0]; 
        matlabbatch{14}.spm.util.voi.roi{3}.mask.image = {fullfile(dir_base,'wfu_mask_LINS.nii')};
        matlabbatch{14}.spm.util.voi.expression = 'i1 & i2 & i3';
        
        matlabbatch{15} = matlabbatch{7};
        matlabbatch{15}.spm.util.voi.name = 'rAI';
        matlabbatch{15}.spm.util.voi.roi{1}.sphere.centre = [34 22 0]; 
        matlabbatch{15}.spm.util.voi.roi{3}.mask.image = {fullfile(dir_base,'wfu_mask_RINS.nii')};
        matlabbatch{15}.spm.util.voi.expression = 'i1 & i2 & i3';
        
        matlabbatch{16} = matlabbatch{7};
        matlabbatch{16}.spm.util.voi.name = 'lDLPFC';
        matlabbatch{16}.spm.util.voi.roi{1}.sphere.centre = [-46 34 32]; 
        matlabbatch{16}.spm.util.voi.expression = 'i1 & i2';
        
        matlabbatch{17} = matlabbatch{7};
        matlabbatch{17}.spm.util.voi.name = 'rDLPFC';
        matlabbatch{17}.spm.util.voi.roi{1}.sphere.centre = [42 38 32]; 
        matlabbatch{17}.spm.util.voi.expression = 'i1 & i2';
        


        
        spm_jobman('run',matlabbatch);
    end  

        % DCM specification
        %--------------------------------------------------------------
    if specify_DCM
        
        glmdir = fullfile(subj_dir,'glm_corrected');
        cd(glmdir)
        load('VOI_PCC_1.mat');
        DCM.xY(1) = xY;
        load('VOI_mPFC_1.mat');
        DCM.xY(2) = xY;
        load('VOI_lHP_1.mat');
        DCM.xY(3) = xY;
        load('VOI_rHP_1.mat');
        DCM.xY(4) = xY;
        load('VOI_lAMG_1.mat');
        DCM.xY(5) = xY;
        load('VOI_rAMG_1.mat');
        DCM.xY(6) = xY;
        load('VOI_dACC_1.mat');
        DCM.xY(7) = xY;
        load('VOI_lAI_1.mat');
        DCM.xY(8) = xY;
        load('VOI_rAI_1.mat');
        DCM.xY(9) = xY;
        load('VOI_lDLPFC_1.mat');
        DCM.xY(10) = xY;
        load('VOI_rDLPFC_1.mat');
        DCM.xY(11) = xY;
 

        DCM.v = length(DCM.xY(1).u); % number of time points
        DCM.n = length(DCM.xY);      % number of regions
        DCM.Y.dt  = RT;
        DCM.Y.X0  = DCM.xY(1).X0;

        for i = 1:DCM.n
            DCM.Y.y(:,i)  = DCM.xY(i).u;
            DCM.Y.name{i} = DCM.xY(i).name;
        end

        Y  = DCM.Y;                             % responses
        v  = DCM.v;                             % number of scans
        n  = DCM.n;                             % number of regions

        DCM.Y.Q    = spm_Ce(ones(1,n)*v);
        DCM.U.u    = zeros(v,1);
        DCM.U.name = {'null'};         

        DCM.a = ones(n,n);
        DCM.b  = zeros(n,n,0);
        DCM.c  = zeros(n,0);
        DCM.d = zeros(n,n,0);

        DCM.TE     = 0.0331;
        DCM.delays = repmat(RT,DCM.n,1);

        DCM.options.nonlinear  = 0;
        DCM.options.two_state  = 0;
        DCM.options.stochastic = 0;
        DCM.options.analysis   = 'CSD';
        DCM.options.induced    = 1;
        DCM.options.maxnodes       = 16; % number of modes that we want to use
        %     DCM.options.nograph    = 1;


        str = sprintf('DCM_%s',name_subj(s0).name);
        DCM.name = str;
        save(fullfile(glmdir,str),'DCM');

    end

    % DCM estimation
    %-----------------------------------------------------------------------
    if estimate_DCM
        glmdir = fullfile(subj_dir,'glm_corrected');
        cd(glmdir);
        str = sprintf('DCM_%s',name_subj(s0).name);
        DCM = spm_dcm_fmri_csd(fullfile(glmdir,str));
        save(fullfile(glmdir,str),'DCM');
       
    end
end



fprintf('This took %s', duration([0, 0, toc(TIME)]));