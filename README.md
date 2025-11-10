# CAT12-T1-MRI-data-preprocessing
This repository stored scripts for T1 MRI data preprocessing using CAT12 toolbox in SPM
Overview

This pipeline performs structural MRI (T1-weighted) preprocessing using the CAT12 toolbox within SPM (Statistical Parametric Mapping).
It prepares T1-weighted anatomical images for further voxel-based morphometry (VBM), cortical thickness analysis, or other morphometric studies.

Main preprocessing steps include:

- Bias field correction

- Segmentation into GM, WM, and CSF

- Spatial normalization to MNI space

- Modulation (optional)

- Smoothing (optional)

## Requirements
- Software

- MATLAB (R2020a or later recommended)

- SPM12 (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)

- CAT12 Toolbox (https://neuro.uni-jena.de/cat/)

## Preprocessing Steps
### 1. Launch MATLAB and Add Paths
'''
addpath('/path/to/spm12');
addpath('/path/to/cat12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');
'''

### 2. Run CAT12 Batch Script

Example MATLAB script (scripts/preprocess_cat12.m), which you can find it in the file.

### 3. Outputs

After running CAT12, you will find the following outputs in the /output/directory/mri/ directory:

- mwp1*	% Modulated normalized gray matter
- mwp2*	% Modulated normalized white matter
- mwcsf*	% Modulated normalized CSF
- y_*	% Deformation field
- p*	% Quality control and parameter files
- report/	%CAT12 quality assurance reports

### 4. Quality Control

CAT12 automatically generates a report (catreport_sub-XX.pdf) for each subject, containing:

- Image quality metrics

- Segmentation results

- Surface reconstruction check

You can also inspect the normalized GM maps using SPM’s “Check Reg” function:
'''
spm_check_registration('mwp1sub-01_T1w.nii');

'''
### 5. Optional: Group-Level Analysis (VBM)

After preprocessing, you can use the smoothed modulated GM images (smwp1*.nii) for VBM analysis.

Example scripts included three main steps: 1)sooth the GM files (mwp1*.nii);2) build the design; 3) estimate the model, and 4) create common contrasts — with optional covariates (Age, Sex, TIV).

'''
% ================================================================
% VBM TWO-GROUP ANALYSIS (CAT12/SPM12)
% Author: <you>
% Date: <today>
%
% Pipeline:
%   1) Read subject table (CSV)
%   2) Ensure smoothed, modulated GM images exist (optional smoothing)
%   3) Build two-sample t-test factorial design (+ optional covariates)
%   4) Estimate model
%   5) Define contrasts
%   6) (Optional) write out a design summary
%
% Requires:
%   - MATLAB
%   - SPM12 on path
%   - (Recommended) CAT12 on path (preprocessing already done)
% ================================================================

%% ===================== CONFIG =====================
% --- Paths ---
SPM_DIR = '/path/to/spm12';         % e.g., '/opt/spm12'
CAT12_DIR = '/path/to/cat12';       % optional but recommended
CSV_FILE = '/path/to/subjects.csv'; % see CSV FORMAT section below
OUT_DIR  = '/path/to/derivatives/vbm_stats';

% --- Images & smoothing ---
% Column in CSV that points to *modulated*, MNI-normalized GM images.
% This can be 'mwp1' (unsmoothed) or already-smoothed 'smwp1'.
IMG_COLUMN_NAME = 'gm_path';   % column with full path to each subject's GM image
DO_SMOOTH   = true;            % set false if your images are already smoothed
SMOOTH_FWHM = [8 8 8];         % typical VBM kernel (in mm)
SMOOTH_PREFIX = 'sm';          % output prefix if smoothing

% --- Groups ---
% Two-sample t-test with GROUP values 1 and 2 (e.g., 1=Controls, 2=Patients)
GROUP_COLUMN_NAME = 'group';   % numeric: 1 or 2 (preferred)
GROUP1_LABEL = 'Controls';
GROUP2_LABEL = 'Patients';

% --- Optional covariates ---
% If present as columns in CSV, they will be added (mean-centered overall).
COVS = { ...
  struct('name','Age', 'column','age', 'add',true), ...
  struct('name','Sex', 'column','sex', 'add',true), ... % code sex as 0/1 or 1/2
  struct('name','TIV', 'column','tiv', 'add',true) ...
};

% --- Masking (optional explicit mask) ---
% Leave empty '' to use implicit masking only
EXPLICIT_MASK = '';  % e.g., '/path/to/group_mask.nii'

% --- QA/robustness ---
ABORT_IF_MISSING_IMAGE = true;  % stop if any image path invalid
WRITE_DESIGN_SUMMARY   = true;  % write text summary of design & covs

% --- Contrasts to create ---
MAKE_GROUP_MEANS = true;  % also create [1 0] and [0 1]
% ================================================================

%% ===================== SETUP =====================
addpath(SPM_DIR);
if exist(CAT12_DIR, 'dir'), addpath(CAT12_DIR); end
spm('defaults','FMRI');
spm_jobman('initcfg');

if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end

% Read CSV (uses readtable; expects headers)
T = readtable(CSV_FILE, 'TextType', 'string');
assert(all(ismember({IMG_COLUMN_NAME, GROUP_COLUMN_NAME}, T.Properties.VariableNames)), ...
       'CSV must contain columns "%s" and "%s".', IMG_COLUMN_NAME, GROUP_COLUMN_NAME);

% Coerce required columns
img_paths = string(T.(IMG_COLUMN_NAME));
groups    = double(T.(GROUP_COLUMN_NAME));

% Filter valid subjects (groups 1 or 2, and image exists)
is_g12 = ismember(groups, [1 2]);
exists_img = arrayfun(@(p) isfile(char(p)), img_paths);
valid = is_g12 & exists_img;

if any(~exists_img)
    missing = find(~exists_img);
    warning('Missing image for %d subject(s). Example: %s', numel(missing), char(img_paths(missing(1))));
    if ABORT_IF_MISSING_IMAGE, error('Aborting due to missing images.'); end
end
if ~any(valid), error('No valid subjects with groups {1,2} and existing images.'); end

T = T(valid, :);
img_paths = img_paths(valid);
groups    = groups(valid);

fprintf('N=%d valid subjects (Group1=%d, Group2=%d)\n', ...
    height(T), sum(groups==1), sum(groups==2));

%% ===================== SMOOTHING (optional) =====================
if DO_SMOOTH
    fprintf('Smoothing %d images with FWHM = [%g %g %g]...\n', numel(img_paths), SMOOTH_FWHM);
    smoothed_paths = strings(size(img_paths));
    for i = 1:numel(img_paths)
        inFile = char(img_paths(i));
        [p,n,e] = fileparts(inFile);
        outFile = fullfile(p, [SMOOTH_PREFIX n e]);
        if ~isfile(outFile)
            spm_smooth(inFile, outFile, SMOOTH_FWHM);
        end
        smoothed_paths(i) = string(outFile);
    end
    img_paths = smoothed_paths;
else
    % check they look smoothed already
    looks_smoothed = startsWith(cellstr(string(img_paths)), SMOOTH_PREFIX);
    if ~all(looks_smoothed)
        warning('DO_SMOOTH==false but some images do not start with "%s". Proceeding anyway.', SMOOTH_PREFIX);
    end
end

% Sanity: split by group
scans1 = cellstr(img_paths(groups==1));
scans2 = cellstr(img_paths(groups==2));

assert(~isempty(scans1) && ~isempty(scans2), 'Each group must have at least one subject.');

%% ===================== COVARIATES =====================
cov_defs = struct('name',{},'values',{},'iCC',{});
for k = 1:numel(COVS)
    cfg = COVS{k};
    if cfg.add && ismember(cfg.column, T.Properties.VariableNames)
        v = double(T.(cfg.column));
        if any(isnan(v))
            warning('Covariate %s has NaNs; removing rows with NaN in that covariate.', cfg.name);
            keep = ~isnan(v);
            T = T(keep, :);
            img_paths = img_paths(keep);
            groups = groups(keep);
            v = v(keep);
            % recompute scans per group
            scans1 = cellstr(img_paths(groups==1));
            scans2 = cellstr(img_paths(groups==2));
        end
        % mean-center overall (SPM can also center; we’ll pass raw and set iCC=1)
        cov_defs(end+1).name   = cfg.name; %#ok<SAGROW>
        cov_defs(end).values   = v(:);
        cov_defs(end).iCC      = 1;  % overall mean
    elseif cfg.add
        warning('Covariate column "%s" not found in CSV; skipping %s.', cfg.column, cfg.name);
    end
end

% Consistency: covariate length equals total subjects
N = numel(img_paths);
for k = 1:numel(cov_defs)
    assert(numel(cov_defs(k).values)==N, 'Covariate %s length mismatch.', cov_defs(k).name);
end

%% ===================== FACTORIAL DESIGN =====================
matlabbatch = {};

matlabbatch{1}.spm.stats.factorial_design.dir = {OUT_DIR};
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = scans1';
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = scans2';

% Independence & variance (defaults are usually fine for two-sample VBM)
% i = independence; v = variance
matlabbatch{1}.spm.stats.factorial_design.des.t2.ind = 1;  % groups independent
matlabbatch{1}.spm.stats.factorial_design.des.t2.var = 1;  % assume unequal variance (Welch) often safer

% Covariates
if isempty(cov_defs)
    matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
else
    for k = 1:numel(cov_defs)
        matlabbatch{1}.spm.stats.factorial_design.cov(k).c      = cov_defs(k).values;
        matlabbatch{1}.spm.stats.factorial_design.cov(k).cname  = cov_defs(k).name;
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCFI   = 1;   % no interaction
        matlabbatch{1}.spm.stats.factorial_design.cov(k).iCC    = cov_defs(k).iCC; % overall mean
    end
end

% Masking
if ~isempty(EXPLICIT_MASK)
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {EXPLICIT_MASK};
else
    matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
end
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;  % implicit mask
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;

% Global normalization (omit; modulation already accounts for size)
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

%% ===================== MODEL ESTIMATION =====================
matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(OUT_DIR, 'SPM.mat')};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

%% ===================== CONTRASTS =====================
matlabbatch{3}.spm.stats.con.spmmat = {fullfile(OUT_DIR,'SPM.mat')};
ci = 1;

% Two-sample: [1 -1] = Group1>Group2, [-1 1] = Group2>Group1
matlabbatch{3}.spm.stats.con.consess{ci}.tcon.name    = sprintf('%s > %s', GROUP1_LABEL, GROUP2_LABEL);
matlabbatch{3}.spm.stats.con.consess{ci}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{ci}.tcon.sessrep = 'none';
ci = ci + 1;

matlabbatch{3}.spm.stats.con.consess{ci}.tcon.name    = sprintf('%s > %s', GROUP2_LABEL, GROUP1_LABEL);
matlabbatch{3}.spm.stats.con.consess{ci}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{ci}.tcon.sessrep = 'none';
ci = ci + 1;

if MAKE_GROUP_MEANS
    matlabbatch{3}.spm.stats.con.consess{ci}.tcon.name    = sprintf('%s mean', GROUP1_LABEL);
    matlabbatch{3}.spm.stats.con.consess{ci}.tcon.weights = [1 0];
    matlabbatch{3}.spm.stats.con.consess{ci}.tcon.sessrep = 'none';
    ci = ci + 1;

    matlabbatch{3}.spm.stats.con.consess{ci}.tcon.name    = sprintf('%s mean', GROUP2_LABEL);
    matlabbatch{3}.spm.stats.con.consess{ci}.tcon.weights = [0 1];
    matlabbatch{3}.spm.stats.con.consess{ci}.tcon.sessrep = 'none';
    ci = ci + 1;
end

matlabbatch{3}.spm.stats.con.delete = 0;

%% ===================== RUN =====================
spm_jobman('run', matlabbatch);

%% ===================== SUMMARY FILE (optional) =====================
if WRITE_DESIGN_SUMMARY
    sumfile = fullfile(OUT_DIR,'design_summary.txt');
    fid = fopen(sumfile,'w');
    fprintf(fid, 'VBM two-sample t-test\n');
    fprintf(fid, 'Output dir: %s\n', OUT_DIR);
    fprintf(fid, 'N total: %d | %s: %d | %s: %d\n', ...
        N, GROUP1_LABEL, numel(scans1), GROUP2_LABEL, numel(scans2));
    fprintf(fid, 'Smoothing: %s | FWHM = [%g %g %g]\n', string(DO_SMOOTH), SMOOTH_FWHM);
    if ~isempty(EXPLICIT_MASK)
        fprintf(fid, 'Explicit mask: %s\n', EXPLICIT_MASK);
    else
        fprintf(fid, 'Explicit mask: (none)\n');
    end
    if ~isempty(cov_defs)
        fprintf(fid, 'Covariates:\n');
        for k = 1:numel(cov_defs)
            fprintf(fid, '  - %s (mean-centered overall)\n', cov_defs(k).name);
        end
    else
        fprintf(fid, 'Covariates: (none)\n');
    end
    fprintf(fid, '\nFirst few scans by group:\n');
    for i = 1:min(5, numel(scans1)), fprintf(fid, 'G1: %s\n', scans1{i}); end
    for i = 1:min(5, numel(scans2)), fprintf(fid, 'G2: %s\n', scans2{i}); end
    fclose(fid);
    fprintf('Wrote %s\n', sumfile);
end

fprintf('Done. Open SPM -> Results to review thresholded maps.\n');

'''
#### CSV FORMAT (example)

Save as subjects.csv (headers required).

gm_path: full path to the modulated GM image for each subject. Use your CAT12 output (e.g., mwp1sub-XX_T1w.nii). If you let the script smooth, it will write smwp1*.nii next to it and use that.

group: 1 for Controls, 2 for Patients (you can rename labels in the config).

Optional columns: age, sex (0/1 or 1/2), tiv (Total Intracranial Volume from CAT12 reports or export).

'''
subject_id,gm_path,group,age,sex,tiv
sub-01,/proj/derivatives/cat12/sub-01/mwp1sub-01_T1w.nii,1,24,1,1536000
sub-02,/proj/derivatives/cat12/sub-02/mwp1sub-02_T1w.nii,1,26,0,1490000
sub-10,/proj/derivatives/cat12/sub-10/mwp1sub-10_T1w.nii,2,25,1,1512000
sub-11,/proj/derivatives/cat12/sub-11/mwp1sub-11_T1w.nii,2,29,0,1483000
'''

### References

- Ashburner J, Friston KJ. SPM12 Manual. Wellcome Centre for Human Neuroimaging.

- Gaser C, Dahnke R. (2016). CAT – A Computational Anatomy Toolbox for the Analysis of Structural MRI Data.
Hbm 2016, 336–348.

- https://neuro.uni-jena.de/cat/
