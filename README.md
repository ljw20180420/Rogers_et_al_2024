# This is a permanent fork.

There is no plan for pr, so this is a permanent fork of git@github.com:sarogers9/Rogers_et_al_2024.git. This repo fixes some bug. This does not mean the original code base and paper are wrong. The author does some "trials and errors" in the script and did not (un)comment or modify everything back to make the script run smoothly. This repo also makes the code more readable and tidy.

# Rogers_et_al_2024
Original data processing pipeline for all data (raw dF/F traces, TCA model outputs) presented in Figures 2-8 and Supp Figures 2-7 Rogers et al 2024. Wrapper, custom functions, and data files.

## Table of contents
- [Data](#data)
- [Custom functions](#custom-functions)
- [Code](#code)
- [Data viz](#data-viz)
## Data
Data types are .txt, .csv, .mat. The data structure rogers2024.mat contains all datasets, group, and task/trial information. Load into the beginning of all subsequent code.
*ORIGINAl DATA FILES TOO BIG TO SHARE HERE* To download the data structure please use [this google drive link](https://drive.google.com/drive/folders/1naMBzlREXzaXH6uyjMA6XUr7c9STwVQ5?usp=sharing) and request permissions from author to access until further notice. Same documentation applies. 
### 1. Psilocybin TFC behavior
* data for figure 1 (no miniscope), files organized by cohort (.txt files)
* percent time freezing for each animal in each session (hab, acq, ext1, ext2, ext3) `R{round number}B{box number}{session}.csv`
* rows are trial, columns are trial period (baseline, stimulus, trace, ITI)
### 2. Behavior
* files organized by miniscope treatment group
* 15 Hz ezTrack output of freezing status of mouse per frame `M{animal number}_{session}.csv`
### 3. dFoF
* files organized by miniscope treatment group
* deconvolved 20 Hz single cell calcium data for each animal in each session `MM{animal number}{session}.csv`
* longitudinal registration of neurons for each animal `MM{animal number}long.csv`
* cell properties (spatial coordinates and size) of neurons for each animal for each session `MM{animal number}{session-props}.csv`
### 4. Cell pooled, session concatenated activity matrices
* files named by treatment (Psil, Sal,or Nonshock) `totApool{treamtent}.csv'
* output of section 1 the code, used in sections 13+
* 1 second-summed dFoF in each 60 second trial in each session (rows) of each longitudinall registered cell (columns) in a treatment group. 
### 5. TCA models
* files organized by treatment, named for mouse `td{animal number}.csv`
* output of tensortools python kit and example TCA code written by [Williams et al. 2018](https://github.com/neurostatslab/tensortools)
* columns are components for a model of rank 5. if t is number of time samples per trial, c is cell, and T is trials, first row is empty, and using Matlab indexing, rows 2:t+1 are temporal factor weights for each time, t+2:t+c+1 are neuron factor weights for each neuron, and t+c+2:t+c+T+1 are trial factor weights for each trial
* used in code sections 8-12
### 6. Threshold testing
* fake tensor of 1 simulated neuron, with identical within and across trial structure as real data. `placeboMat2.mat`
* average neuron factor weights from 100 iterations of TCA on the placebo tensor for x neurons, where x is the number of neurons for each animal, yielding a null hypothesis threshold for each animal. animals (rows) 1-14 are psilocybin, 15-21 are saline. `fMeans.csv`
* created in TCA code written by [Williams et al. 2018](https://github.com/neurostatslab/tensortools)
  
## Custom functions
Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
### crossSeshAuto
this function aligns longitudinally registered datasets
#### inputs   
* L - the cell registration file from IDPS
* d - a cell containing the cell trace files from each session in rows
* locs - properties file of last session (optional)
* nSesh - the number of sessions
#### outputs
* longRegistered - a cell of matrices of longitudinally registered cell traces from each session
* coordinates - matrix of coordinates and size of longitudinally measured cells (optional)
### aligntraces
this function aligns recorded neural activity to stimuli
#### inputs
* traces - a cell with neural traces from each session in rows
* times - a cell with stimulus times from each session in rows
* nSesh - number of sessions
* dt - sampling rate of video
* start - amount of time prior to stimulus start to include
* fin - amount of time after stimulus start to include
#### outputs 
* alignedTraces - a cell with matrices of aligned neural traces concatenated in columns
### percBaselineAvg
%this function takes finds the average trace with respect to baseline activity
#### inputs
* seshMat - cell traces for the given set of trials (time in rows, cells in columns)
* trialLength - how many seconds per trial
* timeBack - time back before stim onset to be considered baseline
* nTrials - number of trials
#### outputs
* percAvgTraces - average traces over trials in terms of % baseline
### isSig
this function checks whether cells are significantly modulated by certain stimuli
#### inputs    
* traces  - cell traces over trial of interest
* baselineTime - array of times considered baseline
* targetTime - array of times of interest
#### outputs
* changedCells - cell containing cells that are up, down, and not regulated during target time
### split_data
this function splits data into training and test sets
#### input   
* X - data
* y - classes
* test_size - % of data to be test data
#### output  
* X_test - test data
* Y_test - test classes
* X_train - training data
* Y_train - test data
### barWithError
this function makes bar graphs with SEM error bars
#### inputs   
* data - your data matrix, with variables to plot in columns and observations/replicates in rows
#### outputs
* mu - mean
* sem - standard error
* bar plot
### permTest
this function tests if means are different such that p<0.05 compared to shuffled data. 
#### inputs
* traces - data matrix, with variables to plot in columns and observations/replicates in rows
* baselineTime - array with null hypothesis start and end time
* targetTime - array with target start and end time
* numIt - number of iterations
#### outputs
* changedCells - size 3x1 cell with significantly {1} upregulated, {2} downregulated, {3} nonresponsive cell indices
### plotPersist
this function plots fractions of cells by animal
#### inputs
* dat - 3-4D data matrix of size: animals x sessions or trials x stimuli x response direction.
* groups - cell of size: number of groups x 1 size cell with animal IDs in each cell
* groupNames - cell with group names
* conditions - array or cell of size length of response direction x 1
* labels - cell of size stimuli x 1
#### outputs
* figure with subplots of animal groups in columns, stimuli in rows, and response directions in lines
### plot2metrics
this function plots bar graphs with SEM error bars over time of 2 different data sets
#### inputs
* data1 - data to be plotted in row 1, variables in columns and replicates in rows
* data2 - data to be plotted in row 2, variables in columns and replicates in rows
* titles - cell of x tick labels
* groupNames - names of groups
* groupMembers - cell of size: number of groups x 1 size cell with animal IDs in each cell
* ylab1 - ylabel for row 1
* ylab2 - ylabel for row 2
* ylim1 - y axis limits for row 1
* ylim2 - y axis limits for row 2
#### outputs
figure with subplots of animal groups in columns, data types in rows
### plotTCAmodel
this function creates TCA model plots with temporal, neural, and trial factor subplots in columns and component subplots in rows
#### inputs
* timeFactor - matrix of time factor weights for each component for each component in columns at each time in rows
* neuronFactor - matrix of neuron factor weights for each component for each component in columns at each cell in rows
* trialFactor - matrix of trial factor weights for each component for each component in columns at each trial in rows
* colorCode - array of length trials with containing the variable for desired color code
#### outputs
Figure of TCA model

## In Vivo Imaging Analysis Pipeline for Rogers et al., 2024, Longitudinal Analysis Pipeline 
Each part depends on previous parts and will not run without one another unless specified. Custom function and data folders must be downloaded and added to path. For the analysis in Fig. 3 and Supp. Fig. 2, please see the file `Rogers_et_al_2024_Figs_3_2S.m` In the following sections you will:
### 1. Load the data of all animals in a given experimental group
* Data: Behavior, dFoF, Cell pooled, session concatenated activity matrices
### 2. Extract the traces of longitudinally registered cells from all files, align to user-defined stimulus times or times of interest
* create tensors for [tensortools](https://github.com/neurostatslab/tensortools)
* downsample (from 20 to 15Hz) calcium data for freezing decoding in Code section 7
* create concatenated activity matrices in Data section 4
#### 3. Plot spatial coordinates of longitudinally registered cells
#### 4. Extract and plot average traces in terms of % baseline of all neurons  averaged over all trials in each session
#### 5. Identify cells that are upregulated or downregulated in response to stimuli in a given session according to their average trace
#### 6. Measure the number of sessions cells tend to maintain stimulus responsiveness
#### 7. Calculate freezing encoding in each neuron in each session and compare within animals across groups
#### 8. Load & plot TCA models, identifying session-dominant components and correlate to freezing
#### 9. Calculate and plot normalized trial factors.
#### 10. Calculate and plot normalized temporal factors.
#### 11. Assess cumulative distribution of neuron weights
#### 12. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and their overlaps
* if recreating primary analysis, set threshTest to 0
* if recreating supplement with custom null hypothesis neuron factors (supp fig 6) set threshTest to 1
#### 13. Use a Fisher linear decoder to predict group based on average activity (z-score relative to Acquisition) of each ensemble over time
#### 14. Extract average activity (z-score relative to Acquisition) of each cell in the ensemble during Extinction 1, Extinction 3
* if recreating primary analysis, set threshTest to 0
* if recreating supplement with custom null hypothesis neuron factors (supp fig 6) set threshTest to 1 and uncomment bottom section of code
#### 15. Plot example traces from specific ensembles.
#### 16. Plot average dFoF of non-intersecting Acq-dominant and Ext1/3-dominant ensembles and their ratios over seconds, trials

## Linear Non-linear Firing Rate Model for Rogers et al., 2024 Fig. 8
The model is described fully in the Methods section of the paper and every variable defined in code comments. This model performs the whole TFC acquisition and extinction task, looping through 8 different possible circuit architectures (bigCond == 1 being the concluded upon architecture) and 100 different "doses" of psilocybin inhibition and allows the units to learn the task over training. Outputs are the raw unit activities over the whole task, z-score with respect to acq of acq-dom units during Ext1 and ext3-dom units during Ext3, and the mean squared error of the real data with respect to these outputs.

## Data viz
For aesthetic purposes and ease of organizing statistical comparisons, most data was exported from these sections in csv or txt files in imported into Prism. In this file are all data figures and associated statistical tests

Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
