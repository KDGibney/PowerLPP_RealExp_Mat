function PLE_vs_NEU_within(test_mode, subjects, trials, iterations, ...
    datadir, outputdir)

%Matlab function for within-subjects power calculations of the LPP  
%contrasting pleasant vs neutral picture categories

%authored by Kyla Gibney: kdgibney.work@gmail.com

%there are some redundancies with this function and the caller function 
%call_by_contrast because sometimes I want to test this function by itself, 
%without using the caller funciton

%%
%INS

%test_mode: option to enable test mode. 1 = test mode, 0 = not.

%subjects: number of subjects to be sampled

%trials: number of trials to be sampled

%iterations: number of times to repeat the simulation

%datadir: directory where the data is stored

%outputdir: directory where results files will be written

%%
%OUTS

%this function does not return any variables to the workspace, but does
%write results files to outputdir

%%
%DEFAULTS

%if you don't pass in any input areguments to the function, the function
%will default to these

if ~exist('test_mode','var')
    test_mode = 0; %if you don't pass in a test_mode variable,
    %it defaults to zero
end


if ~exist('subjects','var')
    subjects = 10; %if you don't pass in a subjects variable,
    %it defaults to zero
end

if ~exist('trials','var') %if you don't pass in a trials variable,
    %it defaults to ten
    trials = 10;
end

if ~exist('iterations','var') %if you don't pass in an iterations 
    %variable, it defaults to one hundred
    iterations = 100;
end

if ~exist('datadir','var') %if you didn't specify a datadir,
    %use the path to the default data dir here
    datadir = 'C:\path_to_data_dir\';
end

if ~exist('outputdir','var') %if you didn't specify a datadir,
    %use the path to the default output dir here
    outputdir = 'C:\path_to_output_dir\';
end

%%
%BEGIN FUNCTION

%data for each condition stored in a seperate .csv file on the datadir

%get PLE (pleasant images) data
PLE_file = [datadir, 'PLE_table.csv'];
PLE_table = readtable(PLE_file);
PLE_table = sortrows(PLE_table,1); %sort rows by subject ID

%same for each subcategory
PH_file = [datadir, 'PH_table.csv']; %PH = pleasant, high arousing
PH_table = readtable(PH_file);
PH_table = sortrows(PH_table, 1);
unique_PH = unique(PH_table.subID); %get all subjecIDs
num_trials_per_cat = length(find(unique_PH(1) == PH_table.subID)); %get 
%the number of trials per subcategory

%same for the remining subcategories
PO_file = [datadir, 'PO_table.csv']; %PO = pleasant objects
PO_table = readtable(PO_file);
PO_table = sortrows(PO_table, 1);

PL_file = [datadir, 'PL_table.csv']; %PL = pleasant, low arousing
PL_table = readtable(PL_file);
PL_table = sortrows(PL_table, 1);

%get NEU (neutral images) data
NEU_file = [datadir, 'NEU_table.csv'];
NEU_table = readtable(NEU_file);
NEU_table = sortrows(NEU_table, 1);


unique_subjid_PLE = unique(PLE_table.subID);
unique_subjid_NEU = unique(NEU_table.subID);
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get total number of NEU trials

total_subjects = length(unique_subjid_NEU); %get total # subjects

%verify that we have the same subjects in each category
if ~isequal(unique_subjid_PLE,unique_subjid_NEU)
    disp('different sets of subjects!')
end

%create indices of neutral and pleasant trials
%these will be diffetent sizes because NEU is a parent category and in PLE
%we are using subcategories
trial_index_NEU = repmat(transpose(1:num_NEU_trials),...
    [total_subjects 1]); 
trial_index_PLE = repmat(transpose(1:num_trials_per_cat),...
    [total_subjects 1]); 

for i = 1:iterations %tick through each iteration
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects %tick through each subject
        
        random_subjects(i,:) = datasample(unique_subjid_PLE,subjects, ...
            'REPLACE',false); %randomly sample subjects without 
        %replacement
        
        for k = 1:trials %tick though each trial
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials, ...
                trials, 'REPLACE',true); %randomly sample trials 
            %with replacement because we don't have enough trials in the 
            %subcategories to sample without replacement
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) ...
                == trial_index_NEU & random_subjects(i,j) ...
                == NEU_table.subID); %find the index of the trials 
            %we just sampled
            
            random_LPP_NEU(i,j,k) = ...
                NEU_table.LPP(find_trial_index_NEU(i,j,k));
            %find the LPP that matches that index
            
            %same but for the pleasant subcategories
            random_trials_PH(i,j,:) = datasample(1:num_trials_per_cat,...
                trials, 'REPLACE',true); 
            find_trial_index_PH(i,j,k) = find(random_trials_PH(i,j,k) ...
                == trial_index_PLE & random_subjects(i,j) == ...
                PH_table.subID);
            random_LPP_PH(i,j,k) = ...
                PH_table.LPP(find_trial_index_PH(i,j,k));
            
            random_trials_PO(i,j,:) = datasample(1:num_trials_per_cat,...
                trials, 'REPLACE',true);
            find_trial_index_PO(i,j,k) = find(random_trials_PO(i,j,k) ...
                == trial_index_PLE & random_subjects(i,j) ...
                == PO_table.subID);
            random_LPP_PO(i,j,k) = ...
                PO_table.LPP(find_trial_index_PO(i,j,k));
            
            random_trials_PL(i,j,:) = datasample(1:num_trials_per_cat,...
                trials, 'REPLACE', true);
            find_trial_index_PL(i,j,k) = find(random_trials_PL(i,j,k) ...
                == trial_index_PLE & random_subjects(i,j) ...
                == PL_table.subID);
            random_LPP_PL(i,j,k) = ...
                PL_table.LPP(find_trial_index_PL(i,j,k));
        end
        
        %get the mean LPP accross the sampled trials for each iteration, 
        %each subject, each category
        mean_within_PH(i,j) = mean(random_LPP_PH(i,j,:));
        mean_within_PO(i,j) = mean(random_LPP_PO(i,j,:));
        mean_within_PL(i,j) = mean(random_LPP_PL(i,j,:));
        mean_within_NEU(i,j) = mean(random_LPP_NEU(i,j,:));
    end
end


%%
%now reshape the files so they can write to a table
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects,1]);
reshape_LPP_PH = reshape(mean_within_PH', [iterations*subjects,1]);
reshape_LPP_PO = reshape(mean_within_PO', [iterations*subjects,1]);
reshape_LPP_PL = reshape(mean_within_PL', [iterations*subjects,1]);
reshape_LPP_NEU = reshape(mean_within_NEU', [iterations*subjects,1]);

%verify that the dimensions match
if length(reshape_LPP_PH) ~= length(reshape_LPP_NEU)
    disp('number of observations per category inconsistent!')
end

if length(reshape_LPP_PH) ~= length(reshape_random_subjects)
    disp('# subjects and # observations inconsistent!')
end

%get the total # of data points
num_observations = iterations*subjects*2;

index1 = 1:2:num_observations;%create an index going from 1 to the end of 
%the data in steps of 2
index2 = 2:2:num_observations;%create an index going from 2 to the end of 
%the data in steps of 2

for q = 1:length(index1) %tick through each iteration
    
    subjectsvec(index1(q)) = reshape_random_subjects(q); %get subjects at 
    %each iteration
    subjectsvec(index2(q)) = reshape_random_subjects(q);
    
    LPP_PHvsNEU(index1(q)) = reshape_LPP_PH(q); %get LPP for each category
    %at each iteration
    LPP_PHvsNEU(index2(q)) = reshape_LPP_NEU(q);
    
    LPP_POvsNEU(index1(q)) = reshape_LPP_PO(q);
    LPP_POvsNEU(index2(q)) = reshape_LPP_NEU(q);
    
    LPP_PLvsNEU(index1(q)) = reshape_LPP_PL(q);
    LPP_PLvsNEU(index2(q)) = reshape_LPP_NEU(q);
    
    condition(index1(q)) = 1; %assign the conditions for each iteration
    condition(index2(q)) = 2;
end

%transpose into column vectors
subjectsvec = subjectsvec';
LPP_PHvsNEU = LPP_PHvsNEU';
LPP_POvsNEU = LPP_POvsNEU';
LPP_PLvsNEU = LPP_PLvsNEU';
condition = condition';


%create a test directory if necessary
if test_mode == 1
    outputdir = ['C:\path_to_output_dir\'];
    mkdir(outputdir)
else
    mkdir(outputdir)
end

%make an output directory
within_results_dir = [outputdir, 'real_results_within\'];
mkdir(within_results_dir)


%create strings for the file name
iterations_string = num2str(iterations);
trials_string = num2str(trials);

if trials < 10
    trials_string = ['0', trials_string];
end

subjects_string = num2str(subjects);

if subjects < 100
    subjects_string = ['0', subjects_string];
end

%write results to a table
PHvsNEU_table = table(subjectsvec,LPP_PHvsNEU,condition,...
    'VariableNames', {'subjectID','LPP','condition'});
%condition 1 is PH, condition 2 is NEU

%assign a file name
PHvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'PHvsNEU_', iterations_string, ...
    '_iterations.csv'];

%write the table to csv
writetable(PHvsNEU_table, PHvsNEU_filename);

%same for the rest of the subcategories
POvsNEU_table = table(subjectsvec,LPP_POvsNEU,condition,...
    'VariableNames', {'subjectID','LPP','condition'}); %condition 1 is PO, 
%condition 2 is NEU
POvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'POvsNEU_', iterations_string, ...
    '_iterations.csv'];
writetable(POvsNEU_table, POvsNEU_filename);

PLvsNEU_table = table(subjectsvec,LPP_PLvsNEU,condition, ...
    'VariableNames', {'subjectID','LPP','condition'}); %condition 1 is PL, 
%condition 2 is NEU
PLvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'PLvsNEU_', iterations_string, ...
    '_iterations.csv'];
writetable(PLvsNEU_table, PLvsNEU_filename);
end