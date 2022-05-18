function UNP_vs_NEU_within(test_mode, subjects, trials, iterations, datadir, ...
    outputdir)

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

if ~exist('test_mode','var')%if you don't pass in a test_mode variable,
    %it defaults to zero
    test_mode = 0;
end


if ~exist('subjects','var')%if you don't pass in a subjects variable,
    %it defaults to zero
    subjects = 10;
end

if ~exist('trials','var')%if you don't pass in a trials variable,
    %it defaults to ten
    trials = 10;
end

if ~exist('iterations','var')%if you don't pass in an iterations 
    %variable, it defaults to one hundred
    iterations = 100;
end

if ~exist('datadir','var')%if you didn't specify a datadir,
    %use the path to the default data dir here
    datadir = 'C:\path_to_data_dir\';
end

if ~exist('outputdir','var')%if you didn't specify a datadir,
    %use the path to the default output dir here
    outputdir = 'C:\path_to_output_dir\';
end

%%
%BEGIN FUNCTION

%data for each condition stored in a seperate .csv file on the datadir


%get UNP (unpleasant images) data
UNP_file = [datadir, 'UNP_table.csv'];
UNP_table = readtable(UNP_file);
UNP_table = sortrows(UNP_table,1); %sort rows by subject ID

%same for the subcategories
UH_file = [datadir, 'UH_table.csv']; %UH = unpleasant, high-arousing images
UH_table = readtable(UH_file);
UH_table = sortrows(UH_table, 1);
unique_UH = unique(UH_table.subID); %get the subject IDs
num_trials_per_cat = length(find(unique_UH(1) == UH_table.subID)); %get the
%number of trials per subcategory

UO_file = [datadir, 'UO_table.csv']; %UO = unpleasant objects
UO_table = readtable(UO_file);
UO_table = sortrows(UO_table, 1);

UL_file = [datadir, 'UL_table.csv']; %UL = unpleasant, low-arousing
UL_table = readtable(UL_file);
UL_table = sortrows(UL_table, 1);

%get NEU (neutral images) data
NEU_file = [datadir, 'NEU_table.csv'];
NEU_table = readtable(NEU_file);
NEU_table = sortrows(NEU_table, 1);


unique_subjid_UNP = unique(UNP_table.subID); %get the subIDs in the NEU and
%UNP categories
unique_subjid_NEU = unique(NEU_table.subID);
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get the number of trials in the NEU category

total_subjects = length(unique_subjid_NEU); %get total number of subjects

%verify that we have the same nubmer of subjects in each category
if ~isequal(unique_subjid_UNP,unique_subjid_NEU)
    disp('different sets of subjects!')
end

%create an index of all the trials in each category
%these will be different sizes because we are using subcategories for UNP
%and the parent category for NEU
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]); 
trial_index_UNP = repmat(transpose(1:num_trials_per_cat),[total_subjects 1]); 

for i = 1:iterations %tick though each iteration
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects %tick though each subject
        
        random_subjects(i,:) = datasample(unique_subjid_UNP,subjects, ...
            'REPLACE',false); %randomly sample subjects without 
        %replacement
        
        for k = 1:trials %tick through each trial
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials,trials, ...
                'REPLACE',true); %randomly sample neutral trials with 
            %replacement becauase we don't have enough trials in the
            %subcategories to sample without
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) == ...
                trial_index_NEU & random_subjects(i,j) == NEU_table.subID);
            %get the index of all the trials we sampled
            
            random_LPP_NEU(i,j,k) = NEU_table.LPP(find_trial_index_NEU(i,j,k));
            %find the LPP that matches that index
            
            %same thing but with the subcategories
            random_trials_UH(i,j,:) = datasample(1:num_trials_per_cat,trials, ...
                'REPLACE', true);
            find_trial_index_UH(i,j,k) = find(random_trials_UH(i,j,k) == ...
                trial_index_UNP & random_subjects(i,j) == UH_table.subID);
            random_LPP_UH(i,j,k) = UH_table.LPP(find_trial_index_UH(i,j,k));
            
            random_trials_UO(i,j,:) = datasample(1:num_trials_per_cat,trials, ...
                'REPLACE',true);
            find_trial_index_UO(i,j,k) = find(random_trials_UO(i,j,k) == ...
                trial_index_UNP & random_subjects(i,j) == UO_table.subID);
            random_LPP_UO(i,j,k) = UO_table.LPP(find_trial_index_UO(i,j,k));
            
            random_trials_UL(i,j,:) = datasample(1:num_trials_per_cat,trials, ...
                'REPLACE',true);
            find_trial_index_UL(i,j,k) = find(random_trials_UL(i,j,k) == ...
                trial_index_UNP & random_subjects(i,j) == UL_table.subID);
            random_LPP_UL(i,j,k) = UL_table.LPP(find_trial_index_UL(i,j,k));
        end
        %get the average LPP accross all trials for each iteration, each
        %subject, each category
        mean_within_UH(i,j,:) = mean(random_LPP_UH(i,j,:));
        mean_within_UO(i,j,:) = mean(random_LPP_UO(i,j,:));
        mean_within_UL(i,j,:) = mean(random_LPP_UL(i,j,:));
        mean_within_NEU(i,j,:) = mean(random_LPP_NEU(i,j,:));
    end
end
     

%%
%now reshape the matrices so they can write to a table
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects,1]);

reshape_LPP_UH = reshape(mean_within_UH', [iterations*subjects,1]);
reshape_LPP_UO = reshape(mean_within_UO', [iterations*subjects,1]);
reshape_LPP_UL = reshape(mean_within_UL', [iterations*subjects,1]);
reshape_LPP_NEU = reshape(mean_within_NEU', [iterations*subjects,1]);

%verfiy that the dimensions match
if length(reshape_LPP_UH) ~= length(reshape_LPP_NEU)
    disp('number of observations per category inconsistent!')
end

if length(reshape_LPP_UH) ~= length(reshape_random_subjects)
    disp('# subjects and # observations inconsistent!')
end

%get total number of observations
num_observations = iterations*subjects*2;

index1 = 1:2:num_observations; %create an index going from 1 to the end of 
%the data in steps of 2
index2 = 2:2:num_observations; %create an index going from 2 to the end of 
%the data in steps of 2

for q = 1:length(index1) %tick through each iteration
    
    subjectsvec(index1(q)) = reshape_random_subjects(q); %get subjects at 
    %each iteration
    subjectsvec(index2(q)) = reshape_random_subjects(q);
    LPP_UHvsNEU(index1(q)) = reshape_LPP_UH(q); %get LPP for each category,
    %each iteration
    LPP_UHvsNEU(index2(q)) = reshape_LPP_NEU(q);
    LPP_UOvsNEU(index1(q)) = reshape_LPP_UO(q);
    LPP_UOvsNEU(index2(q)) = reshape_LPP_NEU(q);
    LPP_ULvsNEU(index1(q)) = reshape_LPP_UL(q);
    LPP_ULvsNEU(index2(q)) = reshape_LPP_NEU(q);
    
    condition(index1(q)) = 1; %get condition assignment for each iteration
    condition(index2(q)) = 2;
end

%transpose into column vectors
subjectsvec = subjectsvec';
LPP_UHvsNEU = LPP_UHvsNEU';
LPP_UOvsNEU = LPP_UOvsNEU';
LPP_ULvsNEU = LPP_ULvsNEU';
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

%make strings for the file name
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
UHvsNEU_table = table(subjectsvec,LPP_UHvsNEU,condition,...,
    'VariableNames', {'subjectID','LPP','condition'}); %condition 1 is UH, 
%condition 2 is NEU

%assign a file name
UHvsNEU_filename = [within_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'UHvsNEU_', ...
    iterations_string, '_iterations.csv'];

%write table to csv
writetable(UHvsNEU_table, UHvsNEU_filename);

%same thing for the rest of the subcategories
UOvsNEU_table = table(subjectsvec,LPP_UOvsNEU,condition,...
    'VariableNames', {'subjectID','LPP','condition'}); %condition 1 is UO, 
%condition 2 is NEU
UOvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'UOvsNEU_', iterations_string, ...
    '_iterations.csv'];
writetable(UOvsNEU_table, UOvsNEU_filename);

ULvsNEU_table = table(subjectsvec,LPP_ULvsNEU,condition,...
    'VariableNames', {'subjectID','LPP','condition'}); %condition 1 is UL, 
%condition 2 is NEU
ULvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'ULvsNEU_', iterations_string, ...
    '_iterations.csv'];
writetable(ULvsNEU_table, ULvsNEU_filename);
end