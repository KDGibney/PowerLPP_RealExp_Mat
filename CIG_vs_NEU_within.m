function CIG_vs_NEU_within(test_mode, subjects, trials, iterations, datadir, ...
    outputdir)

%Matlab function for within-subjects power calculations of the LPP  
%contrasting cigarette vs neutral picture categories

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
    %it defaults to ten
    subjects = 10;
end

if ~exist('trials','var')%if you don't pass in a trials variable,
    %it defaults to ten
    trials = 10;
end

if ~exist('iterations','var')%if you don't pass in an iterations variable,
    %it defaults to one hundred
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


%get CIG (cigarette images) data
CIG_file = [datadir, 'CIG_table.csv'];
CIG_table = readtable(CIG_file);
CIG_table = sortrows(CIG_table,1); %sort rows by subject ID


%same for NEU (neutral images) data
NEU_file = [datadir, 'NEU_table.csv'];
NEU_table = readtable(NEU_file);
NEU_table = sortrows(NEU_table, 1);

unique_subjid_CIG = unique(CIG_table.subID); %get subject IDs
unique_subjid_NEU = unique(NEU_table.subID);

num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get number of NEU trials
num_CIG_trials = length(find(unique_subjid_CIG(1) == CIG_table.subID));
%number of CIG trials
total_subjects = length(unique_subjid_NEU); %get total # subjects

%we should have the same subjects in each category
if ~isequal(unique_subjid_CIG,unique_subjid_NEU)
    disp('different sets of subjects!')
end

%create an index of trials for each category
%these should be the same size because they are both parent categories
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]); 
trial_index_CIG = repmat(transpose(1:num_CIG_trials),[total_subjects 1]); 

for i = 1:iterations %tick through each iteration
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects %tick through each subject
        %we sample 1x subjects for a within-subjects design
        
        random_subjects(i,:) = datasample(unique_subjid_CIG,subjects,...
            'REPLACE',false); %randomly sample subjects without replacement
        
        for k = 1:trials %tick throug each trial
            %we need 1x the number of trials for each category
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials,...
                trials,'REPLACE',true); %randomly sample neutral trials 
            %with replacement, because we don't have enough trials in the
            %subcategories for sampling without replacement
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) == ...
                trial_index_NEU & random_subjects(i,j) == NEU_table.subID);
            %find index of neutral trials
            
            random_LPP_NEU(i,j,k) = NEU_table.LPP(find_trial_index_NEU(i,j,k));
            %get LPP that matches that index
            
            %same as before but with CIG category
            random_trials_CIG(i,j,:) = datasample(1:num_CIG_trials,trials, ...
                'REPLACE',true);
            find_trial_index_CIG(i,j,k) = find(random_trials_CIG(i,j,k) == ...
                trial_index_CIG & random_subjects(i,j) == CIG_table.subID);
            random_LPP_CIG(i,j,k) = CIG_table.LPP(find_trial_index_CIG(i,j,k));          
        end
        
        %get the average accross trials for each category, each subject,
        %each iteration
        mean_within_CIG(i,j,:) = mean(random_LPP_CIG(i,j,:));
        mean_within_NEU(i,j,:) = mean(random_LPP_NEU(i,j,:));
    end
end


%%
%now reshape the files so they can write to a spreadsheet
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects,1]);

reshape_LPP_CIG = reshape(mean_within_CIG', [iterations*subjects,1]);

reshape_LPP_NEU = reshape(mean_within_NEU', [iterations*subjects,1]);

%check to make sure our dimension are consistent
if length(reshape_LPP_CIG) ~= length(reshape_LPP_NEU)
    disp('number of observations per category inconsistent!')
end

if length(reshape_LPP_CIG) ~= length(reshape_random_subjects)
    disp('# subjects and # observations inconsistent!')
end

%get total # data points
num_observations = iterations*subjects*2;

index1 = 1:2:num_observations; %create an index going from 1 to the end of 
%the data in steps of 2
index2 = 2:2:num_observations; %create an index starting at 2 to the end of
%the data in steps of 2

for q = 1:length(index1) %tick through each iteration
    
    subjectsvec(index1(q)) = reshape_random_subjects(q); %get subjects at 
    %each iteration
    subjectsvec(index2(q)) = reshape_random_subjects(q);
    
    LPP_CIGvsNEU(index1(q)) = reshape_LPP_CIG(q); %get CIG LPPs
    LPP_CIGvsNEU(index2(q)) = reshape_LPP_NEU(q); %get NEU LPPs
    
    condition(index1(q)) = 1; %get condition markers
    condition(index2(q)) = 2;
end

%transpose into column vectors
subjectsvec = subjectsvec';
LPP_CIGvsNEU = LPP_CIGvsNEU';
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

%create strings for the file names
iterations_string = num2str(iterations);
trials_string = num2str(trials);

if trials < 10
    trials_string = ['0', trials_string];
end

subjects_string = num2str(subjects);
if subjects < 100
    subjects_string = ['0', subjects_string];
end

%write results into a table
CIGvsNEU_table = table(subjectsvec,LPP_CIGvsNEU,condition,'VariableNames', ...
    {'subjectID','LPP','condition'}); %condition 1 is CIG, condition 2 is NEU

%assign a filename
CIGvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'CIGvsNEU_', iterations_string, ...
    '_iterations.csv'];

%write table to csv
writetable(CIGvsNEU_table, CIGvsNEU_filename);
end