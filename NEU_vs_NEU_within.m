function NEU_vs_NEU_within(test_mode, subjects, trials, iterations, datadir, ...
    outputdir)

%Matlab function for within-subjects power calculations of the LPP  
%contrasting neutral vs neutral picture categories

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


%get NEU (neutral images) data
NEU_file = [datadir, 'NEU_table.csv'];
NEU_table = readtable(NEU_file);
NEU_table = sortrows(NEU_table, 1); %sort rows by subject ID


unique_subjid_NEU = unique(NEU_table.subID); %get all subjects
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get total number of trials
total_subjects = length(unique_subjid_NEU); %get total number of subjects

%get an index of all the neutral trials
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]);

for i = 1:iterations %tick through each iteration
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects %tick through each subject
        random_subjects(i,:) = datasample(unique_subjid_NEU,subjects,...
            'REPLACE',false); %randomly sample subjects without replacement
        
        for k = 1:trials*2 %tick through each trial
            %we need 2x trials to simulate 2 experimental conditions
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials,...
                trials*2,'REPLACE',true); %randomly sample 2x trials with 
                %replacement we don't have enough trials in the 
                %subcategories to sample without replacement
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) ...
                == trial_index_NEU & random_subjects(i,j) == ...
                NEU_table.subID);
                %find the index of the trials we randomly sampled
            
            random_LPP_NEU(i,j,k) = ...
                NEU_table.LPP(find_trial_index_NEU(i,j,k));
                %find the LPP at that index
        end
        
        %assign LPPs to conditions
        random_LPP_NEU_condition1(i,j,:) = random_LPP_NEU(i,j,1:trials);
        random_LPP_NEU_condition2(i,j,:) = random_LPP_NEU(i,j,(trials+1):(trials*2));
        
        %get the average for each iteration, each subject, each condition 
        %accross all the sampled trials
        mean_within_NEU1(i,j) = mean(random_LPP_NEU_condition1(i,j,:));
        mean_within_NEU2(i,j) = mean(random_LPP_NEU_condition2(i,j,:));
    end
end
   

%%
%now reshape the files so they can write to a spreadsheet
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects,1]);
reshape_LPP_NEU_1 = reshape(mean_within_NEU1', [iterations*subjects,1]);
reshape_LPP_NEU_2 = reshape(mean_within_NEU2', [iterations*subjects,1]);

%check that our dimensions are consistent
if size(reshape_LPP_NEU_1) ~= size(reshape_LPP_NEU_2)
    disp('group 1 and 2 are different sizes!')
end

%get total # of data points
num_observations = iterations*subjects*2;

index1 = 1:2:num_observations;%create an index going from 1 to the end of 
%the data in steps of 2
index2 = 2:2:num_observations;%create an index going from 2 to the end of 
%the data in steps of 2

for q = 1:length(index1) %tick through each iteration
    
    subjectsvec(index1(q)) = reshape_random_subjects(q); %get subjects at 
    %each iteration
    subjectsvec(index2(q)) = reshape_random_subjects(q);
    LPP_NEUvsNEU(index1(q)) = reshape_LPP_NEU_1(q); %get LPPs
    LPP_NEUvsNEU(index2(q)) = reshape_LPP_NEU_2(q);
    condition(index1(q)) = 1; %get condition indicators
    condition(index2(q)) = 2;
end

%transpose into column vectors
subjectsvec = subjectsvec';
LPP_NEUvsNEU = LPP_NEUvsNEU';
condition = condition';

%create a test directory if necessary
if test_mode == 1
    outputdir = 'C:\path_to_output_dir\';
    mkdir(outputdir)
else
    mkdir(outputdir)
end

%make an output directory
within_results_dir = [outputdir, 'real_results_within\'];
mkdir(within_results_dir)


%create strings for the filename
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
NEUvsNEU_table = table(subjectsvec,LPP_NEUvsNEU,condition,...
    'VariableNames', {'subjectID','LPP','condition'}); 
    %condition 1 is NEU, condition 2 is NEU

%assign a filename
NEUvsNEU_filename = [within_results_dir, subjects_string, '_subjects_', ...
    trials_string, '_trials_', 'NEUvsNEU_', iterations_string, ...
    '_iterations.csv'];

%write table to csv
writetable(NEUvsNEU_table, NEUvsNEU_filename);
end