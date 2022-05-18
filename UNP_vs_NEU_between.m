function UNP_vs_NEU_between(test_mode, subjects, trials, iterations, ...
    datadir, outputdir)

%Matlab function for between-subjects power calculations of the LPP  
%contrasting unpleasant vs neutral picture categories

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

%same but for each subcategory
UH_file = [datadir, 'UH_table.csv']; %UH = unpleasant, high-arousing images
UH_table = readtable(UH_file);
UH_table = sortrows(UH_table, 1);
unique_UH = unique(UH_table.subID); 
num_trials_per_cat = length(find(unique_UH(1) == UH_table.subID));
%get number of trials per subcategory

%same for the rest of the subcategories
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

unique_subjid_UNP = unique(UNP_table.subID);
unique_subjid_NEU = unique(NEU_table.subID);
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get number of trials in neutral categoriy
total_subjects = length(unique_subjid_NEU); %get total number of subjects

% check to make sure we have the same subjects in each category
if ~isequal(unique_subjid_UNP,unique_subjid_NEU)
    disp('different sets of subjects!')
end

%get an index of trials for each category
%NEU and subcategories will be different - there are fewer trials in the 
%subcategories than their are in the parent categories
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]); 
trial_index_UNP = repmat(transpose(1:num_trials_per_cat),[total_subjects 1]); 

for i = 1:iterations %tick through all iterations
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print the loop that we're on
    
    for j = 1:subjects*2 %tick though all subjects
        %we're sampling 2x subjects to simulate 2 experimental groups
        
        random_subjects(i,:) = datasample(unique_subjid_UNP,subjects*2, ...
            'REPLACE',false); %randomly sample 2x subjects 
        %without replacement
        
        for k = 1:trials*2 %tick through the trials
            %we sample 2x NEU trials to simulate NEU vs NEU vs UH vs NEU,
            %etc
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials,...
                trials*2,'REPLACE',true); %we sample with replacement, 
            %because the subcategories only have 16 trials each and we 
            %sample up to 15 trials, so we don't have enough trials with 
            %without replacement algorithms to work
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) ...
                == trial_index_NEU & random_subjects(i,j) == ...
                NEU_table.subID); %get in index of all the trials we 
            %randomly sampled
            random_LPP_NEU(i,j,k) = ...
                NEU_table.LPP(find_trial_index_NEU(i,j,k)); %get the LPPs 
            %at that index
            for l = 1:trials %tick through all trials
                %we sample 1x trials for each subcategory (UH, etc)
                
                random_trials_UH(i,j,:) = datasample(1:num_trials_per_cat,...
                    trials,'REPLACE',true); %randomle sample trials with
                %replacement
                
                find_trial_index_UH(i,j,l) = ...
                    find(random_trials_UH(i,j,l) == trial_index_UNP ...
                    & random_subjects(i,j) == UH_table.subID); %find trial 
                %index
                random_LPP_UH(i,j,l) = ...
                    UH_table.LPP(find_trial_index_UH(i,j,l)); %get LPP 
                %at that index
                
                %same for the rest of the subcategories
                random_trials_UO(i,j,:) = ...
                    datasample(1:num_trials_per_cat,...
                    trials,'REPLACE',true);
                find_trial_index_UO(i,j,l) = ...
                    find(random_trials_UO(i,j,l) == ...
                    trial_index_UNP & random_subjects(i,j) == ...
                    UO_table.subID);
                random_LPP_UO(i,j,l) = ...
                    UO_table.LPP(find_trial_index_UO(i,j,l));
                
                random_trials_UL(i,j,:) = ...
                    datasample(1:num_trials_per_cat,...
                    trials,'REPLACE',true);
                find_trial_index_UL(i,j,l) = ...
                    find(random_trials_UL(i,j,l) == ...
                    trial_index_UNP & random_subjects(i,j) ...
                    == UL_table.subID);
                random_LPP_UL(i,j,l) = ...
                    UL_table.LPP(find_trial_index_UL(i,j,l));
            end   
        end
    end
    %assign LPPs to groups/conditions
    random_LPP_UH_group1(i,:,:) = random_LPP_UH(i,1:subjects,:);
    random_LPP_UO_group1(i,:,:) = random_LPP_UO(i,1:subjects,:);
    random_LPP_UL_group1(i,:,:) = random_LPP_UL(i,1:subjects,:);
    random_LPP_NEU_group1(i,:,:) = random_LPP_NEU(i,1:subjects,...
        1:trials);
    random_LPP_NEU_group2a(i,:,:) = random_LPP_NEU(i, ...
        (subjects+1):(subjects*2),1:trials);
    random_LPP_NEU_group2b(i,:,:) = random_LPP_NEU(i, ...
        (subjects+1):(subjects*2), ...
        (trials+1):(trials*2));
end

% check to ensure group/condition sizes are consistent
if size(random_LPP_UH_group1) ~= size(random_LPP_NEU_group2b)
    disp('groups are different sizes!')
end

%get dimensions of groups by condition
[a, b, ~] = size(random_LPP_UH_group1);

for aa = 1:a %tick through dimension 1 (iterations)
    
    for ab = 1:b %tick through dimension 2 (subjects)
        
        %get average LPP accross all trials of each condition for each
        %subject, each iteration
        mean_within_subject_UH(aa,ab,:) = ...
            mean(random_LPP_UH_group1(aa,ab,:));
        mean_within_subject_UO(aa,ab,:) = ...
            mean(random_LPP_UO_group1(aa,ab,:));
        mean_within_subject_UL(aa,ab,:) = ...
            mean(random_LPP_UL_group1(aa,ab,:));
        mean_within_subject_NEU1(aa,ab,:) = ...
            mean(random_LPP_NEU_group1(aa,ab,:));
        mean_within_subject_NEU2a(aa,ab,:) = ...
            mean(random_LPP_NEU_group2a(aa,ab,:));
        mean_within_subject_NEU2b(aa,ab,:) = ...
            mean(random_LPP_NEU_group2b(aa,ab,:));
    end
end

%get the contrasts for the conditions
UH_minus_NEU = mean_within_subject_UH - mean_within_subject_NEU1;
UO_minus_NEU = mean_within_subject_UO - mean_within_subject_NEU1;
UL_minus_NEU = mean_within_subject_UL - mean_within_subject_NEU1;

NEU_minus_NEU = mean_within_subject_NEU2a - mean_within_subject_NEU2b;

%reshape the data into a format that is easier to read
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects*2,1]);

reshape_UH_minus_NEU = reshape(UH_minus_NEU', ...
    [iterations*subjects, 1]);
reshape_UO_minus_NEU = reshape(UO_minus_NEU', ...
    [iterations*subjects, 1]);
reshape_UL_minus_NEU = reshape(UL_minus_NEU', ...
    [iterations*subjects, 1]);

reshape_NEU_minus_NEU = reshape(NEU_minus_NEU', ...
    [iterations*subjects, 1]);

%create an array for the group assignment
grouponearray = ones([iterations*subjects, 1]);
grouptwoarray = repmat(2, [iterations*subjects, 1]);
group_array = [grouponearray; grouptwoarray];

%stack the groups on top of eachother
UH_vs_NEU_array = [reshape_UH_minus_NEU; reshape_NEU_minus_NEU];
UO_vs_NEU_array = [reshape_UO_minus_NEU; reshape_NEU_minus_NEU];
UL_vs_NEU_array = [reshape_UL_minus_NEU; reshape_NEU_minus_NEU];

%create a test directory if necessary
if test_mode == 1
    outputdir = 'C:\path_to_output_dir\';
    mkdir(outputdir)
else
    mkdir(outputdir)
end

%create an output directory
between_results_dir = [outputdir, 'real_results_between\'];
mkdir(between_results_dir)

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

%write to table format
UH_vs_NEU_between_table = table(reshape_random_subjects, ...
    UH_vs_NEU_array, group_array,'VariableNames', ...
    {'subjectID','LPP','group'});
%assign a filename
UH_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'UHvsNEU_between_', ...
    iterations_string, '_iterations.csv'];
%write table to csv
writetable(UH_vs_NEU_between_table, UH_vs_NEU_between_filename)

%same for the rest of the subcategories
UO_vs_NEU_between_table = table(reshape_random_subjects, ...
    UO_vs_NEU_array, group_array,'VariableNames', ...
    {'subjectID','LPP','group'});
UO_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'UOvsNEU_between_', ...
    iterations_string, '_iterations.csv'];
writetable(UO_vs_NEU_between_table, UO_vs_NEU_between_filename)


UL_vs_NEU_between_table = table(reshape_random_subjects, ...
    UL_vs_NEU_array, group_array,'VariableNames', ...
    {'subjectID','LPP','group'});
UL_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'ULvsNEU_between_', ...
    iterations_string, '_iterations.csv'];
writetable(UL_vs_NEU_between_table, UL_vs_NEU_between_filename)
end