function PLE_vs_NEU_between(test_mode, subjects, trials, iterations, ...
    datadir, outputdir)

%Matlab function for between-subjects power calculations of the LPP  
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


%get PLE (pleasant images) data
PLE_file = [datadir, 'PLE_table.csv'];
PLE_table = readtable(PLE_file);
PLE_table = sortrows(PLE_table,1); %sort rows by subject ID

%get PLE data by subcategory
PH_file = [datadir, 'PH_table.csv']; %PH = pleasant, high arousing
PH_table = readtable(PH_file);
PH_table = sortrows(PH_table, 1);
unique_PH = unique(PH_table.subID);
num_trials_per_cat = length(find(unique_PH(1) == PH_table.subID));
%get number of trials per subcategory

%same for PO 
PO_file = [datadir, 'PO_table.csv']; %PO = pleasant objects
PO_table = readtable(PO_file);
PO_table = sortrows(PO_table, 1);

%sample for PL
PL_file = [datadir, 'PL_table.csv']; %PL = pleasant, low arousing
PL_table = readtable(PL_file);
PL_table = sortrows(PL_table, 1);

%get NEU (neutral images) data
NEU_file = [datadir, 'NEU_table.csv'];
NEU_table = readtable(NEU_file);
NEU_table = sortrows(NEU_table, 1); %sort rows by subject ID


unique_subjid_PLE = unique(PLE_table.subID);
unique_subjid_NEU = unique(NEU_table.subID);
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get number of NEU trials

total_subjects = length(unique_subjid_NEU); %get total # of subjects

% check to make sure subjects are consistent
if ~isequal(unique_subjid_PLE,unique_subjid_NEU)
    disp('different sets of subjects!')
end

%make indices of the trials
%these will be different between PLE and NEU becuase PLE is by
%subcategories
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]); 
trial_index_PLE = repmat(transpose(1:num_trials_per_cat),[total_subjects 1]); 

for i = 1:iterations %tick through all iterations
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects*2 %tick through all subjects
        %we're sampling 2x subjects to simulate 2 experimental groups
        
        random_subjects(i,:) = datasample(unique_subjid_PLE,subjects*2, ...
            'REPLACE',false); %randomly sample subjects without replacement
        
        for k = 1:trials*2 %tick through all trials
            %we're sampling 2x trials to simulate to experimental
            %conditions
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials,trials*2,...
                'REPLACE',true); %randomly sample NEU trials with replacement
            %we sample with replacement for the trials in the real
            %experiments because the subcategories only go up to 16 trials,
            %and we sample up to 15 trials, so we don't have enough trials
            %to sample without replacement
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) == ...
                trial_index_NEU & random_subjects(i,j) == NEU_table.subID);
            %find the index of the neutral trials
            
            random_LPP_NEU(i,j,k) = NEU_table.LPP(find_trial_index_NEU(i,j,k));
            %get the LPPs that match those indices
            
            for l = 1:trials %tick through all trials
                %we're sampling 1x the number of trials from the
                %subcategories, to simulate NEU vs NEU vs PH vs NEU, etc
                
                random_trials_PH(i,j,:) = datasample(1:num_trials_per_cat, ...
                    trials,'REPLACE',true); %sample PH trials with 
                %replacement
                
                find_trial_index_PH(i,j,l) = find(random_trials_PH(i,j,l) == ...
                    trial_index_PLE & random_subjects(i,j) == PH_table.subID);
                %get the index of the PH trials
                
                random_LPP_PH(i,j,l) = PH_table.LPP(find_trial_index_PH(i,j,l));
                %get the LPPs that match that index
                
                %do the same for PO
                random_trials_PO(i,j,:) = datasample(1:num_trials_per_cat,...
                    trials,'REPLACE',true);
                find_trial_index_PO(i,j,l) = find(random_trials_PO(i,j,l) == ...
                    trial_index_PLE & random_subjects(i,j) == PO_table.subID);
                random_LPP_PO(i,j,l) = PO_table.LPP(find_trial_index_PO(i,j,l));
                
                %same thing again but for PL
                random_trials_PL(i,j,:) = datasample(1:num_trials_per_cat,...
                    trials,'REPLACE',true);
                find_trial_index_PL(i,j,l) = find(random_trials_PL(i,j,l) == ...
                    trial_index_PLE & random_subjects(i,j) == PL_table.subID);
                random_LPP_PL(i,j,l) = PL_table.LPP(find_trial_index_PL(i,j,l));
                
            end
                 
        end
              
    end
    
    %assign LPPs to groups and conditions
    random_LPP_PH_group1(i,:,:) = random_LPP_PH(i,1:subjects,:);
    random_LPP_PO_group1(i,:,:) = random_LPP_PO(i,1:subjects,:);
    random_LPP_PL_group1(i,:,:) = random_LPP_PL(i,1:subjects,:);
    
    random_LPP_NEU_group1(i,:,:) = random_LPP_NEU(i,1:subjects,...
        1:trials);
    random_LPP_NEU_group2a(i,:,:) = random_LPP_NEU(i, ...
        (subjects+1):(subjects*2),1:trials);
    random_LPP_NEU_group2b(i,:,:) = random_LPP_NEU(i, ...
        (subjects+1):(subjects*2), ...
        (trials+1):(trials*2));
end

% check to make sure groups/conditions are the same size
if size(random_LPP_PH_group1) ~= size(random_LPP_NEU_group2b)
    disp('groups are different sizes!')
end

%get the dimensions of the groups/conditions
[a, b , ~] = size(random_LPP_PH_group1);

for aa = 1:a %tick through dimension 1 (iterations)
    for ab = 1:b %tick though dimension 2 (subjects)
        %get the average for each iteration & subject accross all the
        %trials of each condition
        mean_within_subject_PH(aa,ab,:) = ...
            mean(random_LPP_PH_group1(aa,ab,:));
        mean_within_subject_PO(aa,ab,:) = ...
            mean(random_LPP_PO_group1(aa,ab,:));
        mean_within_subject_PL(aa,ab,:) = ...
            mean(random_LPP_PL_group1(aa,ab,:));
        mean_within_subject_NEU1(aa,ab,:) = ...
            mean(random_LPP_NEU_group1(aa,ab,:));
        mean_within_subject_NEU2a(aa,ab,:) = ...
            mean(random_LPP_NEU_group2a(aa,ab,:));
        mean_within_subject_NEU2b(aa,ab,:) = ...
            mean(random_LPP_NEU_group2b(aa,ab,:));
    end
end

%contrast the conditions
PH_minus_NEU = mean_within_subject_PH - mean_within_subject_NEU1;
PO_minus_NEU = mean_within_subject_PO - mean_within_subject_NEU1;
PL_minus_NEU = mean_within_subject_PL - mean_within_subject_NEU1;
NEU_minus_NEU = mean_within_subject_NEU2a - mean_within_subject_NEU2b;

%reshape the data into a format that is easier to read 
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects*2,1]);
reshape_PH_minus_NEU = reshape(PH_minus_NEU', ...
    [iterations*subjects, 1]);
reshape_PO_minus_NEU = reshape(PO_minus_NEU', ...
    [iterations*subjects, 1]);
reshape_PL_minus_NEU = reshape(PL_minus_NEU', ...
    [iterations*subjects, 1]);
reshape_NEU_minus_NEU = reshape(NEU_minus_NEU', ...
    [iterations*subjects, 1]);

%create an array for group assignment
grouponearray = ones([iterations*subjects, 1]);
grouptwoarray = repmat(2, [iterations*subjects, 1]);
group_array = [grouponearray; grouptwoarray];

%stack the groups on top of eachother for each contrast
PH_vs_NEU_array = [reshape_PH_minus_NEU; reshape_NEU_minus_NEU];
PO_vs_NEU_array = [reshape_PO_minus_NEU; reshape_NEU_minus_NEU];
PL_vs_NEU_array = [reshape_PL_minus_NEU; reshape_NEU_minus_NEU];

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

%write the results to table format
PH_vs_NEU_between_table = table(reshape_random_subjects, PH_vs_NEU_array, ...
    group_array,'VariableNames', {'subjectID','LPP','group'});

%assign a filename
PH_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'PHvsNEU_between_', ...
    iterations_string, '_iterations.csv'];

%write results to table
writetable(PH_vs_NEU_between_table, PH_vs_NEU_between_filename)

%same thing for PO
PO_vs_NEU_between_table = table(reshape_random_subjects, PO_vs_NEU_array, ...
    group_array,'VariableNames', {'subjectID','LPP','group'});
PO_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'POvsNEU_between_', ...
    iterations_string, '_iterations.csv'];
writetable(PO_vs_NEU_between_table, PO_vs_NEU_between_filename)

%sample thing for PL
PL_vs_NEU_between_table = table(reshape_random_subjects, PL_vs_NEU_array, ...
    group_array,'VariableNames', {'subjectID','LPP','group'});
PL_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'PLvsNEU_between_', ...
    iterations_string, '_iterations.csv'];
writetable(PL_vs_NEU_between_table, PL_vs_NEU_between_filename)
end