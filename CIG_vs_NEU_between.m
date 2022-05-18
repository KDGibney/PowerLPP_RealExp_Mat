function CIG_vs_NEU_between(test_mode, subjects, trials, ...
    iterations, datadir, outputdir)

%Matlab function for between-subjects power calculations of the LPP  
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


if ~exist('subjects','var') %if you don't pass in a subjects variable,
    %it defaults to ten
    subjects = 10;
end

if ~exist('trials','var') %if you don't pass in a trials variable,
    %it defaults to ten
    trials = 10;
end

if ~exist('iterations','var') %if you don't pass in an iterations variable,
    %it defaults to one hundred
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


%get CIG (cigarette images) data
CIG_file = [datadir, 'CIG_table.csv'];
CIG_table = readtable(CIG_file);
CIG_table = sortrows(CIG_table,1); %sort rows by subject ID


%get NEU (neutral images) data
NEU_file = [datadir, 'NEU_table.csv'];
NEU_table = readtable(NEU_file);
NEU_table = sortrows(NEU_table, 1); %sort rows by subject ID

unique_subjid_CIG = unique(CIG_table.subID); %get unique subjects
unique_subjid_NEU = unique(NEU_table.subID);
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%get number of trials per category
num_CIG_trials = length(find(unique_subjid_CIG(1) == CIG_table.subID));

total_subjects = length(unique_subjid_NEU); %find total # subjects

%make sure we have the same # subjects in each category
if ~isequal(unique_subjid_CIG,unique_subjid_NEU)
    disp('different sets of subjects!')
end

%make an index of trials for each category
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]); 
trial_index_CIG = repmat(transpose(1:num_CIG_trials),[total_subjects 1]); 

for i = 1:iterations %tick through all iterations
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects*2 %iterate through all subjects
        %we're sampling 2x subjects to simulate 2 experimental groups
        
        random_subjects(i,:) = datasample(unique_subjid_CIG,subjects*2,...
            'REPLACE',false);
        %randomly sample 2x subjects without replacement
        
        for k = 1:trials*2 %tick through all trials
            %we sample 2x NEU trials to model a difference between 
            %conditions between groups (CIG vs NEU vs NEU vs NEU)
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials,...
                trials*2, 'REPLACE',true);
            %randomly sample neutral trials with replacement
            %we sample WITH replacement for the trials in the real
            %experiments, because in subcategories (pleasant high, etc) we
            %only have 16 trials total but we sample up to 15, so sampling
            %with replacement is necessary
            
            find_trial_index_NEU(i,j,k) = find(random_trials_NEU(i,j,k) ...
                == trial_index_NEU & random_subjects(i,j) ...
                == NEU_table.subID);
            %find the indices of all the trials we sampled in NEU_table
            
            random_LPP_NEU(i,j,k) = NEU_table.LPP(find_trial_index_NEU(i,j,k)); 
            %find the NEU LPPs that match that index
            
            for l = 1:trials %tick through trials
                %we sample 1x trials for the CIG category
                
                random_trials_CIG(i,j,:) = datasample(1:num_CIG_trials, ...
                    trials,'REPLACE',true); %randomly sample CIG trials
                %with replacement
                
                find_trial_index_CIG(i,j,l) = ...
                    find(random_trials_CIG(i,j,l) == ...
                    trial_index_CIG & random_subjects(i,j) == ...
                    CIG_table.subID); %find indices of all the CIG trials
                    %we sampled in CIG_table
                
                random_LPP_CIG(i,j,l) = ...
                    CIG_table.LPP(find_trial_index_CIG(i,j,l));
                    %find all CIG LPPs          
            end 
        end
    end
    
    %assign LPPs to groups and conditions
    
    %one group for cigarette images vs neutral images
    random_LPP_CIG_group1(i,:,:) = random_LPP_CIG(i,1:subjects,:); 
    random_LPP_NEU_group1(i,:,:) = random_LPP_NEU(i,1:subjects,...
        1:trials);
    
    %one group for neutral vs neutral images
    random_LPP_NEU_group2a(i,:,:) = random_LPP_NEU(i, ...
        (subjects+1):(subjects*2),1:trials);
    random_LPP_NEU_group2b(i,:,:) = random_LPP_NEU(i, ...
        (subjects+1):(subjects*2), ...
        (trials+1):(trials*2));
end

% check to make sure sizes of groups are consistent
if size(random_LPP_CIG_group1) ~= size(random_LPP_NEU_group2b)
    disp('groups are different sizes!')
end

%get size of a group by condition
[a, b, ~] = size(random_LPP_CIG_group1);


for aa = 1:a %tick through dimension 1
    for ab = 1:b %tick through dimension 2
        %get the means accross all the trials for conditions 1 and 2 of
        %each group
        mean_within_subject_CIG(aa,ab,:) = ...
            mean(random_LPP_CIG_group1(aa,ab,:)); 
        mean_within_subject_NEU1(aa,ab,:) = ...
            mean(random_LPP_NEU_group1(aa,ab,:));
        mean_within_subject_NEU2a(aa,ab,:) = ...
            mean(random_LPP_NEU_group2a(aa,ab,:));
        mean_within_subject_NEU2b(aa,ab,:) = ...
            mean(random_LPP_NEU_group2b(aa,ab,:));
    end
end

%contrast cigarette vs neutral
CIG_minus_NEU = mean_within_subject_CIG - mean_within_subject_NEU1;

%contrast neu vs neu
NEU_minus_NEU = mean_within_subject_NEU2a - mean_within_subject_NEU2b;

%reshape the results into an easier to read format 
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects*2,1]);
reshape_CIG_minus_NEU = reshape(CIG_minus_NEU', ...
    [iterations*subjects, 1]);
reshape_NEU_minus_NEU = reshape(NEU_minus_NEU', ...
    [iterations*subjects, 1]);

%make arrays for group designation
grouponearray = ones([iterations*subjects, 1]);
grouptwoarray = repmat(2, [iterations*subjects, 1]);
group_array = [grouponearray; grouptwoarray];

%stack the contrasts on top of each other
CIG_vs_NEU_array = [reshape_CIG_minus_NEU; reshape_NEU_minus_NEU];

%create a test directory if necessary
if test_mode == 1
    outputdir = ['C:\path_to_output_dir\'];
    mkdir(outputdir)
else
    mkdir(outputdir)
end

%make an output directory
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

%write the results to a table
CIG_vs_NEU_between_table = table(reshape_random_subjects, CIG_vs_NEU_array, ...
    group_array,'VariableNames', {'subjectID','LPP','group'});

%assign a filename
CIG_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'CIGvsNEU_between_', ...
    iterations_string, '_iterations.csv'];

%write the table to csv
writetable(CIG_vs_NEU_between_table, CIG_vs_NEU_between_filename)

end

