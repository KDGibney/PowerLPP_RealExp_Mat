function NEU_vs_NEU_between(test_mode, subjects, trials, iterations, ...
    datadir, outputdir)

%Matlab function for between-subjects power calculations of the LPP  
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

unique_subjid_NEU = unique(NEU_table.subID); %get all the subIDs
num_NEU_trials = length(find(unique_subjid_NEU(1) == NEU_table.subID));
%find the number of trials in the NEU condition
total_subjects = length(unique_subjid_NEU); %get total # subjects in the 
%NEU condition

%create an index of trials
trial_index_NEU = repmat(transpose(1:num_NEU_trials),[total_subjects 1]); 

for i = 1:iterations %tick through all the iterations
    
    print_string = ['outer loop number ', num2str(i)];
    disp(print_string) %print which loop we're on
    
    for j = 1:subjects*2 %tick through all subjects
        %we're sampling 2x subjects to simulate 2 experimental groups
        
        random_subjects(i,:) = datasample(unique_subjid_NEU,subjects*2, ...
            'REPLACE',false);%randomly sample 2x subjects without 
            %replacement
            
        for k = 1:trials*2 %tick through all trials
            %we're sampling 2x trials to simulate to conditions
            
            random_trials_NEU(i,j,:) = datasample(1:num_NEU_trials, ...
                trials*2, 'REPLACE',true); 
            %randomly sample neutral trials with replacement
            %we sample WITH replacement for the trials in the real
            %experiments, because in subcategories (pleasant high, etc) we
            %only have 16 trials total but we sample up to 15, so sampling
            %with replacement is necessary
            
            find_trial_index_NEU(i,j,k) = ...
                find(random_trials_NEU(i,j,k) == trial_index_NEU ...
                & random_subjects(i,j) == NEU_table.subID);
                %find the index of all the neutral trials we sampled
            
            random_LPP_NEU(i,j,k) = ...
                NEU_table.LPP(find_trial_index_NEU(i,j,k));
                %get the LPPs that match those indices
        end
    end
    
    %assign groups and conditions
    random_LPP_NEU_group1C1(i,:,:) = random_LPP_NEU(i,1:subjects,1:trials);
    random_LPP_NEU_group1C2(i,:,:) = random_LPP_NEU(i,1:subjects,(trials+1):(trials*2));
    random_LPP_NEU_group2C1(i,:,:) = random_LPP_NEU(i,(subjects+1):(subjects*2),1:trials);
    random_LPP_NEU_group2C2(i,:,:) = random_LPP_NEU(i,(subjects+1):(subjects*2),(trials+1):(trials*2));
end

% check to make sure the sizes of the groups/conditions are
%consistent
if size(random_LPP_NEU_group1C1) ~= size(random_LPP_NEU_group2C2)
    disp('groups are different sizes!')
end

%get the dimensions of the results
[dim1, dim2, ~] = size(random_LPP_NEU_group1C1);

for aa = 1:dim1 %tick through first dimension (iterations)
    
    for bb = 1:dim2 %tick through second dimension (subjects)
        
        %get the average of each subject, each iteration accross all trials
        %of each condition
        avg_LPP_NEU_group1C1(aa,bb) = mean(random_LPP_NEU_group1C1(aa,bb,:));
        avg_LPP_NEU_group1C2(aa,bb) = mean(random_LPP_NEU_group1C2(aa,bb,:));
        avg_LPP_NEU_group2C1(aa,bb) = mean(random_LPP_NEU_group2C1(aa,bb,:));
        avg_LPP_NEU_group2C2(aa,bb) = mean(random_LPP_NEU_group2C2(aa,bb,:));
    end
end

%get the contrasts for each group
NEU_minus_NEU_1 = avg_LPP_NEU_group1C1 - avg_LPP_NEU_group1C2;
NEU_minus_NEU_2 = avg_LPP_NEU_group2C1 - avg_LPP_NEU_group2C2;

%reshape the data so it's easier to read 
reshape_random_subjects = reshape(random_subjects', ...
    [iterations*subjects*2,1]);
reshape_NEU_minus_NEU1 = reshape(NEU_minus_NEU_1', ...
    [iterations*subjects, 1]);
reshape_NEU_minus_NEU2 = reshape(NEU_minus_NEU_2', ...
    [iterations*subjects, 1]);

%greate an array for group assignment
grouponearray = ones([iterations*subjects, 1]);
grouptwoarray = repmat(2, [iterations*subjects, 1]);
group_array = [grouponearray; grouptwoarray];

%stack the two groups on top of eachother
NEU_vs_NEU_array = [reshape_NEU_minus_NEU1; reshape_NEU_minus_NEU2];

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

%write the results to a table
NEU_vs_NEU_between_table = table(reshape_random_subjects, ...
    NEU_vs_NEU_array, group_array,'VariableNames', ...
    {'subjectID','LPP','group'});

%assign a filename
NEU_vs_NEU_between_filename = [between_results_dir, subjects_string, ...
    '_subjects_', trials_string, '_trials_', 'NEUvsNEU_between_', ...
    iterations_string, '_iterations.csv'];

%write table to csv
writetable(NEU_vs_NEU_between_table, NEU_vs_NEU_between_filename)
end