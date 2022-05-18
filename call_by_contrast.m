function call_by_contrast(test_mode, optstable_filename, opts_table,  ...
    datadir, outputdir)

%Matlab caller function to wrap around loop_LPP_by_contrast within and 
%between authored by Kyla Gibney: kdgibney.work@gmail.com

%%
%INS

%test_mode: option to enable test mode. 1 = test mode, 0 = not.

%optstable_filename: the filename of the .csv file containting the 
%opts table. Use this if you do not already have the opts table already in
%the workspace

%opts_table: the opts table itself

%datadir: the directory containing the raw data

%outputdir: the directory to write results to


%%
%OUTS:

%there are no outs for this caller function. loop_LPP_by_contrast
%within and between for each contrast return their own outs

%%
%DEFAULTS

%if you do not pass in any input arguments, the function will default to
%these

if test_mode == 1 %set opts for running in test mode
    subjects = 10; %10 subjects (between: 10 subjects per group)
    trials = 5; %5 trial per condition
    effectsize = 0; %null effect size (zero microvolts)
    iterations = 100; %100 iterations
    opts_table = table(subjects, trials, effectsize, iterations, ...
        'VariableNames', {'subjects', 'trials', ...
        'effectsize','iterations'}); %create opts table for 
    %running in test mode
end

if ~exist('test_mode','var') %if you don't pass in a test_mode variable,
    %it defaults to zero
    test_mode = 0;
end

if ~exist('optstable_filename','var')
    optstable_filename = ['C:\path_to_opts_table']; %if you don't pass in an
    %opts file, use the path to the default opts file here
end

if ~exist('opts_table','var')
    opts_table = readtable(optstable_filename); %if you didn't pass in an 
    %opts table, read table from default file above
end
 
if ~exist('datadir','var')
    datadir = 'C:\path_to_data_dir\'; %if you didn't specify a datadir,
    %use the path to the default data dir here
end

if ~exist('outputdir','var')
    outputdir = 'C:\path_to_output_dir\'; %if you didn't specify a datadir,
    %use the path to the default output dir here
end
 

%begin to loop through simulations
for i = 1:size(opts_table,1)
    
    %parse experiment info for within-subjects simulations   
    disp_string_within{i} = ['within-subjects experiment # ', num2str(i), ...
        ' out of ', num2str(height(opts_table))];
    
    %display experiment # out of total # within-subjects iterations
    disp(disp_string_within{i});
    
    %parse # subjects, trials, & iterations
    subjects(i) = opts_table.subjects(i);
    trials(i) = opts_table.trials(i);
    iterations(i) = opts_table.iterations(i);
    
    %call the within subjects functions for each contrast
    PLE_vs_NEU_within(test_mode, subjects(i),trials(i), ...
        iterations(i), datadir, outputdir) %within-subjects simulation 
    %for pleasant vs neutral images
    
    UNP_vs_NEU_within(test_mode, subjects(i),trials(i), ...
        iterations(i), datadir, outputdir) %within-subjects simulation 
    %for unpleasant vs neutral iamges
    
    CIG_vs_NEU_within(test_mode, subjects(i),trials(i), ...
        iterations(i), datadir, outputdir) %within-subjects simulation 
    %for cigarette images vs neutral images
    
    NEU_vs_NEU_within(test_mode, subjects(i),trials(i), ...
        iterations(i), datadir, outputdir) %within-subjects simulation 
    %for neutral vs neutral images (control/null contrast)

    %display experiment # out of total # between-subjects iterations
    disp_string_between{i} = ['between-subjects experiment # ', ...
    num2str(i), ' out of ', num2str(height(opts_table))];
    disp(disp_string_between{i})

    %call the between-subjects functions for each contrast
    PLE_vs_NEU_between(test_mode, subjects(i),trials(i), ...
        iterations(i), datadir, outputdir) %between-subjects simulation 
    %for pleasant vs neutral images
    
    UNP_vs_NEU_between(test_mode, subjects(i), trials(i), ...
        iterations(i), datadir, outputdir)%between-subjects simulation 
    %for unpleasant vs neutral images
    
    CIG_vs_NEU_between(test_mode, subjects(i), trials(i), ...
        iterations(i), datadir, outputdir)%between-subjects simulation 
    %for cigarette vs neutral images
    
    NEU_vs_NEU_between(test_mode, subjects(i), trials(i), ...
        iterations(i), datadir, outputdir)%between-subjects simulation 
    %for neutral vs neutral images
    
end

end