loop_anova_withinsubjects_realexp <- function(datadir = 'C:/path_to_data_dir',
                                               outputdir = 'C:/path_to_output_dir',
                                               testdatadir = 'C:/path_to_test_data_dir',
                                               testoutputdir = 'C:/path_to_test_output_dir',
                                               test_mode = 0) {
  
  #R function for repeated-measures ANOVAS for monte carlo simulations 
  #for powerlpp real experiments
  #Authored by Kyla Gibney and George Kypriotakis
  #Direct correspondence to: kdgibney.work@gmail.com
  
  
  ###############################################################################################################
  #INS
  
  #datadir: directory that the data to be analyzed is stored in
  
  #outpurdir: directory that the results of the ANOVAs will be stored in
  
  #testdatadir: if test_mode = 1, function will look here for data to be processed
  
  #testoutputdir: if test_mode = 1, function will write results of ANOVAS here
  
  #test_mode: dummy variable to trigger test mode.
  
  ##############################################################################################################
  #OUTS
  
  #this function does not return anything to the workspace, but will write a results file
  #(results_r_betweensubjects.csv') to datadir
  
  #############################################################################################################
  #TEST MODE
  
  #IF YOU ENABLE TEST MODE, RESTULS WILL WRITE TO A TEST DIRECTORY
  
  if(test_mode == 1) {
    datadir = testdatadir
    outputdir = testoutputdir
  }
  
  ###########################################################################################################
  #BEGIN FUNCTION
  
  results_file <- paste(outputdir,'/results_r_withinsubjects_realexp.csv',sep='') #create filename
  
  sink(results_file) # this diverts R output to the connection provided by the file to write output
  
  files <- list.files(datadir, pattern="*.csv", full.names = T) #list all files in datadir
  
  num_files = seq_along(files) #total number of files in datadir
  
  
  for (i in num_files){ #tick through each file in datadir
  
  	data.i <- files[i] #current file being processed
  
  	data.current <- read.csv(data.i)  #read the current data file
    
  	parsefilename = strsplit(data.i, "_subjects")
  	parsesubjects = strsplit(parsefilename[[1]][1],'/real_results_within/')
  	subjects = as.numeric(parsesubjects[[1]][2]) #get number of subjects from file name
  
  	parsetrialnumber = strsplit(data.i, "_trials")	
  	parsetrials = strsplit(parsetrialnumber[[1]][1],"_subjects_")
  	trials = as.numeric(parsetrials[[1]][2]) #get number of trials from file name
  
  	parsecontrast = strsplit(data.i,"_iterations")
  	parseiterations = strsplit(parsecontrast[[1]][1],'_trials_')
  	parseiterations2 = strsplit(parseiterations[[1]][2],"_")
  	simulation_contrast = parseiterations2[[1]][1] #get contrast from file name
  	simulation_iterations = as.numeric(parseiterations2[[1]][2]) #get iterations from file name
  
  
  	num_steps = subjects*2#number of subjects to be taken at a time
  	numrows = nrow(data.current) #total number of rows in the current dataset
  	index1 = seq(1,numrows,num_steps) #create an index starting at one in steps of 10 until 
  	#the end of the dataset
  	index2 = seq(num_steps,numrows,num_steps)#create an index starting at 10 in steps of 10 
  	#until the end of the dataset
  	iterations = numrows/num_steps #calculate the number of iterations based on the file size
  	#this should match the iterations we got from the file name
    
  	for (j in 1:iterations){ #tick through each iterations of the current file
    
  		Rep_M_Anova <- aov(LPP[index1[j]:index2[j]] ~ factor(condition[index1[j]:index2[j]]) + Error(factor(subjectID[index1[j]:index2[j]])), data=data.current) # here you can run your model,
  		# y (LPP) predicted by x (condition) plus error (subjectID)
  		# Error= experimental design - here repeated measures).  We can change that to any model.
  		sum_test = unlist(summary(Rep_M_Anova)) #list the results of the ANOVA
  
      #these are all the outputs we want printed in our results file  
  		print(paste(subjects,"subjects"))
  		print(paste(trials,"trials"))
  		print(paste(simulation_contrast))
  		print(paste(simulation_iterations,"iterations"))
  		print(paste("F val"))
  		print(paste(sum_test[12]))
  		print(paste("P val"))
  		print(paste(sum_test[14]))
  
  }}
  
  sink() #ending the diversion
}