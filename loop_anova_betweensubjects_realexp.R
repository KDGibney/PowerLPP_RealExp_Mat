loop_anova_betweensubjects_realexp <- function(datadir = 'C:/path_to_data_dir',
                                                   outputdir = 'C:/path_to_output_dir',
                                                   testdatadir = 'C:/path_to_test_data_dir',
                                                   testoutputdir = 'C:/path_to_test_output_dir',
                                                   test_mode = 0) {
  
  #R function for ANOVAs of between-subject monte carlo simulations for powerlpp real experiments
  #Authored by Kyla Gibney and George Kypriotakis
  #Direct correspondence to: kdgibney.work@gmail,com
  
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
  results_file <- paste(outputdir,'/results_r_betweensubjects_realexp.csv',sep='') #create filename
  
  sink(results_file) # this diverts R output to the connection provided by the file to write output
  
  files <- list.files(datadir, pattern="*.csv", full.names = T) #list all files in datadir
  
  num_files = seq_along(files) #total number of files in datadir
  
  
  for (i in num_files){ #loop through all files
  
  	data.i <- files[i] #current file being processed
  
  	data.current <- read.csv(data.i)  #read the current data file
    
  	parsefilename = strsplit(data.i, "_subjects")
  	parsesubjects = strsplit(parsefilename[[1]][1],"real_results_between/")
  	subjects = as.numeric(parsesubjects[[1]][2]) #get number of subjects from file name
  
  	parsetrialnumber = strsplit(data.i, "_trials")	
  	parsetrials = strsplit(parsetrialnumber[[1]][1],"_subjects_")
  	trials = as.numeric(parsetrials[[1]][2]) #get number of trials from file name
  
  	parsecontrast = strsplit(data.i,"_iterations")
  	parseiterations = strsplit(parsecontrast[[1]][1],'_trials_')
  	parseiterations2 = strsplit(parseiterations[[1]][2],"_between_")
  	simulation_contrast = parseiterations2[[1]][1] #get contrast from file name
  	simulation_iterations = as.numeric(parseiterations2[[1]][2]) #get iteratiosn from file name
  	
  
  
  	num_steps = subjects #number of subjects to be taken at a time
  	numrows = nrow(data.current) #total number of rows in our data
  	halfpoint = numrows/2 #the end of group one
  	start = halfpoint+1 #the start of group 2
  	finish = numrows #the end of group 2
  	LPP_group_one <- data.current$LPP[1:halfpoint] #assign LPPs to groups
  	LPP_group_two <- data.current$LPP[start:finish]
  	group_one = data.current$group[1:halfpoint] #organize group assignment labels
  	group_two = data.current$group[start:finish]
  	size1 <- length(LPP_group_one) #get the size of each group
  	size2 <- length(LPP_group_two)
  	index1 = seq(1,size1,num_steps) #create an index starting at one in steps of 10 until 
  	#the end of the dataset
  	index2 = seq(num_steps,size2,num_steps) #create an index starting at 10 in steps of 10 
  	#until the end of the dataset
  	iterations = numrows/(num_steps*2) #determine the iterations mathematically
  	#this should match the iterations we got from the file name
  
    
  	for (j in 1:simulation_iterations){ #tick through every iteration
  		combined_LPPs <- c(LPP_group_one[index1[j]:index2[j]],LPP_group_two[index1[j]:index2[j]])
  		#concatenate the LPPs of group 1 and 2, iteration 1 together and so on
  		combined_groups <- c(group_one[index1[j]:index2[j]],group_two[index1[j]:index2[j]])
  		#concatenate the groups labels of group 1 and 2, iteration 1 and so on
  		
  		lppdata <- data.frame(combined_LPPs,combined_groups) #convert to a data frame
  		anova = aov(combined_LPPs ~ combined_groups, data=lppdata) #y (LPP) predicted by x 
  		#(group assignment)		
  		names = unlist(summary(anova)) #list the results of the anova
  
      #these are all the outputs we want printed in our results file
  		print(paste(subjects,"subjects"))
  		print(paste(trials,"trials"))
  		print(paste(simulation_contrast,"contrast"))
  		print(paste(simulation_iterations,"iterations"))
  		print("F value")
  		print(paste(names[7]))
  		print("P value")
  		print(paste(names[9]))
      }
    }
  
  sink() #end the diversion
  }