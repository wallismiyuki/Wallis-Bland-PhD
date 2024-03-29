#######################################################
# This is an introductory example script for the 
# R package "BioGeoBEARS" by Nick Matzke
# 
# All scripts are copyright Nicholas J. Matzke, 
# please cite if you use. License: GPL-3
# http://cran.r-project.org/web/licenses/GPL-3
# 
# I am happy to answer questions at matzke@nimbios.org, but
# I am more happy to answer questions on the 
# BioGeoBEARS google group
#
# The package is designed for ML and Bayesian inference
# of 
# 
# (a) ancestral geographic ranges, and 
# 
# (b) perhaps more importantly, models for the 
#     evolution of geographic range across a phylogeny.
#
# The example below implements and compares:
# 
# (1) The standard 2-parameter DEC model implemented in 
#     the program LAGRANGE (Ree & Smith 2008); users will
#     notice that the ML parameter inference and log-
#     likelihoods are identical
#
# (2) A DEC+J model implemented in BioGeoBEARS, wherein
#     a third parameter, j, is added, representing the 
#     relative per-event weight of founder-event / jump
#     speciation events at cladogenesis events.  The 
#     higher j is, the more probability these events have,
#     and the less probability the standard LAGRANGE
#     cladogenesis events have.
#
# (3) Some standard model-testing (LRT and AIC) is 
#     implemented at the end so that users may compare models
#
# (4) The script does similar tests of a DIVA-like model (Ronquist 1997)
#     and a BAYAREA-like model (Landis, Matzke, Moore, & Huelsenbeck, 2013)
# 
#######################################################

#######################################################
# Installing BioGeoBEARS from GitHub latest version
#######################################################
# CUT 2018: The old instructions to source() online upgrade .R files have been deleted,
#         all updates are now on the GitHub version of the package, version 1.1+
#######################################################
# Paste the stuff below, INSIDE the single-quote (') marks
# but NOT the single-quote marks themselves.
#install_cmds_that_work_as_of_2023='

# Installation-from-scratch commands
#install.packages("devtools")
#install.packages("ape")
#install.packages("optimx")
#install.packages("GenSA")
#install.packages("rexpokit")   
#install.packages("cladoRcpp")
#install.packages("snow")
#install.packages("MultinomialCI")

#library(devtools)
#devtools::install_github(repo="nmatzke/BioGeoBEARS")

# NOTE: If you get a message like this
# * select "2. CRAN packages only" on "3. None"
# * If you get asked about "binaries" vs. "source", choose "binaries" 
#   (binaries are precompiled and easy to install; installing from source
#    requires that your computer have the correct compilers, which can be
#    challenging if you are not fairly expert)
#' # END install_cmds_that_work_as_of_2023

#######################################################
#######################################################

#######################################################
# SETUP -- libraries/BioGeoBEARS updates
#######################################################

# Load the package (after installation, see above).
library(ape)
library(optimx)   # optimx seems better than R's default optim()
library(GenSA)    # GenSA seems better than optimx (but slower) on 5+ parameters, 
                  # seems to sometimes fail on simple problems (2-3 parameters)
library(rexpokit)
library(cladoRcpp)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)
library(BioGeoBEARS)

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

# Set your working directory for output files
# default here is your home directory ("~")
# Change this as you like
wd = np("/Users/wallis/Documents/GitHub/Files/Tester/M0_tests")
setwd(wd)
file_list = list.files(wd)
file_list
numbers = c(2,3,17)

# Double-check your working directory with getwd()
getwd()
Example = NULL
N_ranges = NULL
Max_Range_Size = NULL
Tip_Nums = NULL
Model = NULL
Start_times = NULL
End_times = NULL
Time_taken = NULL


#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory 
# inside the installed package.
#
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# 
# Each computer operating system might install BioGeoBEARS in a different place, 
# depending on your OS and settings. 
# 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the 
# path to the format appropriate for your system (e.g., Mac/Linux use "/", but 
# Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your 
# script, since the plot_BioGeoBEARS_results function needs a script from the 
# extdata directory to calculate the positions of "corners" on the plot. This cannot
# be made into a straight up BioGeoBEARS function because it uses C routines 
# from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Example files are given below. To run your own data,
# make the below lines point to your own files, e.g.
# trfn = "/mydata/frogs/frogBGB/tree.newick"
# geogfn = "/mydata/frogs/frogBGB/geog.data"

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################
# This is the example Newick file for Hawaiian Psychotria
# (from Ree & Smith 2008)
# "trfn" = "tree file name"
for(i in numbers){

    setwd(wd)
    setwd(file_list[i])

    #trfn = np(paste(addslash(file_list[i]), "tree.newick", sep=""))
    trfn = "tree.newick"
    trfn
    # Look at the raw Newick file:

    moref(trfn)

    # Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
    pdffn = "tree.pdf"
    pdf(file=pdffn, width=9, height=12)

    tr = read.tree(trfn)
    tr
    plot(tr)
    title(paste0("Tree from folder ", file_list[i]))
    axisPhylo() # plots timescale

    dev.off()
    ##cmdstr = paste0("open ", pdffn)
    ##system(cmdstr)

    #######################################################
    # Geography file
    # Notes:
    # 1. This is a PHYLIP-formatted file. This means that in the 
    #    first line, 
    #    - the 1st number equals the number of rows (species)
    #    - the 2nd number equals the number of columns (number of areas)
    #    - after a tab, put the areas in parentheses, with spaces: (A B C D)
    #
    # 1.5. Example first line:
    #    10    4    (A B C D)
    # 
    # 2. The second line, and subsequent lines:
    #    speciesA    0110
    #    speciesB    0111
    #    speciesC    0001
    #         ...
    # 
    # 2.5a. This means a TAB between the species name and the area 0/1s
    # 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
    # 
    # 3. See example files at:
    #    http://phylo.wikidot.com/biogeobears#files
    # 
    # 4. Make you understand what a PLAIN-TEXT EDITOR is:
    #    http://phylo.wikidot.com/biogeobears#texteditors
    #
    # 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
    #
    # 4. All names in the geography file must match names in the phylogeny file.
    #
    # 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
    #
    # 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
    #    i.e. genetically isolated populations.  These may or may not be identical 
    #    with species.  You would NOT want to just use specimens, as each specimen 
    #    automatically can only live in 1 area, which will typically favor DEC+J 
    #    models.  This is fine if the species/lineages really do live in single areas,
    #    but you wouldn't want to assume this without thinking about it at least. 
    #    In summary, you should collapse multiple specimens into species/lineages if 
    #    data indicates they are the same genetic population.
    ######################################################

    # This is the example geography file for Hawaiian Psychotria
    # (from Ree & Smith 2008)
    #geogfn = np(paste(addslash(extdata_dir), "Psychotria_geog.data", sep=""))
    geogfn = "geog_data.txt"

    # Look at the raw geography text file:
    moref(geogfn)

    # Look at your geographic range data:
    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
    tipranges
    
    # Maximum range size observed:
    max(rowSums(dfnums_to_numeric(tipranges@df)))
    n_ranges = ncol(tipranges@df)

    # Set the maximum number of areas any species may occupy; this cannot be larger 
    # than the number of areas you set up, but it can be smaller.
     
    max_range_size = n_ranges
    max_range_size    

    ####################################################
    ####################################################
    # KEY HINT: The number of states (= number of different possible geographic ranges)
    # depends on (a) the number of areas and (b) max_range_size.
    # If you have more than about 500-600 states, the calculations will get REALLY slow,
    # since the program has to exponentiate a matrix of e.g. 600x600.  Often the computer
    # will just sit there and crunch, and never get through the calculation of the first
    # likelihood.
    # 
    # (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
    #
    # To check the number of states for a given number of ranges, try:
    numstates_from_numareas(numareas=11, maxareas=11, include_null_range=TRUE)
    numstates_from_numareas(numareas=11, maxareas=8, include_null_range=TRUE)
    numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
    numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)

    # Large numbers of areas have problems:
    numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)

    # ...unless you limit the max_range_size:
    numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
    ####################################################
    ####################################################

    #######################################################
    #######################################################
    # DEC AND DEC+J ANALYSIS
    #######################################################
    #######################################################
    # NOTE: The BioGeoBEARS "DEC" model is identical with 
    # the Lagrange DEC model, and should return identical
    # ML estimates of parameters, and the same 
    # log-likelihoods, for the same datasets.
    #
    # Ancestral state probabilities at nodes will be slightly 
    # different, since BioGeoBEARS is reporting the 
    # ancestral state probabilities under the global ML
    # model, and Lagrange is reporting ancestral state
    # probabilities after re-optimizing the likelihood
    # after fixing the state at each node. These will 
    # be similar, but not identical. See Matzke (2014),
    # Systematic Biology, for discussion.
    #
    # Also see Matzke (2014) for presentation of the 
    # DEC+J model.
    #######################################################
    #######################################################

    #######################################################
    #######################################################

    #######################################################
    # Run DEC
    #######################################################

    # Intitialize a default model (DEC model)
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()

    # Give BioGeoBEARS the location of the phylogeny Newick file
    BioGeoBEARS_run_object$trfn = trfn

    # Give BioGeoBEARS the location of the geography text file
    BioGeoBEARS_run_object$geogfn = geogfn

    # Input the maximum range size
    BioGeoBEARS_run_object$max_range_size = max_range_size

    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change

    # Set up a time-stratified analysis:
    # 1. Here, un-comment ONLY the files you want to use.
    # 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
    # 3. For example files see (a) extdata_dir, 
    #  or (b) http://phylo.wikidot.com/biogeobears#files
    #  and BioGeoBEARS Google Group posts for further hints)
    #
    # Uncomment files you wish to use in time-stratified analyses:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
    # if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
    # problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
    BioGeoBEARS_run_object$num_cores_to_use = 1
    # (use more cores to speed it up; this requires
    # library(parallel) and/or library(snow). The package "parallel" 
    # is now default on Macs in R 3.0+, but apparently still 
    # has to be typed on some Windows machines. Note: apparently 
    # parallel works on Mac command-line R, but not R.app.
    # BioGeoBEARS checks for this and resets to 1
    # core with R.app)

    # Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
    # I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
    # but the results are imprecise and so I haven't explored it further.
    # In a Bayesian analysis, it might work OK, but the ML point estimates are
    # not identical.
    # Also, I have not implemented all functions to work with force_sparse=TRUE.
    # Volunteers are welcome to work on it!!
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table

    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

    # Set up DEC model
    # (nothing to do; defaults)

    # Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
    BioGeoBEARS_run_object

    # This contains the model object
    BioGeoBEARS_run_object$BioGeoBEARS_model_object

    # This table contains the parameters of the model 
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

    # Run this to check inputs. Read the error messages if you get them!
    BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)

    # For a slow analysis, run once, then set runslow=FALSE to just 
    # load the saved result.
    runslow = TRUE
    resfn = "DEC_Output.Rdata"

    Model <- append(Model, "DEC")
    Example  <- append(Example, file_list[i])
    N_ranges <- append(N_ranges, n_ranges)
    Max_Range_Size = append(Max_Range_Size, max_range_size)
    Tip_Nums = append(Tip_Nums, dim(tipranges@df)[1])

    if (runslow)
        {
        start.time <- Sys.time()
        res = bears_optim_run(BioGeoBEARS_run_object)

        res    
        
        end.time <- Sys.time()
        
        save(res, file=resfn)
        
        
        time.taken <- end.time - start.time
        time.taken
        Start_times <- append(Start_times, start.time)
        End_times <- append(End_times, end.time)
        Time_taken <- append(Time_taken, time.taken)
        write.table(c(file_list[i],"DEC", n_ranges, max_range_size, dim(tipranges@df)[1], time.taken, start.time, end.time), "DEC_values.txt")
        
        speedtests = data.frame(Example, Model, N_ranges, Max_Range_Size, Tip_Nums, Time_taken, Start_times, End_times)
        write.csv(speedtests, "/Users/wallis/Documents/GitHub/Files/Tester/R_Speedtestpt3.csv", row.names = FALSE)

        resDEC = res
        } else {
        # Loads to "res"
        load(resfn)
        resDEC = res
        }




    #######################################################
    # Run DEC+J
    #######################################################
    BioGeoBEARS_run_object = define_BioGeoBEARS_run()
    BioGeoBEARS_run_object$trfn = trfn
    BioGeoBEARS_run_object$geogfn = geogfn
    BioGeoBEARS_run_object$max_range_size = max_range_size
    BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
    BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
    # (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
    #  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
    #  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
    #  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
    # Also: search script on "include_null_range" for other places to change

    # Set up a time-stratified analysis:
    #BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
    #BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
    #BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
    #BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
    #BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
    # See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

    # Speed options and multicore processing if desired
    BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
    BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
    BioGeoBEARS_run_object$use_optimx = TRUE    # if FALSE, use optim() instead of optimx();
    # if "GenSA", use Generalized Simulated Annealing, which seems better on high-dimensional
    # problems (5+ parameters), but seems to sometimes fail to optimize on simple problems
    BioGeoBEARS_run_object$num_cores_to_use = 1
    BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

    # This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
    # (It also runs some checks on these inputs for certain errors.)
    BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

    # Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
    #BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE, fossils_older_than=0.001, cut_fossils=FALSE)
    # The stratified tree is described in this table:
    #BioGeoBEARS_run_object$master_table

    # Good default settings to get ancestral states
    BioGeoBEARS_run_object$return_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
    BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

    # Set up DEC+J model
    # Get the ML parameter values from the 2-parameter nested model
    # (this will ensure that the 3-parameter model always does at least as good)
    dstart = resDEC$outputs@params_table["d","est"]
    estart = resDEC$outputs@params_table["e","est"]
    jstart = 0.0001

    # Input starting values for d, e
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

    # Add j as a free parameter
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
    BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

    BioGeoBEARS_run_object = fix_BioGeoBEARS_params_minmax(BioGeoBEARS_run_object=BioGeoBEARS_run_object)
    check_BioGeoBEARS_run(BioGeoBEARS_run_object)

    resfn = "DEC_J_Output.Rdata"
    runslow = TRUE

    Model <- append(Model, "DECJ")
    Example  <- append(Example, file_list[i])
    N_ranges <- append(N_ranges, n_ranges)
    Max_Range_Size = append(Max_Range_Size, max_range_size)
    Tip_Nums = append(Tip_Nums, dim(tipranges@df)[1])

    if (runslow)
        {
        #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
        start.time <- Sys.time()
        res = bears_optim_run(BioGeoBEARS_run_object)
        

        res    
        
        end.time <- Sys.time()

        save(res, file=resfn)
        time.taken <- end.time - start.time
        time.taken
        Start_times <- append(Start_times, start.time)
        End_times <- append(End_times, end.time)
        Time_taken <- append(Time_taken, time.taken)
        write.table(c(file_list[i],"DECJ", n_ranges, max_range_size, dim(tipranges@df)[1], time.taken, start.time, end.time), "DECJ_values.txt")
        
        speedtests = data.frame(Example, Model, N_ranges, Max_Range_Size, Tip_Nums, Time_taken, Start_times, End_times)
        write.csv(speedtests, "/Users/wallis/Documents/GitHub/Files/Tester/R_speedtestpt3.csv", row.names = FALSE)

        

        resDECj = res
        } else {
        # Loads to "res"
        load(resfn)
        resDECj = res
        }

}


Example
Model
Start_times
End_times
Time_taken



