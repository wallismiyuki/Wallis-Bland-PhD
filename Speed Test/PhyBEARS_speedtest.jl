
#######################################################
# Example inference on a simulated tree & geography dataset
# (from Wallis's ss8_sim_001)
# 2022-11-09
#
# We will compare inference under 
# * SSE likelihood calculator v7 (constant rates through time, very fast)
# * SSE likelihood calculator v12 (changing rates through time, very fast)
#   - Here, to start, we will only change the distances through time
#######################################################

using Interpolations	# for Linear, Gridded, interpolate
using LinearAlgebra  	# for "I" in: Matrix{Float64}(I, 2, 2)
										 	# https://www.reddit.com/r/Julia/comments/9cfosj/identity_matrix_in_julia_v10/
using Sundials				# for CVODE_BDF
using Test						# for @test, @testset
using PhyloBits
using PhyBEARS
using DataFrames
using CSV
using Printf
using DelimitedFiles
using NLopt
using StatsBase
using Combinatorics
using BenchmarkTools

using DifferentialEquations
using PhyBEARS.TimeDep
using PhyBEARS.Uppass
using PhyBEARS.StateSpace
using PhyBEARS.TreePass
using PhyBEARS.SSEs
using PhyBEARS.ModelLikes
using PhyBEARS.Optimizers
using PhyBEARS.BGExample
using PhyBEARS.Parsers

using PhyloBits.PNtypes						# Types: Tree, Node, Edge etc.
using PhyloBits.PNmanipulateNet
using PhyloBits.PNreadwrite				# Helper functions
using PhyloBits.PNreadwrite 			# Reading and writing trees; readTopology
using PhyloBits.PNdescriptive			# show() commands for trees etc.
using PhyloBits.TrUtils# basic utility functions, R-like functions				# for prt() tree tables (tree dataframes, trdfs), bd_liks(), etc.
using PhyloBits.TreeTable

#using PhyBEARS.Parsers



###########

# CREATE THE ROOT STATE/SPEED SAVER

	# x1 = ["test", "totallnl", "rootlnl", "Null", "A", "B", "C", "AB", "BC", "AC", "ABC"]
	# headings = permutedims(x1)
	# writedlm("/Users/wbla447/Desktop/School/WallisCastor/September2023/Tester/Speed_Test/results.csv", headings, ',')

	# x2 = ["test", "speed", "totallnl", "rootlnl"]
	# headings_speed = permutedims(x2)
	# writedlm("/Users/wbla447/Desktop/School/WallisCastor/September2023/Tester/Speed_Test/speeds.csv", headings_speed, ',')

	##### Finding Filenames ######
	# Because not all of these sims worked, we have to loop through not numerically, but only the ones taht work!

	# d = readdir()
	# directorynames = d[isdir.(d)]

	cd("/Users/wbla447/Desktop/Files/Tester/M0_tests_Julia")
	readdir()
	file_list = readdir()
	
	headings = ["Example", "N_ranges", "Model", "Time_taken", "Total_lnl", "Bgb_lnl", "Root_lnl"]
	headings = permutedims(headings)
	writedlm(string("/Users/wbla447/Desktop/Files/Tester/Julia_speedtest_bgb.csv"), headings, ',')
	###### RUN INFERENCE

	#### THINGS TO CHANGE ####
    for i in 2:length(file_list)

        cd("/Users/wbla447/Desktop/Files/Tester/M0_tests_Julia")
        cd(file_list[i])
		pwd()
        print(file_list[i])
		lgdata_fn = "geog_data.txt" #geography file
		trfn = "tree.newick"

		### Errors:
		# Frog - prt(tr): node not defined
		# Homs: - 130: BoundsError: attempt to access 8-element Vector{Vector{Float64}} at index [-99999]
		# Podocarp: - 122: BoundsError: attempt to access 2-element interpolate((::Vector{Float64},), ::Vector{Vector{Float64}}, Gridded(Linear())) with element type Vector{Float64} at index [101.4275786173015]

	##### RUN INFERENCE ONCE ABOVE IS CHANGED
		# This simulation has 50 living species

		tr = readTopology(trfn)
		trdf = prt(tr);
		oldest_possible_age = 10000.0

		geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
		include_null_range = true
		numareas = Rncol(geog_df)-1

        if numareas < 10
            continue
        end

		max_range_size = numareas
		n = numstates_from_numareas(numareas, max_range_size, include_null_range)

		# DEC-type SSE model on Hawaiian Psychotria
		# We are setting "j" to 0.0, for now -- so, no jump dispersal
		bmo = construct_BioGeoBEARS_model_object(); 
		#bmo.type[bmo.rownames .== "j"] .= "free";
		bmo.est[bmo.rownames .== "birthRate"] .= ML_yule_birthRate(tr);
		bmo.est[bmo.rownames .== "deathRate"] .= 0.0;
		bmo.est[bmo.rownames .== "d"] .= 0.034;
		bmo.est[bmo.rownames .== "e"] .= 0.028;
		bmo.est[bmo.rownames .== "a"] .= 0.0;
		bmo.est[bmo.rownames .== "j"] .= 0.0;
		bmo.est[bmo.rownames .== "u"] .= 0.0;
		bmo.est[bmo.rownames .== "x"] .= 0.0;
		bmo.est[bmo.rownames .== "xv"] .= 0.0;
		bmo.max[bmo.rownames .== "xv"] .= 0.0;

		bmo_updater_v1!(bmo);


		# Set up the model
		inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=max_range_size, include_null_range=true, bmo=bmo);
		(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

		#######################################################
		# Read in and parse distances and area-of-areas
		#######################################################
		files.times_fn = ""
		files.distances_fn = ""
		files.area_of_areas_fn = ""

		# Construct interpolators, times-at-which-to-interpolate QC
		p = p_Ds_v5;
		interpolators = files_to_interpolators(files, setup.numareas, setup.states_list, setup.v_rows, p.p_indices.Carray_jvals, p.p_indices.Carray_kvals, trdf; oldest_possible_age=100.0);

		p_Es_v12 = (n=p_Ds_v5.n, params=p_Ds_v5.params, p_indices=p_Ds_v5.p_indices, p_TFs=p_Ds_v5.p_TFs, uE=p_Ds_v5.uE, terms=p_Ds_v5.terms, setup=inputs.setup, states_as_areas_lists=inputs.setup.states_list, use_distances=true, bmo=bmo, interpolators=interpolators);

		# Add Q, C interpolators
		p_Es_v12 = p = PhyBEARS.TimeDep.construct_QC_interpolators(p_Es_v12, p_Es_v12.interpolators.times_for_SSE_interpolators);

		# Solve the Es
		prob_Es_v12 = DifferentialEquations.ODEProblem(PhyBEARS.SSEs.parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, Es_tspan, p_Es_v12);
		sol_Es_v12 = solve(prob_Es_v12, solver_options.solver, save_everystep=solver_options.save_everystep, abstol=solver_options.abstol, reltol=solver_options.reltol);

		p = p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

		# Solve the Ds
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = PhyBEARS.TreePass.iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^5, return_lnLs=true)


		#######################################################
		# Maximum likelihood inference
		#######################################################
		bmo.type[bmo.rownames .== "x"] .= "fixed"
		bmo.type[bmo.rownames .== "xv"] .= "fixed"
		bmo.type[bmo.rownames .== "birthRate"] .= "free"
		bmo.type[bmo.rownames .== "deathRate"] .= "fixed"
		pars = bmo.est[bmo.type .== "free"]
		parnames = bmo.rownames[bmo.type .== "free"]
		func = x -> func_to_optimize_v12(x, parnames, inputs, p_Ds_v12; returnval="bgb_lnL", printlevel=1)
		# pars = [0.4, 0.001, -0.001, 4.0, 0.25, 0.0]
		pars = bmo.est[bmo.type .== "free"]
		func(pars)
		function func2(pars, dummy_gradient!)
			return func(pars)
		end # END function func2(pars, dummy_gradient!)


		using NLopt
		opt = NLopt.Opt(:LN_BOBYQA, length(pars))
		ndims(opt)
		opt.algorithm
		algorithm_name(opt::Opt)
		opt.min_objective = func2;
		lower = bmo.min[bmo.type .== "free"];
		upper = bmo.max[bmo.type .== "free"];
		opt.lower_bounds = lower::Union{AbstractVector,Real};
		opt.upper_bounds = upper::Union{AbstractVector,Real};
		opt.ftol_abs = 0.0001 # tolerance on log-likelihood
		#opt.ftol_rel = 0.01 # tolerance on log-likelihood
		opt.xtol_abs = 0.00001 # tolerance on parameters
		#opt.xtol_rel = 0.001 # tolerance on parameters

        print(file_list[i])

		@time (optf,optx,ret) = NLopt.optimize!(opt, pars)
		#######################################################
		# d=0.09035,	e=0.00116,	x=-0.53656,	xv=7.41251,	birthRate=0.25956,	deathRate=0.21085,	Julia_sum_lq=-175.6323, rootstates_lnL=13.6898,	Julia_total_lnLs1=-161.9425, bgb_lnL=-44.7349
		# (-44.73494901514184, [0.09035419093152777, 0.001159414738073889, -0.5365638714798167, 7.41251181819682, 0.2595586392951748, 0.21085253647485827], :ROUNDOFF_LIMITED)
		# d=0.09035,	e=0.00116,	x=-0.53656,	xv=7.41251,	birthRate=0.25956,	deathRate=0.21085,	Julia_sum_lq=-175.6323, rootstates_lnL=13.6898,	Julia_total_lnLs1=-161.9425, bgb_lnL=-44.7349
		# (-44.73494901514184, [0.09035419093152777, 0.001159414738073889, -0.5365638714798167, 7.41251181819682, 0.2595586392951748, 0.21085253647485827], :ROUNDOFF_LIMITED)



		# Get the inputs & res:
		pars = optx;
		inputs.bmo.est[inputs.bmo.type .== "free"] .= pars;
		bmo_updater_v1!(inputs.bmo);
		res = inputs.res;

		# Solution, under best ML parameters
		p_Ds_v5_updater_v1!(p_Ds_v12, inputs);
		p_Es_v12 = TimeDep.construct_QC_interpolators(p_Ds_v12, p_Ds_v12.interpolators.times_for_SSE_interpolators);

		# Solve the Es
		prob_Es_v12 = DifferentialEquations.ODEProblem(parameterized_ClaSSE_Es_v12_simd_sums, p_Es_v12.uE, inputs.Es_tspan, p_Es_v12)
		# This solution is an interpolator
		sol_Es_v12 = solve(prob_Es_v12, inputs.solver_options.solver, save_everystep=inputs.solver_options.save_everystep, abstol=inputs.solver_options.abstol, reltol=inputs.solver_options.reltol);
		p_Ds_v12 = (n=p_Es_v12.n, params=p_Es_v12.params, p_indices=p_Es_v12.p_indices, p_TFs=p_Es_v12.p_TFs, uE=p_Es_v12.uE, terms=p_Es_v12.terms, setup=p_Es_v12.setup, states_as_areas_lists=p_Es_v12.states_as_areas_lists, use_distances=p_Es_v12.use_distances, bmo=p_Es_v12.bmo, interpolators=p_Es_v12.interpolators, sol_Es_v12=sol_Es_v12);

		# Calculate the Ds, and final log-likelihood etc.
		(total_calctime_in_sec, iteration_number, Julia_sum_lq, rootstates_lnL, Julia_total_lnLs1, bgb_lnL) = iterative_downpass_nonparallel_ClaSSE_v12!(res; trdf=trdf, p_Ds_v12=p_Ds_v12, solver_options=inputs.solver_options, max_iterations=10^6, return_lnLs=true)

		example = file_list[i]
		N_ranges = numareas
		Model = "Julia"
		Time_taken = string(total_calctime_in_sec)
		Full_lnl = string(Julia_total_lnLs1)
		Bears_lnl = string(bgb_lnL)
		Root_lnl = string(rootstates_lnL)
		test_speedtable = [example, N_ranges, Model, Time_taken, Full_lnl, Bears_lnl, Root_lnl]
		vert3 = permutedims(test_speedtable)
		CSV.write(string("/Users/wbla447/Desktop/Files/Tester/Julia_speedtest_bgb.csv"), Tables.table(vert3), append=true)

	end

# 		Rnames(res)
# 		rootstateprobs = round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)

# 		# 0.0
# 		#  0.0
# 		#  0.0
# 		#  0.0
# 		#  0.296
# 		#  0.028
# 		#  0.292
# 		#  0.384

#         # results object before uppass
#         rn(res)
#         res.uppass_probs_at_each_nodeIndex_branchBot
#         res.uppass_probs_at_each_nodeIndex_branchTop

#         res.anc_estimates_at_each_nodeIndex_branchBot
#         res.anc_estimates_at_each_nodeIndex_branchTop
  
#   		# Note that print(b) fails since b doesn't exist

# 		uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=false)
#         # results object after uppass
#         rn(res)
#         res.uppass_probs_at_each_nodeIndex_branchBot
#         res.uppass_probs_at_each_nodeIndex_branchTop

#         res.anc_estimates_at_each_nodeIndex_branchBot
#         res.anc_estimates_at_each_nodeIndex_branchTop


#         # The res tables are hard to read, the vfft function helps
#         # by putting into DataFrame
#         vfft(res.anc_estimates_at_each_nodeIndex_branchBot)
#         vfft(res.anc_estimates_at_each_nodeIndex_branchTop)

#         # View the ancestral range probabilities in R's default node order
#         R_node_order = sort(trdf, :Rnodenums).nodeIndex

#         # Tree table is reordered
#         trdf[R_node_order,:]
#         trdf

#         # Ancestral state probs are reorderd
#         vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_node_order])
#         vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order])

# 		print(vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order]))

# 		r_nodes_table2 = trdf.Rnodenums[R_node_order,:]
# 		r_nodes_table2
		
# 		brlen2 = trdf.brlen[R_node_order,:]
# 		nodename = trdf[R_node_order,:nodeName]
# 		node_age = trdf.node_age[R_node_order,:]
		
# 		testbottom2 = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_node_order])
# 		testtop2 = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order])
		
# 		testbottom2[!, "R_nodes"] .= r_nodes_table2
# 		testbottom2[!, "brlen"] .= brlen2
# 		testbottom2[!, "node_age"] .= node_age
# 		testbottom2[!, "nodename"] .= nodename
# 		testbottom2
		
# 		testtop2[!, "R_nodes"] .= r_nodes_table2
# 		testtop2[!, "brlen"] .= brlen2
# 		testtop2[!, "node_age"] .= node_age
# 		testtop2[!, "nodename"] .= nodename
# 		testtop2

# 		st_rootstateprobs = string.(rootstateprobs)

# 		st_lnl = string(Julia_total_lnLs1)
# 		st_rslnl = string(rootstates_lnL)


# 		example = file_list[i]
# 		N_ranges = numareas
# 		Model = "Julia"
# 		Time_taken = string(total_calctime_in_sec)
# 		test_speedtable = [example, N_ranges, Model, Time_taken]
# 		vert3 = permutedims(test_speedtable)
# 		CSV.write(string("/Users/wallis/Documents/GitHub/Files/Tester/Julia_speedtestpt2.csv"), Tables.table(vert3), append=true)
# 	end
# # # Install modified "castor" package in R
# # # install.packages(pkgs="/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz", lib="/Library/Frameworks/R.framework/Resources/library/", repos=NULL, type="source")

# # # Write model out to text files that can be read in to simulator
# # timepoints = seq(0.0, 10.0, 1.0)
# # # (the best way to do this is to do simulations for a fixed period of time; the number of taxa
# # #  will vary, but have an average)
# # outfns = model_to_text_v12(p_Ds_v12, timepoints; prefix="")