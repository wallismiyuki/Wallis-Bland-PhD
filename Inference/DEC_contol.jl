
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
using PhyBEARS.TimeDep
using PhyBEARS.Uppass

using PhyloBits.PNtypes						# Types: Tree, Node, Edge etc.
using PhyloBits.PNmanipulateNet
using PhyloBits.PNreadwrite				# Helper functions
using PhyloBits.PNreadwrite 			# Reading and writing trees; readTopology
using PhyloBits.PNdescriptive			# show() commands for trees etc.
using PhyloBits.TrUtils						# basic utility functions, R-like functions				# for prt() tree tables (tree dataframes, trdfs), bd_liks(), etc.


#using PhyBEARS.Parsers

# Which folder is this?
version = "V1/"

for group in 2:8

	# Change the working directory as needed
	wd = string("/Users/wbla447/Library/CloudStorage/OneDrive-TheUniversityofAuckland/School/WallisCastor/September2023/split_only/", version, "ss", group)
	cd(wd)



	x = ["Sim", "Null", "A", "B", "C", "AB", "BC", "AC", "ABC", "totallnl", "rootlnl"]
	headings = permutedims(x)
	writedlm(string("/Users/wbla447/Library/CloudStorage/OneDrive-TheUniversityofAuckland/School/WallisCastor/September2023/split_only/", version, "bigrootstateprobs", group, "_control_7000.csv"), headings, ',')

	### LOAD THIS OUT OF THE LOOP? ###

	# Times_3 = repeat([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],3)
	# Land1_3 = repeat(["test1","test1","test2"],inner = 11)
	# Land2_3 = repeat(["test2","test3","test3"],inner = 11)
	# Distances_km_3 = [0, 1000, 2000, 3500, 4100, 4900, 6000, 5800, 7000, 7500, 8000,0, 2000, 3000, 4000, 5100, 5900, 7000, 6800, 8000, 8500, 9000,0, 500, 1000, 2500, 3100, 3900, 5000, 4800, 6000, 6500, 7000] 
	# testdf = DataFrame(Times = repeat([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],3), Land1 = repeat(["test1","test1","test2"],inner = 11), Land2 = repeat(["test2","test3","test3"],inner = 11), Distances_km = [0, 1000, 2000, 3500, 4100, 4900, 6000, 5800, 7000, 7500, 8000,0, 2000, 3000, 4000, 5100, 5900, 7000, 6800, 8000, 8500, 9000,0, 500, 1000, 2500, 3100, 3900, 5000, 4800, 6000, 6500, 7000])


	##### Finding Filenames ######
	# Because not all of these sims worked, we have to loop through not numerically, but only the ones taht work!

	d = readdir()
	directorynames = d[isdir.(d)]



	#### Change Wallis Distance style to Nick's style ######

	# numareas = 3
	# unique_times = 1.0 .* sort(unique(testdf.Times)) # unique times must be declared every time, as we later push the endtime onto the vector
	# distmats = [Array{Float64}(undef, numareas, numareas) for _ = 1:(length(unique_times)+1)]

	# # Fill in the distmats
	# for i in 1:(length(unique_times))
	# 	distmats[i] .= 0.0  # zero-out the distances matrix
	# 	TF = testdf.Times .== unique_times[i]
	# 	tmpdf = testdf[TF,:]
	# 	for j in 1:Rnrow(tmpdf)
	# 		land1 = tmpdf.Land1[j]
	# 		land2 = tmpdf.Land2[j]
	# 		area1_num = extract_first_integer_from_string(land1)
	# 		area2_num = extract_first_integer_from_string(land2)
	# 		distmats[i][area1_num, area2_num] = convert(Float64, testdf.Distances_km[i])
	# 		distmats[i][area2_num, area1_num] = convert(Float64, testdf.Distances_km[i])
	# 	end
	# end


	# # Add one last distmat on the end to cover the bottom of the tree
	# distmats[length(distmats)] .= distmats[length(distmats)-1]
	# push!(unique_times, oldest_possible_age)
	# sort!(unique_times)
	# times = unique_times

	# show(distmats)

	# # Let's divide the distances by the minimum nonzero distance
	# dists_vector = collect(Iterators.flatten(vec.(distmats)))
	# dists_vector2 = dists_vector[dists_vector .> 0.0]
	# maxval = maximum(dists_vector2)

	# for i in 1:length(distmats)
	# 	distmats[i] .= distmats[i] ./ maxval
	# end
	# distmats

	# writedlm("/Users/wbla447/Desktop/School/WallisCastor/8StateFiles/distances.txt", distmats)


	####### CREATE TREE.NEWICK

	# for i in 1 : length(directorynames)

	# 	cd(string("/Users/wbla447/Desktop/School/WallisCastor/8SplitSpreadFiles/GoodSims/GoodSims", group, "_095/", directorynames[i]))
	# 	trfn = "newtree"
	# 	tr = PhyloBits.PNreadwrite.readTopology(trfn)

	# 	# Change the tip names from e.g. "6" to "sp6"
	# 	newnames = deepcopy(tr.names)
	# 	for i in 1:length(newnames)
	# 		newnames[i] = Rpaste0(["sp", tr.names[i]])
	# 	end
	# 	tr.names .= newnames

	# 	# Change the names on the nodes
	# 	for i in 1:length(tr.node)
	# 		if tr.node[i].leaf == true
	# 			tr.node[i].name = Rpaste0(["sp", tr.node[i].name])
	# 		end
	# 	end

	# 	PhyloBits.PNreadwrite.writeTopology(tr, "tree.newick")
	# 	cd(string("/Users/wbla447/Desktop/School/WallisCastor/8SplitSpreadFiles/GoodSims/GoodSims", group, "_095/"))

	# end


	###### RUN INFERENCE

	for i in 1:length(directorynames)
		cd(string("/Users/wbla447/Library/CloudStorage/OneDrive-TheUniversityofAuckland/School/WallisCastor/September2023/split_only/", version ,"ss", group, "/", directorynames[i]))

		
		# This simulation has 50 living species
		trfn = "living_tree_noNodeLabels.newick"
		tr = readTopology(trfn)
		trdf = prt(tr);
		oldest_possible_age = 100.0

		lgdata_fn = "geog_living.data"
		geog_df = Parsers.getranges_from_LagrangePHYLIP(lgdata_fn);
		include_null_range = true
		numareas = Rncol(geog_df)-1
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
		bmo.max[bmo.rownames .== "xv"] .= 10.0;

		bmo_updater_v1!(bmo);


		# Set up the model
		inputs = PhyBEARS.ModelLikes.setup_DEC_SSE2(numareas, tr, geog_df; root_age_mult=1.5, max_range_size=NaN, include_null_range=true, bmo=bmo);
		(setup, res, trdf, bmo, files, solver_options, p_Ds_v5, Es_tspan) = inputs;

		#######################################################
		# Read in and parse distances and area-of-areas
		#######################################################
		files.times_fn = "/Users/wbla447/Library/CloudStorage/OneDrive-TheUniversityofAuckland/School/WallisCastor/April_2023/times.txt"
		files.distances_fn = "/Users/wbla447/Library/CloudStorage/OneDrive-TheUniversityofAuckland/School/WallisCastor/September2023/distances_div_7000.txt"
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
		(optf,optx,ret) = NLopt.optimize!(opt, pars)
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

		Rnames(res)
		rootstateprobs = round.(res.normlikes_at_each_nodeIndex_branchTop[tr.root]; digits=3)

		# 0.0
		#  0.0
		#  0.0
		#  0.0
		#  0.296
		#  0.028
		#  0.292
		#  0.384

        # results object before uppass
        rn(res)
        res.uppass_probs_at_each_nodeIndex_branchBot
        res.uppass_probs_at_each_nodeIndex_branchTop

        res.anc_estimates_at_each_nodeIndex_branchBot
        res.anc_estimates_at_each_nodeIndex_branchTop
  
  		# Note that print(b) fails since b doesn't exist

		try 
            uppass_ancstates_v12!(res, trdf, p_Ds_v12, solver_options; use_Cijk_rates_t=false)
        catch
            continue
        end

        # results object after uppass
        rn(res)
        res.uppass_probs_at_each_nodeIndex_branchBot
        res.uppass_probs_at_each_nodeIndex_branchTop

        res.anc_estimates_at_each_nodeIndex_branchBot
        res.anc_estimates_at_each_nodeIndex_branchTop


        # The res tables are hard to read, the vfft function helps
        # by putting into DataFrame
        vfft(res.anc_estimates_at_each_nodeIndex_branchBot)
        vfft(res.anc_estimates_at_each_nodeIndex_branchTop)

        # View the ancestral range probabilities in R's default node order
        R_node_order = sort(trdf, :Rnodenums).nodeIndex

        # Tree table is reordered
        trdf[R_node_order,:]
        trdf

        # Ancestral state probs are reorderd
        vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_node_order])
        vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order])

		print(vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order]))


		r_nodes_table2 = trdf.Rnodenums[R_node_order,:]
		r_nodes_table2
		
		brlen2 = trdf.brlen[R_node_order,:]
		nodename = trdf[R_node_order,:nodeName]
		node_age = trdf.node_age[R_node_order,:]
		
		testbottom2 = vfft(res.anc_estimates_at_each_nodeIndex_branchBot[R_node_order])
		testtop2 = vfft(res.anc_estimates_at_each_nodeIndex_branchTop[R_node_order])
		
		testbottom2[!, "R_nodes"] .= r_nodes_table2
		testbottom2[!, "brlen"] .= brlen2
		testbottom2[!, "node_age"] .= node_age
		testbottom2[!, "nodename"] .= nodename
		testbottom2
		
		testtop2[!, "R_nodes"] .= r_nodes_table2
		testtop2[!, "brlen"] .= brlen2
		testtop2[!, "node_age"] .= node_age
		testtop2[!, "nodename"] .= nodename
		testtop2
		
		CSV.write("branchbottoms_con_7000.csv", testbottom2)
		CSV.write("branchtops_con_7000.csv", testtop2)
		CSV.write("inferredtable_con_7000.csv", trdf[R_node_order,:])
		
		writedlm("rootstate_probs_con_7000.txt", rootstateprobs)

		st_rootstateprobs = string.(rootstateprobs)
		rootstateprobs_tableprint = pushfirst!(st_rootstateprobs,lpad(i,3,"0"))


		st_lnl = string(Julia_total_lnLs1)
		st_rslnl = string(rootstates_lnL)
		rootstateprobslnl_table = push!(rootstateprobs_tableprint, st_lnl)
		RSP_Tlnl_RSlnl = push!(rootstateprobslnl_table, st_rslnl)
		vert2 = permutedims(RSP_Tlnl_RSlnl)
		CSV.write(string("/Users/wbla447/Library/CloudStorage/OneDrive-TheUniversityofAuckland/School/WallisCastor/September2023/split_only/", version ,"bigrootstateprobs", group ,"_control_7000.csv"), Tables.table(vert2), append=true)

		writedlm("total_lnl_con_7000.txt", Julia_total_lnLs1)
		writedlm("rootstates_lnl_con_7000.txt", rootstates_lnL)
	end
end
# # Install modified "castor" package in R
# # install.packages(pkgs="/GitHub/PhyBEARS.jl/simulator/castor_1.7.2.000004.tar.gz", lib="/Library/Frameworks/R.framework/Resources/library/", repos=NULL, type="source")

# # Write model out to text files that can be read in to simulator
# timepoints = seq(0.0, 10.0, 1.0)
# # (the best way to do this is to do simulations for a fixed period of time; the number of taxa
# #  will vary, but have an average)
# outfns = model_to_text_v12(p_Ds_v12, timepoints; prefix="")



  