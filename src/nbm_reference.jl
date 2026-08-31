function nbm_reference_sample_names()
	# NB: REMOVED NBM-5-CD34 that seems have poor quality
	["NBM-1-MNC", "NBM-2-MNC", "NBM-3-MNC", "NBM-3-CD34", "NBM-4-CD34", "NBM-3-CD19", "NBM-4-CD19", "NBM-5-CD19"]
end

function nbm_reference_raw_counts(; ref_names=nothing, kwargs...)
	ref_names = @something ref_names nbm_reference_sample_names()
	load_samples(ref_names; kwargs...)
end


function nbm_reference_counts(data=nothing; kwargs...)
	data = @something data nbm_reference_raw_counts(; kwargs...)
	annots = load_annotations("ALL_merged_metadata_221027.csv.gz")
	data = SCP.annotate_obs(data, annots)

	# B-lineage only for NBM
	b_cell_lineage = Set(("HSC","LMPP","CLP","Pro-B","Pre-B","Immature B","Naive B","Memory B","Plasmablast"))
	data = SCP.filter_obs("annotations_B_new_broad" => in(b_cell_lineage), data; project_obs_ids=:skip)

	data
end




function nbm_reference_filtered_counts(data=nothing; kwargs...)
	data = @something data nbm_reference_counts(; kwargs...)

	data = SCP.filter_obs("percent.mt" => !ismissing, data) # Start with same filtering as in R analysis
	data = SCP.filter_obs("nCount_RNA" => >=(2000), data) # But make it stricter

	# TODO: Blast only for ALL

	data
end



function nbm_reference_transformed(data=nothing; var_filter_transformed=nothing, kwargs...)
	data = @something data nbm_reference_filtered_counts(; kwargs...)
	data = SCP.sctransform(data)

	if var_filter_transformed !== nothing
		data = SCP.filter_var(var_filter_transformed, data)
	end

	data
end

function nbm_reference_normalized(data=nothing; kwargs...)
	data = @something data nbm_reference_transformed(; kwargs...)
	SCP.normalize_matrix(data, "percent.mt")
end

function nbm_reference_pca(data=nothing; nsv=40, seed=1234, ref_names=nothing, kwargs...)
	data = @something data nbm_reference_normalized(; ref_names, kwargs...)
	SCP.pca(data; nsv, seed) # TODO: Pass some kwargs here?
end

function nbm_reference_force_layout(data=nothing; ndim, ref_names=nothing, seed=nothing, kwargs...)
	data = @something data nbm_reference_pca(; ref_names, kwargs...)

	fl = SCP.force_layout(data; ndim,
	                            seed = @something(seed,6771),
	                            k = 100,
	                            niter = 400,
	                            k_projection = 10)

	# TODO: Rotate such that HSC is at the top of the plot
	annot = "annotations_B_new_broad"
	if ndim==2
		# transform = SCP.find_optimal_coord_transform(fl, annot=>isequal("HSC"), annot=>isequal("T-cells"))
		transform = SCP.find_optimal_coord_transform(fl, annot=>isequal("HSC"), annot=>isequal("Memory B"))
		# match rotation from paper
		transform = fetch!(transform)
		transform = rot2d(pi/9)*transform
	else
		transform = SCP.find_optimal_coord_transform(fl, annot=>isequal("HSC"), annot=>isequal("Memory B"), annot=>isequal("Plasmablast"))
	end

	fl = SCP.transform_coords(fl, transform; keep_var=true)

	fl
end


# Trajectories - NB: These are specific to a given force layout realization

# From paper
# bcell_trajectory_HSC2MemoryB() = BezierTrajectory([-800 -1800 -1200   300 1800;
#                                                    2600   400 -1200 -1200 -300], 2300, 101)

# Slightly adjusted since force layout had breaking changes between paper and current version
bcell_trajectory_HSC2MemoryB() = BezierTrajectory([-800 -1800 -1200   300 1800;
                                                   2600   400 -1200 -1200 -300], 2500, 101)


bcell_trajectory_bridge() = LineTrajectory([-1600 600; 500 -600], 1000)


function nbm_reference_plot(; seed=nothing)
	fl = nbm_reference_force_layout(; ndim=2, seed)
	fig = scatter_categorical_2d(fl, "annotations_B_new_broad"; colors=colortable_annotations_B_new_broad(), )

	ax = current_axis(fig)
	draw_trajectory!(ax, bcell_trajectory_HSC2MemoryB(); nticks=6)
	# draw_trajectory!(ax, bcell_trajectory_bridge(); nticks=6)

	fig
end
