function nbm_reference_sample_names()
	# NB: REMOVED NBM-5-CD34 that seems have poor quality
	["NBM-1-MNC", "NBM-2-MNC", "NBM-3-MNC", "NBM-3-CD34", "NBM-4-CD34", "NBM-3-CD19", "NBM-4-CD19", "NBM-5-CD19"]
end

function nbm_reference_raw_counts(; ref_names=nothing, kwargs...)
	ref_names = @something ref_names nbm_reference_sample_names()
	load_samples(ref_names; kwargs...)
end


function nbm_reference_replacements(sample_names; ref_names=nothing, kwargs...)
	ref_names = @something ref_names nbm_reference_sample_names()

	ref = nbm_reference_raw_counts(; kwargs..., ref_names)
	proj = nbm_reference_raw_counts(; kwargs..., ref_names=sample_names)

	ref=>proj
end




function nbm_reference_counts(data=nothing; kwargs...)
	data = @something data nbm_reference_raw_counts(; kwargs...)
	annots = load_annotations("cell_annotations.csv.gz")
	data = SCP.annotate_obs(data, annots)

	data = SCP.filter_obs("percent.mt" => !ismissing, data) # Start with same filtering as in R analysis
	data = SCP.filter_obs("nCount_RNA" => >=(2000), data) # But make it stricter

	# B-lineage only for NBM
	b_cell_lineage = Set(("HSC","LMPP","CLP","Pro-B","Pre-B","Immature B","Naive B","Memory B","Plasmablast"))
	data = SCP.filter_obs("celltypes" => in(b_cell_lineage), data; project_obs_ids=:skip) # NB: :skip means this filter is ignored during projection

	# Keep NBM cells, but for (projected) ALL cases keep Blast cells
	data = SCP.filter_obs("tumor_celltype"=>in(("Blast","NBM")), data)

	data
end





function nbm_reference_transformed(data=nothing; var_filter_transformed=nothing, kwargs...)
	data = @something data nbm_reference_counts(; kwargs...)
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
	annot = "celltypes"
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

# Slightly adjusted (wider trajectory) since force layout had breaking changes between paper and current SCP version
bcell_trajectory_HSC2MemoryB() = BezierTrajectory([-800 -1800 -1200   300 1800;
                                                   2600   400 -1200 -1200 -300], 2500, 101)


bcell_trajectory_bridge() = LineTrajectory([-1600 600; 500 -600], 1000)


function nbm_reference_plot(; seed=nothing)
	fig = Figure(; size=(1024, 1024))
	ax = Axis(fig[1, 1]; autolimitaspect=true)

	hidedecorations!(ax)
	hidespines!(ax)


	fl = nbm_reference_force_layout(; ndim=2, seed)
	scatter_categorical_2d!(ax, fl, "celltypes"; colors=colortable_celltypes())
	draw_trajectory!(ax, bcell_trajectory_HSC2MemoryB(); nticks=6)
	# draw_trajectory!(ax, bcell_trajectory_bridge(); nticks=6)

	fig
end



function nbm_reference_trajectory_histogram()
	fig = Figure(; size=(1024, 256))
	ax = Axis(fig[1, 1]; xlabel="Inferred time points")

	hideydecorations!(ax)

	fl = nbm_reference_force_layout(; ndim=2)

	# plot_trajectory_histogram(ax, fl, nothing; trajectory=bcell_trajectory_HSC2MemoryB(), σ=1e-2)
	plot_trajectory_histogram(ax, fl, "celltypes"; trajectory=bcell_trajectory_HSC2MemoryB(), σ=1e-2, colors=colortable_celltypes())

	fig
end
