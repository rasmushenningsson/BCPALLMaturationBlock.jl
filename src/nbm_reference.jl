function nbm_reference_sample_names()
	# NB: REMOVED NBM-5-CD34 that seems have poor quality
	["NBM-1-MNC", "NBM-2-MNC", "NBM-3-MNC", "NBM-3-CD34", "NBM-4-CD34", "NBM-3-CD19", "NBM-4-CD19", "NBM-5-CD19"]
end

function nbm_reference_raw_counts(; ref_names=nothing, kwargs...)
	ref_names = @something ref_names nbm_reference_sample_names()
	load_samples(ref_names; kwargs...)
end


function nbm_reference_replacements(proj_names, ref_names=nothing; kwargs...)
	ref_names = @something ref_names nbm_reference_sample_names()

	ref = nbm_reference_raw_counts(; kwargs..., ref_names)
	proj = nbm_reference_raw_counts(; kwargs..., ref_names=proj_names)

	ref=>proj
end



function nbm_reference_counts(data=nothing; kwargs...)
	data = @something data nbm_reference_raw_counts(; kwargs...)

	annots = load_annotations("ALL_merged_metadata_221027.csv.gz")
	data = SCP.annotate_obs(data, annots)

	data = SCP.filter_obs("percent.mt" => !ismissing, data) # Start with same filtering as in R analysis
	data = SCP.filter_obs("nCount_RNA" => >=(2000), data) # But make it stricter

	# TODO: keep only b-cell lineage

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
	                            seed = @something(seed,6712),
	                            k = 100,
	                            niter = 400,
	                            k_projection = 10)

	# TODO: Rotate such that HSC is at the top of the plot

	fl
end
