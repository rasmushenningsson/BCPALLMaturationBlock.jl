function projection_counts(sample_names; kwargs...)
	data = nbm_reference_raw_counts(; ref_names=sample_names, kwargs...)
	annots = load_annotations("ALL_merged_metadata_221027.csv.gz")
	data = SCP.annotate_obs(data, annots)

	# Blast only for ALL
	missing2na = Base.Fix2(coalesce, "NA")
	blast_pattern = startswith(r"HeH|DUX4|ETV6-RUNX1|BCR-ABL1")
	data = SCP.filter_obs("annotations_B_new_broad_ALLcase" => blast_pattern∘missing2na, data)

	data
end


function projection_replacements(sample_names; ref_names=nothing, kwargs...)
	ref_names = @something ref_names nbm_reference_sample_names()

	ref = nbm_reference_counts(; kwargs..., ref_names)
	proj = projection_counts(sample_names; kwargs...)

	ref=>proj
end


function projection_scatter_plot(sample_names)
	replacements = projection_replacements(sample_names)

	fl = nbm_reference_force_layout(; ndim=2)
	fl_proj = SCP.project(fl, replacements)


	fig = scatter_categorical_2d(fl_proj, "annotations_B_new_broad"; colors=colortable_annotations_B_new_broad(), )

	ax = current_axis(fig)
	draw_trajectory!(ax, bcell_trajectory_HSC2MemoryB(); nticks=6)
	# draw_trajectory!(ax, bcell_trajectory_bridge(); nticks=6)

	fig
end
