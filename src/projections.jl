function projection_scatter_plot(sample_names)
	fig = Figure(; size=(1024, 1024))
	ax = Axis(fig[1, 1]; autolimitaspect=true)
	hidedecorations!(ax)
	hidespines!(ax)

	fl = nbm_reference_force_layout(; ndim=2)
	replacements = nbm_reference_replacements(sample_names)
	fl_proj = SCP.project(fl, replacements)


	scatter_2d!(ax, fl; color=colorant"#BFCCE6")
	scatter_categorical_2d!(ax, fl_proj, "celltypes_original"; colors=colortable_celltypes())

	draw_trajectory!(ax, bcell_trajectory_HSC2MemoryB(); nticks=6)
	# draw_trajectory!(ax, bcell_trajectory_bridge(); nticks=6)

	fig
end


function projection_trajectory_histogram(sample_names)
	fig = Figure(; size=(1024, 256))
	ax = Axis(fig[1, 1])

	hideydecorations!(ax)

	fl = nbm_reference_force_layout(; ndim=2)
	replacements = nbm_reference_replacements(sample_names)
	fl_proj = SCP.project(fl, replacements)

	# plot_trajectory_histogram(ax, fl_proj, nothing; trajectory=bcell_trajectory_HSC2MemoryB(), σ=1e-2)
	plot_trajectory_histogram(ax, fl_proj, "celltypes_original"; trajectory=bcell_trajectory_HSC2MemoryB(), σ=1e-2, colors=colortable_celltypes())

	fig
end
