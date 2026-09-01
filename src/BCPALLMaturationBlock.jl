module BCPALLMaturationBlock

export
	fetch!, # Reexport from ReproducibleJobs
	set_annotations_path,
	get_annotations_path,
	set_bcpall_path,
	get_bcpall_path,
	find_samples,
	load_samples,
	load_annotations,
	nbm_reference_raw_counts,
	nbm_reference_counts,
	nbm_reference_transformed,
	nbm_reference_normalized,
	nbm_reference_pca,
	nbm_reference_force_layout,
	nbm_reference_plot,
	nbm_reference_trajectory_histogram,
	projection_counts,
	projection_replacements,
	projection_scatter_plot,
	projection_trajectory_histogram

import SingleCellProjections as SCP
using SingleCellProjections
using ReproducibleJobs
using Preferences
using DataFrames
using CSV
using LinearAlgebra
using Statistics: mean, var

# basic plotting
using Colors
using GLMakie

using KernelDensity


set_annotations_path(path::String) = @set_preferences!("annotations_path"=>expanduser(path))
get_annotations_path() = @load_preference("annotations_path")

set_bcpall_path(path::String) = @set_preferences!("bcpall_path"=>expanduser(path))
get_bcpall_path() = @load_preference("bcpall_path")


include("utils.jl")
include("samples.jl")

include("trajectories.jl")

include("plot_utils.jl")
include("colors.jl")

include("nbm_reference.jl")
include("projections.jl")

end
