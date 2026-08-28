module BCPALLMaturationBlock

export
	set_annotations_path,
	get_annotations_path,
	set_bcpall_path,
	get_bcpall_path,
	fetch! # Reexport from ReproducibleJobs

import SingleCellProjections as SCP
using SingleCellProjections
using ReproducibleJobs
using Preferences
using DataFrames
using CSV
using LinearAlgebra
using Statistics: mean, var


set_annotations_path(path::String) = @set_preferences!("annotations_path"=>expanduser(path))
get_annotations_path() = @load_preference("annotations_path")

set_bcpall_path(path::String) = @set_preferences!("bcpall_path"=>expanduser(path))
get_bcpall_path() = @load_preference("bcpall_path")

end
