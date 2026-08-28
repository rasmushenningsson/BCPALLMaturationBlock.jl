# It could be argued that this should detect changes in the table on disk, and reload the file.
# But that will be very rare in practice. So it's OK to restart in that case.
const sample_table = OncePerProcess{DataFrame}() do
	CSV.read(joinpath(pkgdir(BCPALLMaturationBlock), "data/sample_table.tsv"), DataFrame)
end


"""
	find_samples(; [key=>value...])

Get names of samples. Optionally specify onr or more `key=>value` pairs to restrict the search.

* `key` can be one of $(join(string.('`', names(sample_table()), '`'), ", ", " and ")).
* `value` can be a:
    * `String` - samples matching the value are returned.
    * `Regex` - samples matching the regex pattern are returned.
    * predicate function - samples for which the predicate returns true are returned.

# Examples

Retrieve names of all samples:
```julia
julia> find_samples()
```

Retrieve names of DUX4 samples:
```julia
julia> find_samples(subtype="DUX4")
```

"""
find_samples(; kwargs...) = _find_samples(; kwargs...)
function _find_samples(; kwargs...)
	si = sample_table()
	for (k,v) in kwargs
		if v isa String
			si = filter(k=>isequal(v), si)
		elseif v isa Regex
			si = filter(k=>contains(v), si)
		else
			si = filter(k=>v, si)
		end
	end
	convert.(String, si.sample_name)
end

function get_sample_info(sample_name::String, col)
	si = sample_table()
	ind = findall(==(sample_name), si.sample_name)
	si[only(ind), col]
end
function get_sample_table(sample_names::AbstractVector, col)
	si = sample_table()
	ind = indexin(sample_names, si.sample_name)
	si[ind, col]
end
get_sample_table(sample_name::String, col) = get_sample_table([sample_name], col)
get_sample_table(col) = sample_table()[:, col]


function sample2path(sample_name; bcpall_path=get_bcpall_path())
	si = sample_table()
	i = findfirst(isequal(sample_name), si.sample_name)
	i === nothing && error("sample_name \"$sample_name\" not found.")

	bcpall_path === nothing && error("Use `set_bcpall_path(\"/my/bcpall/path/\")` or use `bcpall_path` kwarg.")
	joinpath(bcpall_path, string(sample_name, ".h5"))
end


function load_samples(sample_names; kwargs...)
	path_kwargs, load_kwargs = split_kwargs(endswith("_path"); kwargs...)

	sample_paths = sample2path.(sample_names; path_kwargs...)
	SCP.load_counts(sample_paths; sample_names, load_kwargs...)
end

function load_annotations(filename; annotations_path=get_annotations_path())
	SCP.load_csv(joinpath(annotations_path, filename))
end
