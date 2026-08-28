function split_kwargs(f; kwargs...)
	a_keys = filter(x->f(string(x)), keys(kwargs))
	a = NamedTuple{a_keys}(values(kwargs))
	b = Base.structdiff(values(kwargs), NamedTuple{a_keys})
	a,b
end
