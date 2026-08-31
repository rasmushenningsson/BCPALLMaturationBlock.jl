function scatter_2d!(ax, job; color = :black)
	matrix = fetch!(SCP.get_matrix(job))
	scatter!(ax, matrix; color, markersize = 4)
	ax
end




function category_colors(annot, colors::Vector{<:Pair})
	unique_annotations = unique(annot)
	unique_annotations_set = Set(unique_annotations)
	colors = filter(x->x[1] in unique_annotations_set, colors)

	categories = first.(colors)
	@assert isempty(setdiff(unique_annotations, categories)) # ensure all categories have colors specified

	# categories, colors
	colors
end

function category_colors(annot, ::Nothing)
	categories = unique(annot)
	colors = distinguishable_colors(length(categories)+2, [colorant"black", colorant"white"])[3:end]

	categories .=> colors
end



function scatter_categorical_2d!(ax, job, annot_name; colors=nothing)
	matrix = fetch!(SCP.get_matrix(job))
	annot = fetch!(SCP.value_column_data(SCP.annotation(SCP.get_obs(job), annot_name)))

	colors = category_colors(annot, colors)
	plots = [scatter!(ax, matrix[:,isequal.(annot, cat)]; markersize=6, color, label=cat) for (cat,color) in colors]

	axislegend(ax, plots .=> Ref((;markersize=16)), first.(colors)) # use a larger markersize in the legend
	ax
end


function draw_trajectory!(ax, trajectory::Trajectory; nticks=0)
	color = colorant"#C0C0C0"
	tick_color = colorant"#808080"

	poly = polygon(trajectory)

	outer_edge = poly[:,1:cld(end,2)-1]
	lines!(ax, outer_edge; color, linewidth=4)

	inner = poly[:,cld(end,2)-1:end]
	lines!(ax, inner; color, linewidth=4, linestyle=:dash)


	arrow_end = Point2f(outer_edge[:,end])
	arrow_dir = Vec2f(arrow_end - Point2f(outer_edge[:,end-1]))
	arrow_dir /= norm(arrow_dir)
	arrow_dir *= 300 # make more room for arrow head
	arrow_start = arrow_end - arrow_dir

	arrows2d!(ax, arrow_start, arrow_dir; tipwidth=10, tiplength=10, tailwidth=0, taillength=0, minshaftlength=0, maxshaftlength=0, color=[color,color])


	if nticks != 0
		timepoints = range(0,1,size(outer_edge,2))

		for i in 1:nticks
			t = (i-1)/(nticks-1)

			j = searchsortedlast(timepoints, t)
			j1 = max(j, 1)
			j2 = min(j+1, size(outer_edge,2))
			α = (t - timepoints[j1]) / step(timepoints)

			tick_pos = Point2f(outer_edge[:,j1].*(1-α) .+ outer_edge[:,j2].*α)

			text!(ax, tick_pos; text=rich("t", subscript(string(t))), fontsize=40,  align=(:center,:center), color=tick_color)
		end
	end


	ax
end




function padded_kde(x; boundary::Tuple{<:Real,<:Real}, npoints::Integer,
                       bandwidth::Real, pad_factor::Real = 6)
    lo, hi = boundary
    dx_target = (hi - lo) / (npoints - 1)
    pad = pad_factor * bandwidth

    # oversample the padded interval so resolution on [lo, hi] is at least
    # as fine as the requested npoints, after we crop/interpolate back down
    npoints_padded = ceil(Int, (hi - lo + 2pad) / dx_target)

    kd = kde(x; boundary = (lo - pad, hi + pad),
             npoints = npoints_padded, bandwidth = bandwidth)

    xrange = range(lo, hi, length = npoints)
    itp = InterpKDE(kd)
    density = [KernelDensity.pdf(itp, t) for t in xrange]

    return xrange, density
end




function plot_trajectory_histogram(ax, job, annot_name=nothing; trajectory::Trajectory, σ=1e-2, colors=nothing)
	matrix = fetch!(SCP.get_matrix(job))
	t = trajectoryprojection(matrix, trajectory) # time value for each cell - or nothing if outside the trajectory

	# remove cells that does not fall inside the trajectory
	mask = t .!== nothing
	t = identity.(t[mask]) # drop nothings

	kde_kwargs = (; boundary=(0.0, 1.0), npoints = 2048, bandwidth=σ)

	if annot_name === nothing
		# gaussian smoothing
		x, y = padded_kde(t; kde_kwargs...)
		band!(ax, x, zeros(length(x)), y)
	else # categories
		annot = fetch!(SCP.value_column_data(SCP.annotation(SCP.get_obs(job), annot_name)))
		annot = annot[mask,:]

		colors = category_colors(annot, colors)

		ylower = nothing
		for (cat,color) in colors
			m = isequal.(annot, cat)
			tc = t[m]

			# gaussian smoothing
			xc, yc = padded_kde(tc; kde_kwargs...)
			yc .*= (length(tc)/length(t)) # rescale to have total area proportional to number of points in category

			ylower = @something ylower zeros(length(xc))
			yupper = ylower .+ yc
			band!(ax, xc, ylower, yupper; color)

			ylower = yupper
		end
	end

	ax
end
