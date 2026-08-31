abstract type Trajectory end

struct LineTrajectory <: Trajectory
	points::Matrix{Float64} # dims × 2 matrix (start and end point)
	width::Float64
end

struct BezierTrajectory <: Trajectory
	# setup
	control::Matrix{Float64} # dims × n matrix (n control points)
	width::Float64

	# discretization
	points::Matrix{Float64} # dims × nparts
	normals::Matrix{Float64} # dims × nparts (scaled by width/2)
end

function BezierTrajectory(control::Matrix, width, nparts::Integer)
	@assert size(control,1)==2
	n = size(control,2)

	t = range(0,1;length=nparts)
	rt = range(1,0;length=nparts)

	# interpolation matrix: n×nparts
	B = zeros(n,nparts)
	dB = zeros(n,nparts) # derivative
	for i in 1:n
		C = binomial(n-1,i-1)
		B[i,:] .= C .* rt.^(n-i) .* t.^(i-1)
		dB[i,:] .= C .* ((i-n).*rt.^max(0,n-i-1).*t.^(i-1) .+ rt.^(n-i).*(i-1).*t.^max(0,i-2)) # max helps with derivatives of t^0
	end


	points = control * B
	tangents = control * dB
	normals = [0 -1; 1 0] * tangents
	d2 = sum(abs2, normals; dims=1)
	d = max.(1e-12, sqrt.(d2))
	normals .*= width./(2.0.*d)

	BezierTrajectory(control,width,points,normals)
end


function polygon(l::LineTrajectory)
	p1 = l.points[:,1]
	p2 = l.points[:,2]
	v = p2 .- p1
	v ./= norm(v)
	n = [-v[2], v[1]]

	d = n*(l.width/2)
	# [p1.+d p2.+d p2.-d p1.-d]
	# [p1.+d p2.+d p2.-d p1.-d p1.+d] # closed path
	[p1.-d p2.-d p2.+d p1.+d p1.-d] # closed path
end

# polygon(b::BezierTrajectory) =
# 	hcat(b.points .+ b.normals, reverse(b.points .- b.normals; dims=2))

function polygon(b::BezierTrajectory)
	# fwd = b.points .+ b.normals
	# rev = reverse(b.points .- b.normals; dims=2)
	# hcat(fwd, rev, fwd[:,1]) # closed path
	fwd = b.points .- b.normals
	rev = reverse(b.points .+ b.normals; dims=2)
	hcat(fwd, rev, fwd[:,1]) # closed path
end



# Outdated code below, needs minor adjustments to make it work


# trajectoryprojection(data::DataMatrix, t::Trajectory) =
# 	trajectoryprojection(obs_coordinates(data), t)


# """
# 	trajectoryprojection(points::Matrix{Float64}, l::LineTrajectory)

# Returns a float in [0,1] for each point, or nothing for points outside the line segment.
# """
# function trajectoryprojection(points::Matrix{Float64}, l::LineTrajectory)
# 	p1 = l.points[:,1]
# 	p2 = l.points[:,2]
# 	v = p2 .- p1

# 	P = points .- p1
# 	t = (P'v)./sum(abs2,v)
# 	dist2 = vec(sum(abs2, P .- v.*t'; dims=1)) # orthogonal distance to the line squared (this formulation should work in any dimension)
# 	Union{Nothing,Float64}[0≤s≤1 && d2<=(l.width/2)^2 ? s : nothing for (s,d2) in zip(t,dist2)]
# end


# _t_or_nothing(t1::Real,t2::Real) = (t1+t2)/2
# _t_or_nothing(t1::Real,::Nothing) = t1
# _t_or_nothing(::Nothing,t2::Real) = t2
# _t_or_nothing(::Nothing,::Nothing) = nothing
# _scale_t(t::Real, j, dt) = (t+j-1)*dt
# _scale_t(::Nothing, j, dt) = nothing

# function _merge_t(t::AbstractVector)
# 	# there should normally be one unique t value (all the other nothing)
# 	# but we allow two consecutive values to be something, in which case we take the mean
# 	ind = findall(!isnothing,t)
# 	isempty(ind) && return nothing
# 	length(ind)==1 && return t[ind[1]]
# 	length(ind)==2 && return (t[ind[1]]+t[ind[2]])/2
# 	error("Cell is inside multiple line segments (ind=$ind)")
# end


# function trajectoryprojection(points::Matrix{Float64}, b::BezierTrajectory)
# 	npoints = size(points,2)
# 	nsegments = size(b.points,2)-1

# 	dt = 1/nsegments

# 	t = Matrix{Union{Nothing,Float64}}(nothing, npoints, nsegments)

# 	# NB: This is an approximate way to get `t`. Use with many small segments.
# 	# TODO: Find correct way to compute `t`.
# 	for j in 1:nsegments
# 		c1 = b.points[:,j]
# 		c2 = b.points[:,j+1]
# 		n1 = b.normals[:,j]
# 		n2 = b.normals[:,j+1]

# 		p1 = c1-n1
# 		p2 = c1+n1
# 		p3 = c2+n2
# 		p4 = c2-n2

# 		# first triangle
# 		u1 = p4-p1
# 		v1 = p2-p1

# 		q1 = points .- p1
# 		ts1 = [u1 v1] \ q1

# 		ε = 1e-6

# 		inside1 = vec(all(>=(-ε), ts1; dims=1) .& (sum(ts1; dims=1).<=(1+ε)))
# 		t1 = ifelse.(inside1, @view(ts1[1,:]), nothing)

# 		# second triangle
# 		u2 = p2-p3
# 		v2 = p4-p3

# 		q2 = points .- p3
# 		ts2 = [u2 v2] \ q2

# 		inside2 = vec(all(>=(-ε), ts2; dims=1) .& (sum(ts2; dims=1).<=(1+ε)))
# 		t2 = ifelse.(inside2, 1.0.-@view(ts2[1,:]), nothing)

# 		t[:,j] .= _scale_t.(_t_or_nothing.(t1,t2), j, dt)
# 	end

# 	# handle points that are within multiple segments
# 	t = _merge_t.(eachrow(t))
# end

# """
# 	obs2time(data::DataMatrix, t::AbstractVector; ntimepoints=101)

# Create a new DataMatrix, where observations are replace by time points.
# Takes a time value `t` for each observation (in `[0,1]` or `nothing`), and merges observations using linear interpolation.
# """
# function obs2time(data::DataMatrix, t::AbstractVector; ntimepoints=101)
# 	P,N = size(data)
# 	obs_ind = findall(!isnothing,t)
# 	# N2 = count(!isnothing,t)
# 	N2 = length(obs_ind)

# 	timepoints = range(0,1,ntimepoints)

# 	# create a sparse matrix for computing linear interpolation
# 	I = zeros(Int, 2N2)
# 	J = zeros(Int, 2N2)
# 	V = zeros(2N2)
# 	for (k,i) in enumerate(obs_ind)
# 		j = searchsortedlast(timepoints, t[i])
# 		j1 = max(j, 1)
# 		j2 = min(j+1, ntimepoints)
# 		α = (t[i] - timepoints[j1]) / step(timepoints)
# 		I[2k-1] = i
# 		I[2k]   = i
# 		J[2k-1] = j1
# 		J[2k]   = j2
# 		V[2k-1] = 1-α
# 		V[2k]   = α
# 	end
# 	S = sparse(I,J,V,N,ntimepoints)

# 	X = data.matrix*S
# 	w = ones(1,N)*S
# 	X ./= max.(w,1e-12)

# 	obs = DataFrame(id=string.("t_",timepoints), time=timepoints, ncells=round.(vec(w);digits=2))

# 	update_matrix(data, X; var=:copy, obs, obs_id_cols=["id"])
# end



# # Experimental feature to get the "orthogonal" coordinate for each point on a trajectory
# # NB: assumes points with t===nothing have been filtered out
# function trajectory_perpendicular_coord(points::Matrix{Float64}, t::Vector{<:Real}, b::BezierTrajectory)
# 	@assert size(points,1)==2

# 	nparts = size(b.points,2)

# 	timepoints = range(0, 1; length=nparts)
# 	dt = (timepoints[2]-timepoints[1])

# 	i1 = max.(1,      1 .+ round.(Int,t.*(nparts-1),RoundDown))
# 	i2 = min.(nparts, 1 .+ round.(Int,t.*(nparts-1),RoundUp))

# 	# @show i1,i2

# 	# α = (t.-getindex.(Ref(timepoints), i1)) ./ dt
# 	α = (t.-timepoints[i1]) ./ dt

# 	m = (1 .- α').*b.points[:,i1] .+ α'.*b.points[:,i2]
# 	n = (1 .- α').*b.normals[:,i1] .+ α'.*b.normals[:,i2]


# 	d = points .- m # dims × npoints
# 	s = sum(d.*n; dims=1) # scalar product between d and normal for each point
# 	s ./= sum(abs2, n; dims=1)
# 	s = vec(s)

# 	# TESTING
# 	# @show dot(d[:,1],n[:,1]) / (norm(d[:,1])*norm(n[:,1]))


# 	# TODO: test code to ensure we can get the original `points` back from `t` and `s`. Doesn't hold for the moment due to approximate `t` estimate (see above).
# 	# reconstructed = m .+ s' .* n
# 	# @show reconstructed ≈ points

# 	s
# end

# function point_on_trajectory(t, s, b::BezierTrajectory)
# 	# allow ϵ outside?
# 	@assert 0 <= t <= 1
# 	@assert -1 <= s <= 1

# 	nparts = size(b.points,2)

# 	timepoints = range(0, 1; length=nparts)
# 	dt = (timepoints[2]-timepoints[1])

# 	i1 = max(1,      1 + round(Int,t*(nparts-1),RoundDown))
# 	i2 = min(nparts, 1 + round(Int,t*(nparts-1),RoundUp))

# 	# @show i1,i2

# 	α = (t-timepoints[i1]) / dt

# 	# midpoint and normal at `t`
# 	m = (1-α)*b.points[:,i1] + α*b.points[:,i2]
# 	n = (1-α)*b.normals[:,i1] + α*b.normals[:,i2]

# 	m + n*s
# end
