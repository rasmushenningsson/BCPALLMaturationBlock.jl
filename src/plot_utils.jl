
function scatter_2d(job)
    matrix = fetch!(SCP.get_matrix(job))
    fig = Figure(; size=(768, 768))
    ax = Axis(fig[1, 1])
    scatter!(ax, matrix; color = :black, markersize = 4)
    fig
end

function scatter_categorical_2d(job, annot_name; bg=nothing, colors=nothing)
    matrix = fetch!(SCP.get_matrix(job))
    annot = fetch!(SCP.value_column_data(SCP.annotation(SCP.get_obs(job), annot_name)))

    fig = Figure(; size=(768, 768))
    ax = Axis(fig[1, 1])

    if bg !== nothing
        bg_matrix = fetch!(SCP.get_matrix(bg))
        scatter!(ax, bg_matrix; color=colorant"#BFCCE6", markersize=2)
    end

    if colors !== nothing
        unique_annotations = unique(annot)
        unique_annotations_set = Set(unique_annotations)
        colors = filter(x->x[1] in unique_annotations_set, colors)

        categories = first.(colors)
        @assert isempty(setdiff(unique_annotations, categories)) # ensure all categories have colors specified

        plots = [scatter!(ax, matrix[:,isequal.(annot, cat)]; markersize=6, color) for (cat,color) in colors]
    else
        categories = unique(annot)
        plots = [scatter!(ax, matrix[:,isequal.(annot, cat)]; markersize=6) for cat in categories]
    end

    axislegend(ax, plots .=> Ref((;markersize=16)), categories) # use a larger markersize in the legend
    fig
end
