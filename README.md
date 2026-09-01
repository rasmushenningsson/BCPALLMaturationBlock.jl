# BCPALLMaturationBlock.jl
The SingleCellProjections.jl analysis performed in the paper ["Single-cell genomics details the maturation block in BCP-ALL and identifies therapeutic vulnerabilities in DUX4-r cases" by Thorsson et al.](https://doi.org/10.1182/blood.2023021705)

## Installation
Add the package:
```julia-repl
julia> using Pkg

julia> Pkg.add(;url="https://github.com/rasmushenningsson/BCPALLMaturationBlock.jl.git")
```

## Data
Samples and annotation files can be downloaded from the [release page](https://github.com/rasmushenningsson/BCPALLMaturationBlock.jl/releases#release-samples_and_annotations).

Then run
```julia-repl
julia> set_bcpall_path("path/to/samples_folder")
julia> set_annotations_path("path/to/annotations_folder")
```
once to set the data locations.

## Usage

```julia
using BCPALLMaturationBlock
```

Generate plots using the `nbm_reference_plot`, `nbm_reference_trajectory_histogram`, `projection_scatter_plot` and `projection_trajectory_histogram` functions.
