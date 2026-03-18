using DataFrames, Statistics, StatsBase, LinearAlgebra, MultivariateStats
using PyPlot, Distributed, Random, CSV, Revise, Distributions, Dates, MultipleTesting

if length(ARGS) < 2
    error("Usage: julia src/julia/run_phase_inference.jl <data_path> <output_path>")
end

data_path = abspath(ARGS[1])
output_path = abspath(ARGS[2])
path_to_cyclops = joinpath(@__DIR__, "cyclops_core.jl")
verbose = get(ENV, "CPI_VERBOSE", "0") == "1"

# The first column of expression.csv is treated as the ordered gene list used
# during model fitting and downstream plotting. This matches the final dataset
# layout used in this project.
expression_data = CSV.read(joinpath(data_path, "expression.csv"), DataFrame)
seed_genes = String.(expression_data[!, 1])

function find_column(df::DataFrame, candidates::Vector{Symbol})
    existing = Set(Symbol.(names(df)))
    for candidate in candidates
        if candidate in existing
            return candidate
        end
    end
    return nothing
end

function extract_phase_hours(pred_df::DataFrame)
    sample_col = find_column(pred_df, [:ID, :Sample_ID])
    isnothing(sample_col) && error("Fit output must contain an ID or Sample_ID column.")

    phase_col = find_column(pred_df, [:Phase_MA, :Phases_MA, :Predicted_Phase_Hours, :Phase])
    isnothing(phase_col) && error("Fit output must contain a supported phase column.")

    raw_phase = [ismissing(value) ? NaN : Float64(value) for value in pred_df[!, phase_col]]
    valid_mask = .!isnan.(raw_phase)

    phase_hours = if phase_col == :Predicted_Phase_Hours
        mod.(raw_phase[valid_mask], 24.0)
    else
        mod.(raw_phase[valid_mask], 2π) .* (24.0 / (2π))
    end

    sample_ids = String.(pred_df[valid_mask, sample_col])
    return sample_ids, phase_hours, phase_col
end


sample_ids_with_collection_times = []
sample_collection_times = []

if ((length(sample_ids_with_collection_times)+length(sample_collection_times))>0) && 
   (length(sample_ids_with_collection_times) != length(sample_collection_times))
    error("Number of sample ids must match number of collection times.")
end

training_parameters = Dict(
    :regex_cont => r".*_C",
    :regex_disc => r".*_D",
    :blunt_percent => 0.975,
    :seed_min_CV => -Inf,
    :seed_max_CV => Inf,
    :seed_mth_Gene => 10000,
    :norm_gene_level => true,
    :norm_disc => false,
    :norm_disc_cov => 1,
    :eigen_reg => true,
    :eigen_reg_disc_cov => 1,
    :eigen_reg_exclude => false,
    :eigen_reg_r_squared_cutoff => 0.6,
    :eigen_reg_remove_correct => false,
    :eigen_first_var => false,
    :eigen_first_var_cutoff => 0.85,
    :eigen_total_var => 0.85,
    :eigen_contr_var => 0.05,
    :eigen_var_override => true,
    :eigen_max => 5,
    :out_covariates => true,
    :out_use_disc_cov => true,
    :out_all_disc_cov => true,
    :out_disc_cov => 1,
    :out_use_cont_cov => false,
    :out_all_cont_cov => true,
    :out_use_norm_cont_cov => false,
    :out_all_norm_cont_cov => true,
    :out_cont_cov => 1,
    :out_norm_cont_cov => 1,
    :init_scale_change => true,
    :init_scale_1 => false,
    :train_n_models => 80,
    :train_μA => 0.001,
    :train_β => (0.9, 0.999),
    :train_min_steps => 1500,
    :train_max_steps => 2050,
    :train_μA_scale_lim => 1000,
    :train_circular => false,
    :train_collection_times => true,
    :train_collection_time_balance => 1.0,
    :cosine_shift_iterations => 192,
    :cosine_covariate_offset => true,
    :align_p_cutoff => 0.05,
    :align_base => "radians",
    :align_disc => false,
    :align_disc_cov => 1,
    :align_other_covariates => false,
    :align_batch_only => false,
    :X_Val_k => 10,
    :X_Val_omit_size => 0.1,
    :plot_use_o_cov => true,
    :plot_correct_batches => true,
    :plot_disc => false,
    :plot_disc_cov => 1,
    :plot_separate => false,
    :plot_color => ["b", "orange", "g", "r", "m", "y", "k"],
    :plot_only_color => true,
    :plot_p_cutoff => 0.05
)

# Each worker needs the same environment and core implementation because the
# fitting step distributes model training across available CPU cores.
Distributed.addprocs(length(Sys.cpu_info()))
@everywhere begin
    using DataFrames, Statistics, StatsBase, LinearAlgebra, MultivariateStats
    using PyPlot, Random, CSV, Revise, Distributions, Dates, MultipleTesting
    include($path_to_cyclops)
end

mkpath(output_path)

println("Running circadian phase inference...")
eigendata, modeloutputs, correlations, bestmodel, parameters = 
    CYCLOPS.Fit(expression_data, seed_genes, training_parameters)

verbose && println("Aligning and saving results...")
CYCLOPS.Align(expression_data, modeloutputs, correlations, bestmodel, parameters, output_path)

# The plotting section is best-effort only. It should not fail the entire run
# if a dataset lacks metadata or if a plotting assumption is not satisfied.
verbose && println("Generating visualizations...")
clock_genes = ["ARNTL", "CLOCK", "PER1", "PER2", "PER3", "CRY1", "CRY2", "NR1D1", "NR1D2", "DBP", "TEF", "HLF"]

try
    fit_candidates = sort(filter(x -> startswith(x, "Fit_Output"), readdir(output_path)))
    isempty(fit_candidates) && error("No Fit_Output*.csv file was produced.")

    fit_file = joinpath(output_path, last(fit_candidates))
    pred_df = CSV.read(fit_file, DataFrame)
    ordered_sample_ids, ordered_phase_hours, used_phase_col = extract_phase_hours(pred_df)
    sample_to_phase = Dict(zip(ordered_sample_ids, ordered_phase_hours))
    verbose && println("Using $(used_phase_col) from fit output for visualization.")
    
    meta_file = joinpath(data_path, "metadata.csv")
    meta_df = nothing
    cell_types = [nothing]
    if isfile(meta_file)
        meta_df = CSV.read(meta_file, DataFrame)
        meta_names = Set(Symbol.(names(meta_df)))
        if :CellType_D in meta_names && :Sample in meta_names
            cell_types = collect(skipmissing(unique(meta_df[!, :CellType_D])))
        end
    end
    
    gene_col = findfirst(x -> occursin("gene", lowercase(string(x))) || x == :Gene_Symbol, names(expression_data))
    isnothing(gene_col) && (gene_col = 1)
    gene_ids = String.(expression_data[!, gene_col])
    
    for cell_type in cell_types
        samples = isnothing(cell_type) ? ordered_sample_ids : String.(meta_df[meta_df.CellType_D .== cell_type, :Sample])
        suffix = isnothing(cell_type) ? "" : "_$(cell_type)"
        
        available_genes = filter(g -> g in gene_ids, clock_genes)
        length(available_genes) == 0 && continue
        
        n_cols = min(3, length(available_genes))
        n_rows = ceil(Int, length(available_genes) / n_cols)
        fig, axes = subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))
        axes = length(available_genes) == 1 ? [axes] : vec(axes)
        
        for (idx, gene) in enumerate(available_genes)
            ax = axes[idx]
            gene_idx = findfirst(==(gene), gene_ids)
            isnothing(gene_idx) && continue
            
            phases, expressions = Float64[], Float64[]
            for sample in samples
                if haskey(sample_to_phase, sample) && sample in names(expression_data)
                    push!(phases, sample_to_phase[sample])
                    push!(expressions, expression_data[gene_idx, sample])
                end
            end
            
            length(phases) == 0 && continue
            ax.scatter(phases, expressions, alpha=0.6, s=30, edgecolors="black", linewidth=0.5)
            
            try
                x_fit = range(0, 24, length=1000)
                baseline, amp = mean(expressions), (maximum(expressions) - minimum(expressions)) / 2
                phase_guess = phases[argmax(expressions)]
                y_fit = baseline .+ amp .* cos.(2π .* (x_fit .- phase_guess) ./ 24)
                ax.plot(x_fit, y_fit, "-", linewidth=2, alpha=0.8)
                ax.set_title("$(gene)\nPhase: $(round(phase_guess, digits=1))h, Amp: $(round(amp, digits=1))", fontsize=10)
            catch
                ax.set_title(gene, fontsize=10)
            end
            
            ax.set_xlabel("Phase (hours)", fontsize=9)
            ax.set_ylabel("Expression", fontsize=9)
            ax.set_xlim(0, 24)
            ax.set_xticks([0, 6, 12, 18, 24])
            ax.grid(true, alpha=0.3)
        end
        
        for idx in (length(available_genes)+1):length(axes)
            axes[idx].axis("off")
        end
        
        suptitle("Clock Genes$(suffix)", fontsize=14, fontweight="bold")
        tight_layout()
        savefig(joinpath(output_path, "clock_genes$(suffix).png"), dpi=300, bbox_inches="tight")
        close()
    end
    
    figure(figsize=(10, 6))
    phase_df = DataFrame(Sample = ordered_sample_ids, Phase_Hours = ordered_phase_hours)
    if !isnothing(meta_df) && :CellType_D in Symbol.(names(meta_df)) && :Sample in Symbol.(names(meta_df))
        merged = leftjoin(phase_df, meta_df, on=:Sample => :Sample)
        for ct in sort(unique(merged.CellType_D))
            ct_values = skipmissing(merged[merged.CellType_D .== ct, :Phase_Hours])
            isempty(ct_values) && continue
            hist(ct_values, bins=24, alpha=0.6, label=string(ct), edgecolor="black")
        end
        legend()
    else
        hist(ordered_phase_hours, bins=24, edgecolor="black")
    end
    xlabel("Predicted Phase (hours)", fontsize=12)
    ylabel("Count", fontsize=12)
    title("Phase Distribution", fontsize=14, fontweight="bold")
    xlim(0, 24)
    grid(true, alpha=0.3)
    savefig(joinpath(output_path, "phase_distribution.png"), dpi=300, bbox_inches="tight")
    close()
    
    verbose && println("Visualization complete!")
catch e
    verbose && println("Visualization skipped: $e")
end

println("Phase inference completed. Results written to $(output_path)")