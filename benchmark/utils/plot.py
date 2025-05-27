import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_benchmark(df_results: pd.DataFrame):
    """Plot the benchmark results.

    Parameters
    ----------
    df_results : pd.DataFrame
        DataFrame containing the benchmark results, with models as columns and metrics as rows.

    Returns
    -------
    fig, axes : tuple
        Figure and axes objects for the plot.
    """
    metric_index = df_results.index.get_level_values(0).unique()
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(
        figsize=(10 * len(metric_index), 10),
        nrows=1,
        ncols=len(metric_index),
    )
    if len(metric_index) == 1:
        axes = [axes]

    for index, ax in zip(metric_index, axes, strict=False):
        sub_df = df_results.loc[index]
        bp = sns.boxplot(
            sub_df.T,
            ax=ax,
            width=0.8,
            dodge=False,
            palette=sns.color_palette("pastel", n_colors=len(sub_df.T.columns)),
            flierprops={"marker": "o", "markersize": 2, "linestyle": "none", "markerfacecolor": "gray"},
        )
        add_median(bp, sub_df.T, color="black")
        ax.tick_params(axis="x", labelrotation=90)
        ax.set_title(index)
        ax.set_xlabel("Models")
        ax.set_ylabel(index)
        ax.grid(False)

    return fig, axes


def add_median(bp, data, color="black"):
    """Add median lines to the boxplot."""
    for i, column in enumerate(data.columns):
        median = data[column].median()
        bp.text(i, median, f"{median:.3f}", color=color, ha="center", va="bottom", fontsize=10)


def plot_benchmark_ncells(df_results: pd.DataFrame):
    df_results = df_results.reorder_levels(
        [1, 2, 0]
    ).sort_index()  # we want to have the index in order: metric, model, n_cells

    metric_index = df_results.index.get_level_values(0).unique()
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(
        figsize=(5 * len(metric_index), 10),
        nrows=1,
        ncols=len(metric_index),
    )
    if len(metric_index) == 1:
        axes = [axes]

    for index, ax in zip(metric_index, axes, strict=False):
        sub_df = df_results.loc[index]
        sub_df = sub_df.T.melt()
        sub_df.columns = ["model", "n_cells", index]
        sns.lineplot(sub_df, x="n_cells", y=index, hue="model", ax=ax)
        ax.set_title(index)

        ax.grid(False)

    return fig, axes


def plot_celltype_correlation_heatmap(ground_truth: pd.DataFrame, prediction: pd.DataFrame, method: str = "pearson"):
    """
    Plot a heatmap of correlations between predicted and ground truth cell type proportions.

    Parameters:
    - ground_truth: pd.DataFrame of shape (n_samples, n_cell_types)
    - prediction: pd.DataFrame of shape (n_samples, n_cell_types)
    - method: correlation method, either 'pearson' or 'spearman'

    Returns:
    - correlation_matrix: pd.DataFrame of shape (n_cell_types_gt, n_cell_types_pred)
    """
    # Validate shapes
    if ground_truth.shape != prediction.shape:
        raise ValueError("ground_truth and prediction must have the same shape (same samples and columns)")

    # Reorder prediction columns to match ground truth
    if set(ground_truth.columns) != set(prediction.columns):
        raise ValueError("ground_truth and prediction must have the same set of columns (cell types)")

    prediction = prediction[ground_truth.columns]

    # Compute correlation matrix
    correlation_matrix = pd.DataFrame(index=ground_truth.columns, columns=prediction.columns)

    for gt_col in ground_truth.columns:
        for pred_col in prediction.columns:
            if method == "pearson":
                corr = ground_truth[gt_col].corr(prediction[pred_col])
            elif method == "spearman":
                corr = ground_truth[gt_col].corr(prediction[pred_col], method="spearman")
            else:
                raise ValueError("method must be 'pearson' or 'spearman'")
            correlation_matrix.loc[gt_col, pred_col] = corr

    # Convert to float
    correlation_matrix = correlation_matrix.astype(float)

    # Plot heatmap
    fig = plt.figure(figsize=(8, 6))
    sns.heatmap(
        correlation_matrix,
        annot=True,
        fmt=".2f",
        cmap="coolwarm",
        center=0,
        xticklabels=prediction.columns,
        yticklabels=ground_truth.columns,
    )
    plt.title(f"Correlation map: Ground Truth vs Prediction ({method.title()})")
    plt.xlabel("Predicted cell types")
    plt.ylabel("Ground truth cell types")
    plt.tight_layout()

    return fig
