from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from tqdm import tqdm
from utils.generate_data import load_pseudobulk_dataset
from utils.metric import correlation, mse
from utils.plot import plot_benchmark, plot_benchmark_ncells

from pydeconv.model import NNLS, OLS, RLR, MixupVI, NuSVR, Scaden, Tape
from pydeconv.model.base_model import NeuralNetworkModel, SignatureBasedModel
from pydeconv.signature_matrix.registry import (
    sig_matrix_cti_granularity_1,
    sig_matrix_cti_granularity_2,
)


FCT_SIGNATURE_MATRIX = {
    "1st_level_granularity": sig_matrix_cti_granularity_1,
    "2nd_level_granularity": sig_matrix_cti_granularity_2,
}


def benchmark_models(
    folder_data: str, model_classes: list, metrics: list, granularity: str, n_cells: int, debug: bool = True
):
    """
    Run the benchmark for the given models on the pseudobulk dataset.

    Parameters
    ----------
    folder_data : str
        Path to the folder containing the pseudobulk dataset.
    model_classes : list
        List of model classes to benchmark.
    metrics : list
        List of metric functions to evaluate the models.
    granularity : str
        Granularity level for the signature matrix.
    n_cells : int
        Number of cells to use in the pseudobulk dataset.
    Returns
    -------
    pd.DataFrame
        DataFrame containing the results of the benchmark, with models as columns and metrics as rows.
    """

    # Load the pseudobulk dataset

    pseudobulk_test, labels = load_pseudobulk_dataset(
        n_sample=400, n_cells=n_cells, granularity=granularity, folder_data=Path(folder_data)
    )

    all_results_celltypes = {}
    all_results_patients = {}
    signature_matrix = FCT_SIGNATURE_MATRIX[granularity]()
    for model_class in tqdm(model_classes, desc="Benchmarking models", position=0, leave=False):
        if issubclass(model_class, NeuralNetworkModel):
            model = model_class(weights_version=f"cti_{granularity}")
        elif issubclass(model_class, SignatureBasedModel):
            model = model_class(signature_matrix=signature_matrix)
        results = model.transform(pseudobulk_test, ratio=False)
        from utils.plot import plot_celltype_correlation_heatmap

        if debug:
            fig = plot_celltype_correlation_heatmap(labels, results)
            fig.savefig(f"benchmark_results/debug_heatmap_{model._name}.png", dpi=300)
            plt.close()
        all_results_celltypes[model._name] = {metric.__name__: metric(results, labels) for metric in metrics}
        all_results_patients[model._name] = {metric.__name__: metric(results.T, labels.T) for metric in metrics}

    all_results_celltypes = pd.concat(
        {model: pd.DataFrame(results).T for model, results in all_results_celltypes.items()}
    )
    all_results_patients = pd.concat(
        {model: pd.DataFrame(results).T for model, results in all_results_patients.items()}
    )
    all_results_celltypes = all_results_celltypes.swaplevel(0, 1).sort_index()
    all_results_patients = all_results_patients.swaplevel(0, 1).sort_index()
    return all_results_celltypes, all_results_patients


def run(
    folder_data: str,
    granularities: list,
    models: list,
    metrics: list,
    list_n_cells: list,
    output_folder: Path = Path(""),
    benchmark_granularities: bool = True,
    benchmark_n_cells: bool = True,
    debug: bool = False,
):
    output_folder.mkdir(parents=True, exist_ok=True)

    if benchmark_granularities:
        print("Benchmarking for different granularities...")

        for granularity in tqdm(granularities, desc="Benchmarking granularities", position=1, leave=False):
            _output_folder = output_folder / granularity
            _output_folder.mkdir(parents=True, exist_ok=True)
            _, df_results_patients = benchmark_models(
                folder_data=folder_data,
                model_classes=models,
                metrics=metrics,
                granularity=granularity,
                n_cells=100,
                debug=debug,
            )

            # Save results and generate plots
            df_results_patients.to_csv(_output_folder / f"benchmark_granularities_results_{granularity}_patients.csv")

            fig, _ = plot_benchmark(df_results_patients)
            fig.suptitle(f"Benchmarking results for {granularity}")
            fig.savefig(_output_folder / f"benchmark_granularities_results_{granularity}_patients.png", dpi=300)
            plt.close()

    if benchmark_n_cells:
        print("Benchmarking for different numbers of cells...")

        for granularity in tqdm(granularities, desc="Benchmarking granularities", position=2, leave=False):
            _output_folder = output_folder / granularity
            _output_folder.mkdir(parents=True, exist_ok=True)

            results_patients = {}

            for n_cells in tqdm(list_n_cells, desc="Benchmarking n_cells", position=1, leave=False):
                _, df_results_patients = benchmark_models(
                    folder_data=folder_data,
                    model_classes=models,
                    metrics=metrics,
                    granularity=granularity,
                    n_cells=n_cells,
                    debug=debug,
                )
                results_patients[n_cells] = df_results_patients

            df_results_patients = pd.concat(results_patients)

            # Save results and generate plots
            df_results_patients.to_csv(_output_folder / f"benchmark_n_cells_results_{granularity}_celltypes.csv")

            fig, _ = plot_benchmark_ncells(df_results_patients)
            fig.suptitle(f"Benchmarking results for {granularity} - Patients")
            fig.savefig(_output_folder / f"benchmark_n_cells_results_{granularity}_patients.png", dpi=300)
            plt.close()


if __name__ == "__main__":
    run(
        granularities=["1st_level_granularity", "2nd_level_granularity"],
        models=[NNLS, OLS, RLR, Tape, Scaden, MixupVI, NuSVR],  # DWLS, WNNLS],
        metrics=[correlation, mse],
        list_n_cells=[10, 25, 50, 75, 100, 125, 150, 200, 250, 300, 400, 500, 750, 1000],
        folder_data=Path("benchmark/data"),
        output_folder=Path("benchmark/benchmark_results"),
        benchmark_granularities=True,
        benchmark_n_cells=False,
        debug=False,
    )
