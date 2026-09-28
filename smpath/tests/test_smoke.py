"""Minimal smoke test for smpath core workflow."""
import os
import tempfile

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import pytest

from smpath.processing import measures_methods


def create_minimal_test_data(tmpdir):
    """Create minimal valid test data."""
    # Create minimal h5ad
    n_spots, n_mz = 10, 5
    np.random.seed(42)
    X = csr_matrix(np.random.rand(n_spots, n_mz) * 1000)
    coords = np.array([[i // 5, i % 5] for i in range(n_spots)])
    mz_values = [100.0, 200.0, 300.0, 400.0, 500.0]

    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"x{i//5}_y{i%5}" for i in range(n_spots)]),
        var=pd.DataFrame({"mz": mz_values}, index=[str(m) for m in mz_values]),
        obsm={"spatial": coords}
    )
    adata_path = os.path.join(tmpdir, "test.h5ad")
    adata.write(adata_path)

    # Create minimal annotation CSV
    annot_df = pd.DataFrame({
        "mz": [100.0, 200.0, 300.0, 400.0, 500.0 ],
        "formula": ["C18H36O", "C40H80O", "H2O", "CO2", "CH4", ],
        "moleculeIds": ["HMDB0000001", "HMDB0000002", "HMDB0000003",
                        "HMDB0000004", "HMDB0000005"],
        "moleculeNames": ["AA", "FA", "Water", "Carbon dioxide", "Methane" ],
        "fdr": [0.2, 0.2, 0.2, 0.2, 0.5],
        "msm": [0.8, 0.7, 0.6, 0.8, 0.7],
        "rawOffSampleProb": [0.1, 0.1, 0.3, 0.2, 0.8],
        "totalIntensity": [1000.0, 1000.0, 1000.0, 1000.0, 1000.0 ]
    })
    annot_path = os.path.join(tmpdir, "test-HMDB.csv")
    annot_df.to_csv(annot_path, index=False)

    # Create minimal pathways file
    path = os.path.join(tmpdir, "test.gmt")
    pathway_data = pd.DataFrame({
        0: ["pathway1", "pathway2"],
        1: ["Arachidonate synthesis", "FAO"],
        2: ["HMDB0000001", "HMDB0000003"],
        3: ["HMDB0000002", "HMDB0000004"]
    })
    pathway_data.to_csv(path, sep='\t', index=False,
                        header=False)

    return tmpdir


class TestImports:
    """Test that core modules import without error."""

    def test_import_data(self):
        import smpath.data as spd
        assert hasattr(spd, 'MetabolomeAnnotations')
        assert hasattr(spd, 'MetabolicPathways')

    def test_import_method(self):
        import smpath.method as spm
        assert hasattr(spm, 'NormalizeIntensities')
        assert hasattr(spm, 'ClassifyMZProtocol')

    def test_import_helpers(self):
        import smpath.helpers as helpers
        assert hasattr(helpers, 'setup_logger')

    def test_import_processing(self):
        import smpath.processing.measures_methods as mm
        import smpath.processing.pathways_methods as pm
        assert hasattr(mm, 'tic_normalization')


class TestMinimalWorkflow:
    """Test minimal data loading and normalization workflow."""

    def test_load_h5ad(self):
        """Test that h5ad files can be loaded."""
        with tempfile.TemporaryDirectory() as tmpdir:
            create_minimal_test_data(tmpdir)
            adata = ad.read_h5ad(os.path.join(tmpdir, "test.h5ad"))
            assert adata.shape == (10, 5)
            assert "mz" in adata.var.columns
            assert "spatial" in adata.obsm

    def test_load_annotations(self):
        """Test that annotations can be loaded."""
        import smpath.data as spd

        with tempfile.TemporaryDirectory() as tmpdir:
            create_minimal_test_data(tmpdir)

            annots = spd.MetabolomeAnnotations(
                name="test",
                data_dir=tmpdir,
                subfolder="",
                annotation_files={"HMDB": "test-HMDB.csv"},
                discard_plant_compounds=False,
                discard_pharmaceuticals=False,
                discard_other_non_animal=False,
                discard_odd_n_carb_lipids=False
            )
            annots.load()

            assert annots.no_translated.shape[0] > 0
            assert "mz" in annots.no_translated.columns
            assert "moleculeIds" in annots.no_translated.columns

    def test_normalize_intensities(self):
        """Test TIC normalization."""
        import smpath.method as spm

        with tempfile.TemporaryDirectory() as tmpdir:
            create_minimal_test_data(tmpdir)
            adata = ad.read_h5ad(os.path.join(tmpdir, "test.h5ad"))

            adata = spm.NormalizeIntensities.run(adata)

            assert "norm" in adata.uns
            assert adata.uns["norm"] == "TotalIonCount"
            assert "raw_X" in adata.layers

    def test_classify_mzs(self):
        import smpath.method as spm
        import smpath.data as spd

        sample_name = "test"
        with tempfile.TemporaryDirectory() as tmpdir:
            create_minimal_test_data(tmpdir)

            annots = spd.MetabolomeAnnotations(
                name=sample_name,
                data_dir=tmpdir,
                subfolder="",
                annotation_files={"HMDB": "test-HMDB.csv"}
            )
            annots.load()
            adata = ad.read_h5ad(os.path.join(tmpdir, "test.h5ad"))
            annots, adata = measures_methods.match_adata_mz_and_annots(
                annots, adata
            )
            adata = spm.NormalizeIntensities.run(adata)
            updated_config = spm.ClassifyMZProtocolConfig(
                msm_cutoff=0.2,
                offSample_cutoff=0.6,
                fdr_cutoff=0.2,
                moran_I_cutoff=-1, # disabled for testing
                moran_padj_cutoff=1, # disabled for testing
                n_neighs_for_graph=3
            )
            os.makedirs(os.path.join(tmpdir, sample_name)) # to tmp save plot
            adata = spm.ClassifyMZProtocol.run(
                adata, annots, sample_name, plot_folder=tmpdir,
                config=updated_config)

            assert np.all(
                np.array(adata.var['final_quality'].tolist()
                ) == np.array(["good", "good", "good", "good", "bad"]))
            assert 'features_n_after_mz_clean_1' in list(
                adata.uns['history'].keys())
            assert 'features_n_after_mz_clean_2' in list(
                adata.uns['history'].keys())
            assert 'moranI' in list(adata.uns.keys())

    def test_compute_assignments_weights(self):
        import smpath.method as spm
        import smpath.data as spd
        sample_name = "test"
        with tempfile.TemporaryDirectory() as tmpdir:
            create_minimal_test_data(tmpdir)
            annots = spd.MetabolomeAnnotations(
                name=sample_name,
                data_dir=tmpdir,
                subfolder="",
                annotation_files={"HMDB": "test-HMDB.csv"}
            )
            annots.load()
            adata = ad.read_h5ad(os.path.join(tmpdir, "test.h5ad"))
            annots, adata = measures_methods.match_adata_mz_and_annots(
                annots, adata
            )
            annots.usable_id = "moleculeIds"
            annots.verify()

            adata = spm.NormalizeIntensities.run(adata)

            pathways = spd.MetabolicPathways(
                data_dir=tmpdir,
                subfolder="",
                file="test.gmt"
            )
            pathways.molecule_id_type = "moleculeIds"  # not KEGG
            pathways.load()

            pathways, annots.df = spm.SyncPathwaysAndAnnots.run(
                pathways.model_copy(),
                annots.df,
                usable_id="moleculeIds",
                min_set_size=2)
            weights_df, coords = spm.ComputeAssignmentsWeights.run(
                adata, annots.df, usable_id="moleculeIds")

            assert weights_df.shape[0] == adata.obs.shape[0]
            assert weights_df.shape[1] == len(annots.df['moleculeIds'].unique())

    def test_compute_pathways_scores(self):
        import smpath.method as spm
        import smpath.data as spd

        n_spots = 50
        rng = np.random.default_rng(seed=42)
        molecules_id = ["HMDB0000001", "HMDB0000002", "HMDB0000003", "HMDB0000004"]
        shapes = [1.0, 1.3, 4.0, 3.7]
        scales = [0.5, 0.5, 2.0, 2.0]
        mat = rng.gamma(shapes, scales, size=(n_spots, 4)) * 500
        computed_weights_df = pd.DataFrame(
            data=mat,
            index=[f"x{i // 5}_y{i % 5}" for i in range(n_spots)],
            columns=molecules_id)

        with tempfile.TemporaryDirectory() as tmpdir:
            create_minimal_test_data(tmpdir)
            pathways = spd.MetabolicPathways(
                data_dir=tmpdir,
                subfolder="",
                file="test.gmt"
            )
            pathways.molecule_id_type = "moleculeIds"  # not KEGG
            pathways.load()

            config_run_scores_wiki = spm.ComputePathwaysScoresConfig(
                type_abbr="WIKI", min_entity=2)

            wiki_scores_df = spm.ComputePathwaysScores.run(
                pathways.to_dict(),
                computed_weights_df,
                config=config_run_scores_wiki)

            assert wiki_scores_df.shape[0] == n_spots
            assert wiki_scores_df.shape[1] == len(pathways.to_dict().keys())
            assert wiki_scores_df.min().min() < 0
            assert wiki_scores_df.max().max() > 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])