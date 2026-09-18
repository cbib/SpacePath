"""Minimal smoke test for smpath core workflow."""
import os
import tempfile

import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
import pytest


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
        "mz": [100.0, 200.0, 300.0],
        "formula": ["H2O", "CO2", "CH4"],
        "moleculeIds": ["HMDB0000001", "HMDB0000002", "HMDB0000003"],
        "moleculeNames": ["Water", "Carbon dioxide", "Methane"],
        "fdr": [0.1, 0.1, 0.1],
        "msm": [0.8, 0.8, 0.8],
        "rawOffSampleProb": [0.1, 0.1, 0.1],
        "totalIntensity": [1000.0, 1000.0, 1000.0]
    })
    annot_path = os.path.join(tmpdir, "test-HMDB.csv")
    annot_df.to_csv(annot_path, index=False)

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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])