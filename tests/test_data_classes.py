"""Tests for core data classes."""
import os
import tempfile
import pytest
import pandas as pd
import numpy as np

import smpath.data as spd


class TestMetabolomeAnnotations:
    """Tests for MetabolomeAnnotations class."""

    def create_valid_csv(self, tmpdir, filename="test.csv"):
        """Create a valid annotation CSV."""
        df = pd.DataFrame({
            "mz": [100.0, 200.0],
            "formula": ["H2O", "CO2"],
            "moleculeIds": ["HMDB0000001", "HMDB0000002"],
            "moleculeNames": ["Water", "Carbon dioxide"],
            "fdr": [0.1, 0.1],
            "msm": [0.8, 0.8],
            "rawOffSampleProb": [0.1, 0.1],
            "totalIntensity": [1000.0, 1000.0]
        })
        path = os.path.join(tmpdir, filename)
        df.to_csv(path, index=False)
        return path

    def test_load_valid_csv(self):
        """Test loading valid annotation file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            self.create_valid_csv(tmpdir, "test-HMDB.csv")

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

            assert annots.no_translated.shape[0] == 2

    def test_load_missing_file_raises(self):
        """Test that missing file raises error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            annots = spd.MetabolomeAnnotations(
                name="test",
                data_dir=tmpdir,
                subfolder="",
                annotation_files={"HMDB": "nonexistent.csv"},
                discard_plant_compounds=False,
                discard_pharmaceuticals=False,
                discard_other_non_animal=False,
                discard_odd_n_carb_lipids=False
            )

            with pytest.raises(FileNotFoundError):
                annots.load()

    def test_load_wrong_delimiter_warns(self):
        """Test that wrong delimiter is detected."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create TSV but tell loader it's CSV
            df = pd.DataFrame({
                "mz": [100.0],
                "formula": ["H2O"],
                "moleculeIds": ["HMDB0000001"],
                "moleculeNames": ["Water"]
            })
            path = os.path.join(tmpdir, "test-HMDB.tsv")
            df.to_csv(path, sep="\t", index=False)

            annots = spd.MetabolomeAnnotations(
                name="test",
                data_dir=tmpdir,
                subfolder="",
                annotation_files={"HMDB": "test-HMDB.tsv"},
                delimiter=","  # Wrong delimiter
            )

            # Should raise or warn about insufficient columns
            with pytest.raises(Exception):
                annots.load()

class TestMetabolicPathways:
    """Tests for MetabolomePathways class."""
    def create_valid_pathways_gmt(self, tmpdir, filename="test.gmt"):
        path = os.path.join(tmpdir, filename)
        pathway_data = pd.DataFrame({
            0: ["pathway1", "pathway2"],
            1: ["Arachidonate synthesis", "FAO"],
            2: ["HMDB0000001", "HMDB0000003"],
            3: ["HMDB0000002", "HMDB0000004"]
        })
        pathway_data.to_csv(path, sep='\t', index=False,
                            header=False)
        return path

    def test_load_valid_pathways(self):
        import smpath.data as spd

        with tempfile.TemporaryDirectory() as tmpdir:
            self.create_valid_pathways_gmt(tmpdir, filename="test.gmt")
            pathways = spd.MetabolicPathways(
                data_dir=tmpdir,
                subfolder="",
                file="test.gmt"
            )
            pathways.molecule_id_type = "moleculeIds"  # not KEGG
            pathways.load()

            assert isinstance(pathways, spd.MetabolicPathways)
            assert isinstance(pathways.to_dict(), dict)
            assert pathways.sets_ids_df.columns[1] == 'pathway_name'
            assert np.all(np.array(pathways.piled.columns.tolist()
                                   ) == np.array(
                ['pathway_id', 'moleculeIds']))

    def test_load_pathways_missing_file_raises(self):
        """Test that missing file raises error."""
        with tempfile.TemporaryDirectory() as tmpdir:
            pathways = spd.MetabolicPathways(
                data_dir=tmpdir,
                subfolder="",
                file="nonexistent.mmt"
            )

            with pytest.raises(FileNotFoundError):
                pathways.load()


class TestMetaboliteIdsMappings:
    """Tests for MetaboliteIdsMappings class."""
    def create_valid_identifiers_mappings(self, tmpdir):
        filenames = ["HMDB-to-KEGG.tsv", "LipidMaps-to-KEGG.tsv",
                     "CoreMetabolome-to-KEGG.tsv"]
        for i in range(len(filenames)):
            path = os.path.join(tmpdir, filenames[i])
            if i == 0:
                data = pd.DataFrame({
                    "HMDB": ["HMDB0000008", "HMDB0000010", "HMDB0000011", "HMDB0000012"],
                    "KEGG": ["C05984", "C05299", "C01089", "C00526"],
                })
            elif i == 1:
                data = pd.DataFrame({
                    "LipidMaps": ["LMFA08040047", "LMFA01140018", "LMFA03060005",
                             "LMFA01050243"],
                    "KEGG": ["C13829", "C08366", "C04805", "C01089"],
                })
            elif i == 2:
                data = pd.DataFrame({
                    "CoreMetabolome": ["HMDB0304827", "HMDB0304940"],
                    "KEGG": ["C18724", "C06333"],
                })

            data.to_csv(path, sep='\t', index=False, header=True)

        return path

    def test_load_identifiers_mappings(self):
        import smpath.data as spd

        with tempfile.TemporaryDirectory() as tmpdir:
            self.create_valid_identifiers_mappings(tmpdir)
            mappings_to_kegg = spd.MetaboliteIdsMappings(
                data_dir=tmpdir,
                subfolder="",
                mapping_files={"HMDB": "HMDB-to-KEGG.tsv",
                               "LipidMaps": "LipidMaps-to-KEGG.tsv"}
            )
            mappings_to_kegg.load()
            assert isinstance(mappings_to_kegg, spd.MetaboliteIdsMappings)
            assert isinstance(mappings_to_kegg.table, pd.DataFrame)
            assert mappings_to_kegg