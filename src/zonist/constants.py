
from typing import List, Dict

#

supported_mappings: List[str] = ["LipidMaps", "HMDB"]  # this list >
# > states the types of ref tables that match all not-kegg IDs with KEGG ones

supported_annotations: List[str] = ["LipidMaps", "HMDB", "KEGG", "custom_kegg"]

columns_ref_map_expected: List[str] = ["not_kegg_id", "KEGG", "database"]
