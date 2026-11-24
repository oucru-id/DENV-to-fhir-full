LOINC = "http://loinc.org"
SNOMED = "http://snomed.info/sct"
UCUM = "http://unitsofmeasure.org"
HGVS = "http://varnomen.hgvs.org"
SEQUENCE_ONTOLOGY = "http://www.sequenceontology.org"

CANONICAL_BASE = "http://terminology.spheres.id"
CS_DENV_SEROTYPE = f"{CANONICAL_BASE}/CodeSystem/denv-serotype"
CS_DENV_LINEAGE = f"{CANONICAL_BASE}/CodeSystem/denv-lineage"
CS_OBSERVATION_CODE = f"{CANONICAL_BASE}/CodeSystem/denvfhir-observation-code"
CS_PIPELINE_PROPERTY = f"{CANONICAL_BASE}/CodeSystem/pipeline-property"

SEROTYPE_DISPLAYS = {
    "DENV-1": "Dengue virus serotype 1",
    "DENV-2": "Dengue virus serotype 2",
    "DENV-3": "Dengue virus serotype 3",
    "DENV-4": "Dengue virus serotype 4",
    "unknown": "Unknown serotype",
}

OBSERVATION_CODE_DISPLAYS = {
    "dengue-classification": "Dengue virus genomic classification",
    "serotype": "Dengue virus serotype",
    "genotype": "Dengue virus genotype",
    "major-lineage": "Dengue virus major lineage",
    "minor-lineage": "Dengue virus minor lineage",
    "clade": "Nextclade lineage code",
    "classification-confidence": "Classification confidence",
    "coinfection-flag": "Possible mixed-serotype co-infection flag",
    "sylvatic-flag": "Sylvatic lineage reference match",
    "serotype-read-fraction": "Fraction of DENV reads supporting primary serotype",
    "viral-consensus-sequence": "Viral consensus genome sequence",
    "sequence-length": "Consensus sequence length",
    "sequence-completeness": "Consensus sequence completeness (non-N bases)",
    "ambiguous-bases": "Ambiguous bases (N content)",
    "gc-content": "GC content",
    "coverage-breadth": "Fraction of reference bases at or above minimum depth",
    "sequencing-qc": "Sequencing and consensus QC status",
}

VARIANT_TYPE_SO = {
    "substitution": {"code": "SO:0001583", "display": "missense_variant"},
    "deletion": {"code": "SO:0001822", "display": "inframe_deletion"},
    "insertion": {"code": "SO:0001821", "display": "inframe_insertion"},
}
