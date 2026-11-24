# Clinical metadata tables

Clinical metadata is optional. Without `--patient_metadata`, DENVtoFHIR emits only PHI-free genomic FHIR bundles (`fhir/genomic/<id>.genomic.fhir.json`) with placeholder patient references `Patient/<fhir_id>-patient`.

When `--patient_metadata` is supplied (CSV or `.xlsx`), DENVtoFHIR builds `Patient`, `Specimen`, `Encounter`, and `ServiceRequest` resources per matched sample, plus shared `Organization`, `Practitioner`, and `PractitionerRole` resources when `--organization_metadata` / `--practitioner_metadata` are provided, and links everything into self-contained `transaction` bundles in `fhir/bundle/<id>.bundle.fhir.json`.

Templates with the expected headers live in `assets/metadata_templates/`:
- `patient_clinical_metadata.csv`
- `organization_metadata.csv`
- `practitioner_metadata.csv`

Keep all real clinical metadata files, `fhir/clinical/`, and `fhir/bundle/` out of version control.
