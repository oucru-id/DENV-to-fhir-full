from denvfhir_terminology import (CS_DENV_LINEAGE, CS_DENV_SEROTYPE, CS_OBSERVATION_CODE, LOINC,
                                  OBSERVATION_CODE_DISPLAYS, SEROTYPE_DISPLAYS)

GENOMICS_SD = "http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition"
FHIR_OBSERVATION_SD = "http://hl7.org/fhir/StructureDefinition/Observation"

KEMKES_SD = "https://fhir.kemkes.go.id/r4/StructureDefinition"
KEMKES_SP = "http://terminology.kemkes.go.id/sp"
KEMKES_TAG = {"system": KEMKES_SP, "code": "genomics", "display": "Genomics"}
KEMKES_SYS_IDS = "http://sys-ids.kemkes.go.id"
KEMKES_NIK = "https://fhir.kemkes.go.id/id/nik"

REPORT_CODE_V1 = {"system": LOINC, "code": "81247-9", "display": "Master HL7 genetic variant reporting panel"}
REPORT_CODE_IG3 = {"system": LOINC, "code": "51969-4", "display": "Genetic analysis report"}
REPORT_TEXT = "Dengue Virus Genomic Analysis Report"


class Profile:
    name = "generic"
    include_nik = False
    kemkes_extensions = False
    specimen_defaults = False

    def __init__(self, identifier_system_base=None):
        self.identifier_system_base = (identifier_system_base or "").rstrip("/") or None

    def meta(self, kind):
        profiles = self._profiles(kind)
        meta = {}
        if profiles:
            meta["profile"] = profiles
        tags = self._tags(kind)
        if tags:
            meta["tag"] = tags
        return meta or None

    def _profiles(self, kind):
        return {
            "variant": [f"{GENOMICS_SD}/variant"],
            "region_studied": [f"{GENOMICS_SD}/region-studied"],
            "classification": [FHIR_OBSERVATION_SD],
            "consensus": [FHIR_OBSERVATION_SD],
            "report": [f"{GENOMICS_SD}/genomic-report"],
        }.get(kind, [])

    def _tags(self, kind):
        return []

    def bundle_meta(self, last_updated):
        return {"lastUpdated": last_updated}

    def observation_code(self, code, extra_codings=None, text=None):
        codings = list(extra_codings or [])
        codings.append({
            "system": CS_OBSERVATION_CODE,
            "code": code,
            "display": OBSERVATION_CODE_DISPLAYS[code],
        })
        concept = {"coding": codings}
        concept["text"] = text or OBSERVATION_CODE_DISPLAYS[code]
        return concept

    def component_code(self, code):
        return {"coding": [{
            "system": CS_OBSERVATION_CODE,
            "code": code,
            "display": OBSERVATION_CODE_DISPLAYS[code],
        }]}

    def serotype_concept(self, serotype, text=None):
        code = serotype if serotype in SEROTYPE_DISPLAYS else "unknown"
        return {
            "coding": [{
                "system": CS_DENV_SEROTYPE,
                "code": code,
                "display": SEROTYPE_DISPLAYS[code],
            }],
            "text": text or serotype,
        }

    def clade_concept(self, clade):
        if not clade or clade in ("unknown", "Indeterminate", "Unassigned"):
            return {"text": clade or "Indeterminate"}
        return {
            "coding": [{
                "system": CS_DENV_LINEAGE,
                "code": clade,
                "display": clade,
            }],
            "text": clade,
        }

    def gene_concept(self, gene_name):
        return {"text": gene_name}

    def report_code(self):
        return {"coding": [dict(REPORT_CODE_IG3), dict(REPORT_CODE_V1)], "text": REPORT_TEXT}

    def identifier(self, kind, value, org_id=None):
        if not self.identifier_system_base or value in (None, ""):
            return None
        return {"system": f"{self.identifier_system_base}/{kind}", "value": str(value)}


class SatusehatProfile(Profile):
    name = "satusehat"
    include_nik = True
    kemkes_extensions = True
    specimen_defaults = True

    _KEMKES = {
        "patient": "Patient",
        "specimen": "Specimen",
        "organization": "Organization",
        "practitioner": "Practitioner",
        "practitioner_role": "PractitionerRole",
        "service_request": "ServiceRequest",
        "encounter": "Encounter",
    }

    def _profiles(self, kind):
        if kind in self._KEMKES:
            return [f"{KEMKES_SD}/{self._KEMKES[kind]}"]
        if kind == "report":
            return [f"{GENOMICS_SD}/genomics-report"]
        return super()._profiles(kind)

    def _tags(self, kind):
        if kind in ("variant", "region_studied", "classification", "consensus", "report", "device", "provenance"):
            return [dict(KEMKES_TAG)]
        return []

    def bundle_meta(self, last_updated):
        return {"lastUpdated": last_updated, "profile": [f"{KEMKES_SD}/Bundle"]}

    def report_code(self):
        return {"coding": [dict(REPORT_CODE_V1)], "text": REPORT_TEXT}

    def identifier(self, kind, value, org_id=None):
        if value in (None, ""):
            return None
        if kind == "organization":
            return {"use": "official", "system": f"{KEMKES_SYS_IDS}/organization", "value": str(value)}
        if not org_id:
            return None
        system_kind = {"service_request": "servicerequest"}.get(kind, kind)
        return {"system": f"{KEMKES_SYS_IDS}/{system_kind}/{org_id}", "value": str(value)}


PROFILES = {"generic": Profile, "satusehat": SatusehatProfile}


def get_profile(name, identifier_system_base=None):
    try:
        return PROFILES[name](identifier_system_base=identifier_system_base)
    except KeyError:
        raise SystemExit(f"unknown --fhir_profile {name!r}; choose one of {sorted(PROFILES)}")
