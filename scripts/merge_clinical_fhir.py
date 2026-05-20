#!/usr/bin/env python3

import json
import argparse
import uuid
import sys
import os  
from datetime import datetime, timezone
from clinical_metadata_parser import (
    load_clinical_metadata, find_matching_sample, get_clinical_value,
    load_organization_metadata, load_practitioner_metadata
)
import base64
import re

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def create_patient_resource(sample_id, clinical_data=None, org_data=None):
    if not clinical_data:
        raise ValueError(f"Clinical data is required for sample {sample_id}")

    org_data = org_data or {}
    org_id = org_data.get('org_id')
    family_name = get_clinical_value(clinical_data, 'family_name')
    given_name = get_clinical_value(clinical_data, 'given_name')
    gender = get_clinical_value(clinical_data, 'gender', 'unknown').lower()
    birth_date = get_clinical_value(clinical_data, 'birth_date')
    nik = get_clinical_value(clinical_data, 'nik')
    province_code = get_clinical_value(clinical_data, 'province_code')
    city_code = get_clinical_value(clinical_data, 'city_code')
    district_code = get_clinical_value(clinical_data, 'district_code')
    village_code = get_clinical_value(clinical_data, 'village_code')
    citizenship_status = get_clinical_value(clinical_data, 'citizenship_status')
    lat = get_clinical_value(clinical_data, 'latitude', None)
    lon = get_clinical_value(clinical_data, 'longitude', None)

    if gender in ['laki-laki', 'pria', 'male', 'm']:
        gender = "male"
    elif gender in ['perempuan', 'wanita', 'female', 'f']:
        gender = "female"
    else:
        gender = "unknown"

    geo_extensions = []
    if lat and lon:
        try:
            geo_extensions = [{
                "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [
                    {"url": "latitude",  "valueDecimal": float(lat)},
                    {"url": "longitude", "valueDecimal": float(lon)}
                ]
            }]
        except ValueError:
            pass

    return {
        "resourceType": "Patient",
        "id": f"{sample_id}-patient",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
        },
        "active": True,
        "name": [
            {
                "use": "official",
                "family": family_name,
                "given": [given_name]
            }
        ],
        "gender": gender,
        "birthDate": birth_date,
        "identifier": [
            {
                "use": "official",
                "system": "https://fhir.kemkes.go.id/id/nik",
                "value": nik
            },
            {
                "use": "usual",
                "type": {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                            "code": "MR",
                            "display": "Medical record number"
                        }
                    ]
                },
                "system": f"http://sys-ids.kemkes.go.id/mr/{org_id}",
                "value": sample_id
            }
        ],
        "extension": [
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                "extension": [
                    {"url": "province", "valueCode": province_code},
                    {"url": "city",     "valueCode": city_code},
                    {"url": "district", "valueCode": district_code},
                    {"url": "village",  "valueCode": village_code}
                ]
            },
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/citizenshipStatus",
                "valueCode": citizenship_status
            }
        ],
        "address": [
            {
                "use": "home",
                "type": "physical",
                "text": get_clinical_value(clinical_data, 'address'),
                "city": get_clinical_value(clinical_data, 'city'),
                "state": get_clinical_value(clinical_data, 'state'),
                "country": "ID",
                "extension": [
                    {
                        "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                        "extension": [
                            {"url": "province", "valueCode": province_code},
                            {"url": "city",     "valueCode": city_code},
                            {"url": "district", "valueCode": district_code},
                            {"url": "village",  "valueCode": village_code}
                        ]
                    },
                    *geo_extensions
                ]
            }
        ]
    }

def create_organization_resource(org_data=None):
    org_data = org_data or {}
    org_id        = org_data.get('org_id', 'unknown-org')
    name          = org_data.get('name', 'Unknown Organization')
    alias         = org_data.get('alias', '')
    type_code     = org_data.get('type_code', '')
    type_display  = org_data.get('type_display', '')
    type_text     = org_data.get('type_text', '')
    phone         = org_data.get('phone', '')
    email         = org_data.get('email', '')
    address_line  = org_data.get('address_line', '')
    city          = org_data.get('city', '')
    state         = org_data.get('state', '')
    country       = org_data.get('country', 'ID')
    province_code = org_data.get('province_code', '')
    city_code     = org_data.get('city_code', '')
    district_code = org_data.get('district_code', '')
    lat           = org_data.get('latitude', None)
    lon           = org_data.get('longitude', None)

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})
    if email:
        telecom.append({"system": "email", "value": email, "use": "work"})

    addr_extensions = []
    if province_code or city_code or district_code:
        code_ext = {"url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode", "extension": []}
        if province_code:
            code_ext["extension"].append({"url": "province", "valueCode": province_code})
        if city_code:
            code_ext["extension"].append({"url": "city",     "valueCode": city_code})
        if district_code:
            code_ext["extension"].append({"url": "district", "valueCode": district_code})
        addr_extensions.append(code_ext)
    if lat and lon:
        try:
            addr_extensions.append({
                "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                "extension": [
                    {"url": "latitude",  "valueDecimal": float(lat)},
                    {"url": "longitude", "valueDecimal": float(lon)}
                ]
            })
        except ValueError:
            pass

    resource = {
        "resourceType": "Organization",
        "id": org_id,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Organization"]
        },
        "identifier": [
            {
                "use": "official",
                "system": "http://sys-ids.kemkes.go.id/organization",
                "value": org_id
            }
        ],
        "active": True,
        "type": [
            {
                "coding": [
                    {
                        "system": "http://terminology.kemkes.go.id/CodeSystem/organization-type",
                        "code": type_code,
                        "display": type_display
                    }
                ],
                "text": type_text
            }
        ],
        "name": name,
        "telecom": telecom,
        "address": [
            {
                "use": "work",
                "type": "physical",
                "line": [address_line] if address_line else [],
                "city": city,
                "state": state,
                "country": country,
                "extension": addr_extensions
            }
        ]
    }
    if alias:
        resource["alias"] = [alias]
    return resource

def create_practitioner_resource(practitioner_data=None):
    practitioner_data = practitioner_data or {}
    pid        = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    nik        = practitioner_data.get('nik', '')
    name       = practitioner_data.get('name', 'Unknown Practitioner')
    phone      = practitioner_data.get('phone', '')
    gender     = practitioner_data.get('gender', 'unknown')
    birth_date = practitioner_data.get('birth_date', '')
    str_kki    = practitioner_data.get('str_kki_number', '')
    qual_start = practitioner_data.get('qualification_period_start', '')

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})

    qualification = []
    if str_kki:
        qual = {
            "code": {
                "coding": [{
                    "system": "https://terminology.kemkes.go.id/v1-0302",
                    "code": "STR-KKI",
                    "display": "Surat Tanda Registrasi Dokter"
                }],
                "text": "Surat Tanda Registrasi Dokter"
            }
        }
        qual["identifier"] = [{"system": "https://fhir.kemkes.go.id/id/str-kki-number", "value": str_kki}]
        if qual_start:
            qual["period"] = {"start": qual_start}
        qualification.append(qual)

    resource = {
        "resourceType": "Practitioner",
        "id": pid,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Practitioner"]
        },
        "active": True,
        "name": [{"use": "official", "text": name}],
        "telecom": telecom,
        "gender": gender
    }
    if nik:
        resource["identifier"] = [{
            "use": "official",
            "system": "https://fhir.kemkes.go.id/id/nik",
            "value": nik
        }]
    if birth_date:
        resource["birthDate"] = birth_date
    if qualification:
        resource["qualification"] = qualification
    return resource

def create_practitioner_role_resource(practitioner_data=None, org_data=None):
    practitioner_data = practitioner_data or {}
    org_data = org_data or {}
    pid          = practitioner_data.get('practitioner_id', 'unknown-practitioner')
    pname        = practitioner_data.get('name', 'Unknown Practitioner')
    phone        = practitioner_data.get('phone', '')
    role_id      = practitioner_data.get('role_id', 'unknown-role')
    role_code    = practitioner_data.get('role_code', '')
    role_display = practitioner_data.get('role_display', '')
    role_text    = practitioner_data.get('role_text', '')
    org_id       = org_data.get('org_id', 'unknown-org')
    org_name     = org_data.get('name', 'Unknown Organization')

    telecom = []
    if phone:
        telecom.append({"system": "phone", "value": phone, "use": "work"})

    return {
        "resourceType": "PractitionerRole",
        "id": role_id,
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/PractitionerRole"]
        },
        "active": True,
        "practitioner": {
            "reference": f"Practitioner/{pid}",
            "display": pname
        },
        "organization": {
            "reference": f"Organization/{org_id}",
            "display": org_name
        },
        "code": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": role_code,
                        "display": role_display
                    }
                ],
                "text": role_text
            }
        ],
        "telecom": telecom
    }

def create_specimen_resource(sample_id, clinical_data=None, practitioner_data=None, org_data=None,
                             spec_type_code='119297000', spec_type_display='Blood specimen',
                             method_code='BLD', method_display='Blood', method_text='Blood collection'):
    practitioner_data = practitioner_data or {}
    org_data = org_data or {}
    pid   = practitioner_data.get('practitioner_id', '')
    pname = practitioner_data.get('name', '')
    org_id = org_data.get('org_id', '')

    if clinical_data:
        given_name  = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
        collection_date = get_clinical_value(clinical_data, 'collection_date', '')
    else:
        patient_display = f"Patient {sample_id}"
        collection_date = ''

    collection = {
        "method": {
            "coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0488",
                        "code": method_code, "display": method_display}],
            "text": method_text
        }
    }
    if pid:
        collection["collector"] = {"reference": f"Practitioner/{pid}", "display": pname}
    if collection_date:
        collection["collectedDateTime"] = collection_date

    identifier = []
    if org_id:
        identifier.append({"system": f"http://sys-ids.kemkes.go.id/specimen/{org_id}", "value": sample_id})

    resource = {
        "resourceType": "Specimen",
        "id": f"{sample_id}-specimen",
        "meta": {"profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Specimen"]},
        "type": {
            "coding": [{"system": "http://snomed.info/sct",
                        "code": spec_type_code, "display": spec_type_display}]
        },
        "subject": {"reference": f"Patient/{sample_id}-patient", "display": patient_display},
        "collection": collection
    }
    if identifier:
        resource["identifier"] = identifier
    return resource


def create_service_request_resource(sample_id, clinical_data=None, practitioner_data=None, org_data=None):
    practitioner_data = practitioner_data or {}
    org_data = org_data or {}
    pid      = practitioner_data.get('practitioner_id', '')
    pname    = practitioner_data.get('name', '')
    role_id  = practitioner_data.get('role_id', '')
    org_id   = org_data.get('org_id', 'unknown-org')

    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"

    resource = {
        "resourceType": "ServiceRequest",
        "id": f"{sample_id}-service-request",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/ServiceRequest"]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/servicerequest/{org_id}",
                "value": f"SR-{sample_id}"
            }
        ],
        "status": "active",
        "intent": "original-order",
        "priority": "routine",
        "category": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "108252007",
                        "display": "Laboratory procedure"
                    }
                ]
            }
        ],
        "code": {
            "coding": [
                {
                    "system": "http://loinc.org",
                    "code": "69548-6",
                    "display": "Genetic variant assessment"
                }
            ],
            "text": "Dengue Genetic Variant Assessment"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "Dengue Testing Encounter"
        },
        "occurrenceDateTime": datetime.now(timezone.utc).isoformat(),
        "reasonCode": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "38362002",
                        "display": "Dengue fever"
                    }
                ],
                "text": "Suspected Dengue virus infection"
            }
        ]
    }
    if pid:
        resource["requester"] = {"reference": f"Practitioner/{pid}", "display": pname}
    if role_id:
        resource["performer"] = [{"reference": f"PractitionerRole/{role_id}", "display": "Laboratory Technician"}]
    return resource

def extract_dengue_classification(observations):
    serotype = "Unknown"
    genotype = "Unknown"
    major_lineage = "Unknown"
    minor_lineage = "Unknown"
    confidence = "low"
    
    for obs in observations:
        if 'consensus' in obs.get('id', ''):
            continue

        code_coding = obs.get('code', {}).get('coding', [])
        is_classification = False
        for c in code_coding:
            if c.get('code') == '31343-7':
                is_classification = True
                break
        
        if not is_classification:
            code_text = obs.get('code', {}).get('text', '').lower()
            if 'dengue virus classification' in code_text:
                is_classification = True

        if is_classification:
            for comp in obs.get('component', []):
                comp_code = comp.get('code', {})
                comp_text = comp_code.get('text', '').lower()
                
                if 'serotype' in comp_text:
                    serotype = comp.get('valueCodeableConcept', {}).get('text', 'Unknown')
                
                if 'genotype' in comp_text:
                    genotype = comp.get('valueCodeableConcept', {}).get('text', 'Unknown')

                if 'major lineage' in comp_text:
                    major_lineage = comp.get('valueCodeableConcept', {}).get('text', 'Unknown')

                if 'minor lineage' in comp_text:
                    minor_lineage = comp.get('valueCodeableConcept', {}).get('text', 'Unknown')

                if 'confidence' in comp_text:
                    confidence = comp.get('valueCodeableConcept', {}).get('text', 'low')
            
            break
    
    return serotype, genotype, major_lineage, minor_lineage, confidence

def create_diagnostic_report(sample_id, observations, clinical_data=None, org_data=None, practitioner_data=None):
    org_data = org_data or {}
    practitioner_data = practitioner_data or {}
    org_id   = org_data.get('org_id', 'unknown-org')
    org_name = org_data.get('name', '')
    pid      = practitioner_data.get('practitioner_id', '')
    pname    = practitioner_data.get('name', '')

    serotype, genotype, major_lineage, minor_lineage, confidence = extract_dengue_classification(observations)
    
    conclusion = f"Dengue virus serotype {serotype} detected. Genotype: {genotype}. Lineage: {major_lineage} / {minor_lineage}. Confidence: {confidence.upper()}"
    
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"
    
    report_id = f"{sample_id}-genomic-report"
    current_time = datetime.now(timezone.utc).isoformat()
    
    # Simplified HTML Content
    html_content = f"""<div xmlns="http://www.w3.org/1999/xhtml">
<h1>Dengue Virus Genomic Analysis Report</h1>
<p><strong>Patient:</strong> {patient_display}</p>
<p><strong>Sample ID:</strong> {sample_id}</p>
<p><strong>Report Date:</strong> {current_time}</p>
<p><strong>Serotype:</strong> {serotype}</p>
<p><strong>Genotype:</strong> {genotype}</p>
<p><strong>Major Lineage:</strong> {major_lineage}</p>
<p><strong>Minor Lineage:</strong> {minor_lineage}</p>
<p><strong>Confidence:</strong> {confidence.upper()}</p>
<p><strong>Conclusion:</strong> {conclusion}</p>
</div>"""
    
    html_base64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')
    
    return {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"],
            "tag": [
                {
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "genomics",
                    "display": "Genomics"
                }
            ]
        },
        "identifier": [
            {
                "system": f"http://sys-ids.kemkes.go.id/diagnostic-report/{org_id}",
                "value": f"DENGUE-GEN-{sample_id}-{datetime.now().strftime('%Y%m%d')}"
            }
        ],
        "status": "final",
        "category": [
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                    "code": "GE",
                    "display": "Genetics"
                }]
            }
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "81247-9",
                "display": "Master HL7 genetic variant reporting panel"
            }],
            "text": "Dengue Virus Genomic Analysis Report"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "Dengue Testing Encounter"
        },
        "effectiveDateTime": current_time,
        "issued": current_time,
        "performer": [
            *([{"reference": f"Organization/{org_id}", "display": org_name}] if org_id else []),
            *([{"reference": f"Practitioner/{pid}",   "display": pname}]    if pid    else [])
        ],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations if obs.get('id')],
        "specimen": [{
            "reference": f"Specimen/{sample_id}-specimen",
            "display": f"Blood specimen from {patient_display}"
        }],
        "conclusion": conclusion,
        "conclusionCode": [
            {"text": f"Serotype {serotype}"},
            {"text": f"Genotype {genotype}"},
            {"text": f"Major Lineage {major_lineage}"},
            {"text": f"Minor Lineage {minor_lineage}"},
            {"text": f"Confidence: {confidence}"}
        ],
        "presentedForm": [
            {
                "contentType": "text/html",
                "language": "en-US",
                "title": "Dengue Virus Genomic Analysis Report",
                "data": html_base64
            }
        ]
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to input FHIR bundle')
    parser.add_argument('--output', required=True, help='Path to output merged FHIR bundle')
    parser.add_argument('--patient_metadata', help='Path to patient/clinical metadata CSV/Excel file')
    parser.add_argument('--organization_metadata', help='Path to organization metadata CSV file')
    parser.add_argument('--practitioner_metadata', help='Path to practitioner metadata CSV file')
    args = parser.parse_args()

    clinical_data = {}
    if args.patient_metadata and os.path.exists(args.patient_metadata):
        clinical_data = load_clinical_metadata(args.patient_metadata)
    else:
        debug_print(f"Patient metadata file not found or not provided")

    org_data = {}
    if args.organization_metadata and os.path.exists(args.organization_metadata):
        org_data = load_organization_metadata(args.organization_metadata)
    else:
        debug_print(f"Organization metadata file not found or not provided")

    practitioner_data = {}
    if args.practitioner_metadata and os.path.exists(args.practitioner_metadata):
        practitioner_data = load_practitioner_metadata(args.practitioner_metadata)
    else:
        debug_print(f"Practitioner metadata file not found or not provided")

    try:
        with open(args.input, 'r') as f:
            fhir_bundle = json.load(f)

        sample_ids = set()
        all_observations = []
        
        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            if resource.get('resourceType') == 'Observation':
                all_observations.append(resource)
                subject_ref = resource.get('subject', {}).get('reference', '')
                
                if subject_ref.startswith('Patient/') or 'patient' in subject_ref:
                    sample_id = subject_ref.replace('Patient/', '').replace('urn:uuid:patient-', '').replace('-patient', '')
                    sample_ids.add(sample_id)

        if not sample_ids:
            filename = os.path.basename(args.input)
            filename_sample_id = filename.replace('.fhir.json', '').replace('_ont', '').replace('_illumina', '')
            sample_ids.add(filename_sample_id)

        matched_samples = {}
        for sample_id in sample_ids:
            sample_clinical_data = find_matching_sample(sample_id, clinical_data)
            if sample_clinical_data:
                matched_samples[sample_id] = sample_clinical_data

        merged_bundle = {
            "resourceType": "Bundle",
            "id": str(uuid.uuid4()),
            "meta": {
                "lastUpdated": datetime.now(timezone.utc).isoformat(),
                "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]
            },
            "type": "transaction",
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "entry": []
        }

        org_resource = create_organization_resource(org_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": org_resource,
            "request": {
                "method": "PUT",
                "url": f"Organization/{org_resource['id']}"
            }
        })

        practitioner_resource = create_practitioner_resource(practitioner_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": practitioner_resource,
            "request": {
                "method": "PUT",
                "url": f"Practitioner/{practitioner_resource['id']}"
            }
        })

        role_resource = create_practitioner_role_resource(practitioner_data, org_data)
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": role_resource,
            "request": {
                "method": "PUT",
                "url": f"PractitionerRole/{role_resource['id']}"
            }
        })

        for sample_id, sample_clinical_data in matched_samples.items():
            patient_resource = create_patient_resource(sample_id, sample_clinical_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": patient_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Patient/{patient_resource['id']}"
                }
            })

            specimen_resource = create_specimen_resource(sample_id, sample_clinical_data, practitioner_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": specimen_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Specimen/{specimen_resource['id']}"
                }
            })

            service_request_resource = create_service_request_resource(sample_id, sample_clinical_data, practitioner_data, org_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": service_request_resource,
                "request": {
                    "method": "PUT",
                    "url": f"ServiceRequest/{service_request_resource['id']}"
                }
            })

        observations_by_sample = {}
        for obs in all_observations:
            subject_ref = obs.get('subject', {}).get('reference', '')
            if 'patient' in subject_ref.lower():
                sample_id = subject_ref.replace('Patient/', '').replace('urn:uuid:patient-', '').replace('-patient', '')
                if sample_id not in observations_by_sample:
                    observations_by_sample[sample_id] = []
                observations_by_sample[sample_id].append(obs)

        for sample_id, sample_observations in observations_by_sample.items():
            sample_clinical_data = matched_samples.get(sample_id)
            
            diagnostic_report = create_diagnostic_report(
                sample_id,
                sample_observations,
                sample_clinical_data,
                org_data,
                practitioner_data
            )
            
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": diagnostic_report,
                "request": {
                    "method": "PUT",
                    "url": f"DiagnosticReport/{diagnostic_report['id']}"
                }
            })
            
        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            resource_type = resource.get('resourceType')
            resource_id = resource.get('id')
            
            entry_with_request = {
                "fullUrl": entry.get('fullUrl', f"urn:uuid:{str(uuid.uuid4())}"),
                "resource": resource,
                "request": {
                    "method": "PUT",
                    "url": f"{resource_type}/{resource_id}" if resource_id else f"{resource_type}"
                }
            }
            merged_bundle['entry'].append(entry_with_request)

        with open(args.output, 'w') as f:
            json.dump(merged_bundle, f, indent=2)

    except Exception as e:
        debug_print(f"Error occurred: {str(e)}")
        import traceback
        debug_print(f"Traceback: {traceback.format_exc()}")

if __name__ == "__main__":
    main()