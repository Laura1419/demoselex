#!/usr/bin/env python3
"""Retrieve literature evidence only; never reads biological genotype/outcome data."""
from __future__ import annotations
import json
from pathlib import Path
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import provenance

# Explicitly disable ALL git probes before constructing the repository recorder.
provenance.get_git_info = lambda: {"disabled": True, "reason": "operator forbids git commands"}
OUT = Path("data/derived/local_ancestral_state/literature")
QUERIES = {
    "ELAI": 'TITLE:"Detecting structure of haplotypes and local ancestry"',
    "RFMix": 'TITLE:RFMix AND FIRST_PDATE:[2013-01-01 TO 2013-12-31]',
    "FLARE": 'TITLE:"Fast, accurate local ancestry inference with FLARE"',
    "LAMP_LD": 'TITLE:"Fast and accurate inference of local ancestry in Latino populations"',
    "Gnomix": '"High-resolution ancestry deconvolution"',
    "ASAFE": 'TITLE:"ASAFE: ancestry-specific allele frequency estimation"',
    "Tractor": 'TITLE:Tractor AND AUTH_FIRST:Atkinson',
    "MIXSCORE": 'TITLE:"Admixture mapping for genetic diseases"',
    "ALDER": 'TITLE:"Inferring admixture histories of human populations using linkage disequilibrium"',
    "CLUES": 'TITLE:"Estimating the full allele frequency trajectory"',
    "Relate": 'TITLE:"A method for genome-wide genealogy estimation for thousands of samples"',
    "tsinfer": 'TITLE:"Inferring whole-genome histories in large population datasets"',
    "ancient_imputation": 'TITLE:"Imputation of ancient human genomes"',
    "ancient_LAI": '(TITLE:"local ancestry" AND ABSTRACT:"ancient DNA")',
    "Lanc_map": 'TITLE:"Inference of locus-specific ancestry in closely related populations"',
    "Ancestry_HMM": 'TITLE:"A hidden Markov model approach for simultaneously estimating local ancestry and admixture time"',
    "MOSAIC": 'TITLE:"Fine-scale inference of ancestry segments without prior knowledge of admixing groups"',
    "time_series": 'TITLE:"Estimating a full selection model from allele frequency time series"',
}


def fetch(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "trajectory-method-study/1.0"})
    with urllib.request.urlopen(req, timeout=25) as response:
        return response.read()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    recorder = provenance.RunRecorder("local-ancestral-state-literature")
    recorder.record_input(__file__)
    recorder.record_parameter("queries", QUERIES)
    index = {"retrieved_utc": datetime.now(timezone.utc).isoformat(), "sources": {}}
    for name, query in QUERIES.items():
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
            {"query": query, "format": "json", "resultType": "core", "pageSize": 5})
        try:
            body = fetch(url)
            path = OUT / f"{name}.json"
            path.write_bytes(body)
            recorder.record_output(path)
            result = json.loads(body)
            rows = result.get("resultList", {}).get("result", [])
            entry = {"query": query, "url": url, "hits": result.get("hitCount"), "papers": []}
            for row in rows:
                info = {key: row.get(key) for key in ("title", "authorString", "pubYear", "doi", "pmcid", "id")}
                pmcid = row.get("pmcid")
                if pmcid:
                    full_url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
                    try:
                        xml = fetch(full_url)
                        xml_path = OUT / f"{name}_{pmcid}.xml"
                        xml_path.write_bytes(xml)
                        recorder.record_output(xml_path)
                        root = ET.fromstring(xml)
                        text_path = OUT / f"{name}_{pmcid}.txt"
                        text_path.write_text("\n\n".join("".join(node.itertext()) for node in root.iter()
                                                        if node.tag in ("article-title", "title", "p")))
                        recorder.record_output(text_path)
                        info["full_text_url"] = full_url
                    except Exception as error:
                        info["full_text_error"] = str(error)
                entry["papers"].append(info)
            index["sources"][name] = entry
            print(name, entry["hits"], [p["title"] for p in entry["papers"]], flush=True)
        except Exception as error:
            index["sources"][name] = {"query": query, "url": url, "error": str(error)}
            print(name, "ERROR", str(error), flush=True)
    index_path = OUT / "index.json"
    index_path.write_text(json.dumps(index, indent=2) + "\n")
    recorder.record_output(index_path)
    recorder.save()


if __name__ == "__main__":
    main()
