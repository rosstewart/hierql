"""
HierQL experiments.

Prerequisites:
  pip install psycopg2-binary rdflib neo4j sparqlwrapper
  brew install souffle-lang/souffle/souffle   # or equivalent

  Run with python3.10 (or whichever interpreter has the packages above).

Usage:
  python3.10 experiments.py --dataset sample
  python3.10 experiments.py --dataset real

For the real dataset, Neo4j must be reachable at NEO4J_URL (see below).
Start it with:
  docker run --rm --name neo4j -p 7474:7474 -p 7687:7687 \\
    --env NEO4J_AUTH=neo4j/password123 neo4j:latest

Set DB_URL to your PostgreSQL 16 connection string if needed.
"""

import argparse
import csv
import re
import gzip
import json
import math
import os
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from pathlib import Path

from hierql_transpiler import transpile

DB_URL     = "postgresql://postgres:password@localhost:5434/bioql"
NEO4J_URL  = "bolt://localhost:7687"
NEO4J_AUTH = ("neo4j", "password123")

REAL_GNOMAD_PATH = "gnomad.exomes.v4.1.1.sites.chrY.vcf.bgz"
REAL_GENCC_PATH  = "gencc-submissions.tsv"
REAL_MONDO_PATH  = "mondo.json"

# ── Sample schema + data ─────────────────────────────────────────────────────

SETUP_SQL_SAMPLE = """
DROP TABLE IF EXISTS gnomAD, Variant, Gene_Disease, Gene, DiseaseEdge, Disease CASCADE;

CREATE TABLE Disease (
    id        SERIAL PRIMARY KEY,
    name      TEXT NOT NULL,
    parent_id INT  REFERENCES Disease(id)
);
CREATE INDEX ON Disease(parent_id);
CREATE INDEX ON Disease(name);

CREATE TABLE DiseaseEdge (
    parent_id INT  REFERENCES Disease(id),
    child_id  INT  REFERENCES Disease(id),
    rel       TEXT NOT NULL,
    PRIMARY KEY (parent_id, child_id, rel)
);
CREATE INDEX ON DiseaseEdge(parent_id, rel);
CREATE INDEX ON DiseaseEdge(child_id, rel);

CREATE TABLE Gene (
    id     SERIAL PRIMARY KEY,
    symbol TEXT NOT NULL
);

CREATE TABLE Gene_Disease (
    gene_id    INT REFERENCES Gene(id),
    disease_id INT REFERENCES Disease(id),
    PRIMARY KEY (gene_id, disease_id)
);

CREATE TABLE Variant (
    id      SERIAL PRIMARY KEY,
    gene_id INT  REFERENCES Gene(id),
    rsid    TEXT
);
CREATE INDEX ON Variant(gene_id);

CREATE TABLE gnomAD (
    id          SERIAL PRIMARY KEY,
    variant_id  INT   REFERENCES Variant(id),
    allele_freq FLOAT NOT NULL
);
CREATE INDEX ON gnomAD(variant_id);
CREATE INDEX ON gnomAD(allele_freq);

INSERT INTO Disease(id, name, parent_id) VALUES
    (1, 'Metabolic Disease',  NULL),
    (2, 'Diabetes Mellitus',  1),
    (3, 'Obesity',            1),
    (4, 'Type 1 Diabetes',    2),
    (5, 'Type 2 Diabetes',    2),
    (6, 'MODY',               2);

INSERT INTO DiseaseEdge(parent_id, child_id, rel) VALUES
    (1, 2, 'is_a'),
    (1, 3, 'is_a'),
    (2, 4, 'is_a'),
    (2, 5, 'is_a'),
    (2, 6, 'is_a'),
    (1, 4, 'part_of'),
    (1, 5, 'part_of');

INSERT INTO Gene(id, symbol) VALUES
    (1, 'INS'), (2, 'HNF1A'), (3, 'GCK'), (4, 'PPARG');

INSERT INTO Gene_Disease(gene_id, disease_id) VALUES
    (1, 4), (2, 6), (3, 6), (4, 5);

INSERT INTO Variant(id, gene_id, rsid) VALUES
    (1, 1, 'rs7903146'), (2, 2, 'rs1169288'),
    (3, 3, 'rs4607517'), (4, 4, 'rs1801282');

INSERT INTO gnomAD(id, variant_id, allele_freq) VALUES
    (1, 1, 0.005), (2, 2, 0.002),
    (3, 3, 0.150), (4, 4, 0.008);
"""

SETUP_SQL_REAL = """
DROP TABLE IF EXISTS gnomAD, Variant, Gene_Disease, Gene, DiseaseEdge, Disease CASCADE;

CREATE TABLE Disease (
    id   TEXT PRIMARY KEY,
    name TEXT NOT NULL
);
CREATE INDEX ON Disease(name);

CREATE TABLE DiseaseEdge (
    parent_id TEXT REFERENCES Disease(id),
    child_id  TEXT REFERENCES Disease(id),
    rel       TEXT NOT NULL,
    PRIMARY KEY (parent_id, child_id, rel)
);
CREATE INDEX ON DiseaseEdge(parent_id, rel);
CREATE INDEX ON DiseaseEdge(child_id, rel);

CREATE TABLE Gene (
    id     SERIAL PRIMARY KEY,
    symbol TEXT NOT NULL UNIQUE
);

CREATE TABLE Gene_Disease (
    gene_id         INT  REFERENCES Gene(id),
    disease_id      TEXT REFERENCES Disease(id),
    classification  TEXT,
    PRIMARY KEY (gene_id, disease_id)
);
CREATE INDEX ON Gene_Disease(disease_id);

CREATE TABLE Variant (
    id      BIGSERIAL PRIMARY KEY,
    gene_id INT REFERENCES Gene(id),
    rsid    TEXT,
    chrom   TEXT NOT NULL,
    pos     INT  NOT NULL,
    ref     TEXT NOT NULL,
    alt     TEXT NOT NULL
);
CREATE INDEX ON Variant(gene_id);
CREATE INDEX ON Variant(chrom, pos);

CREATE TABLE gnomAD (
    id          BIGSERIAL PRIMARY KEY,
    variant_id  BIGINT REFERENCES Variant(id),
    allele_freq FLOAT NOT NULL,
    filt        TEXT
);
CREATE INDEX ON gnomAD(variant_id);
CREATE INDEX ON gnomAD(allele_freq);
"""

# ── Sample benchmark queries ─────────────────────────────────────────────────

SAMPLE_QUERIES = [
    ("Q1", "Variants with AF < 0.001 in gene INS",
     """
     SELECT v.rsid, a.allele_freq
     FROM   Variant v
     JOIN   gnomAD  a ON a.variant_id = v.id
     JOIN   Gene    g ON g.id         = v.gene_id
     WHERE  g.symbol = 'INS'
       AND  a.allele_freq < 0.001
     """),

    ("Q2", "Genes with > 0 variants (relaxed for sample data)",
     """
     SELECT g.symbol, COUNT(*) AS n_variants
     FROM   Gene    g
     JOIN   Variant v ON v.gene_id = g.id
     GROUP  BY g.symbol
     HAVING COUNT(*) > 0
     ORDER  BY n_variants DESC
     """),

    ("Q3", "All Disease descendants of 'Diabetes Mellitus' (UNDER)",
     """
     WITH RECURSIVE sub(id) AS (
         SELECT id FROM Disease WHERE name = 'Diabetes Mellitus'
         UNION
         SELECT d.id FROM Disease d JOIN sub s ON d.parent_id = s.id
     )
     SELECT d.name FROM Disease d JOIN sub s ON d.id = s.id
     """),

    ("Q4", "Ancestors of 'Type 2 Diabetes' (ANCESTORS)",
     """
     WITH RECURSIVE anc(id) AS (
         SELECT parent_id FROM Disease WHERE name = 'Type 2 Diabetes'
         UNION
         SELECT d.parent_id FROM Disease d JOIN anc a ON d.id = a.id
         WHERE d.parent_id IS NOT NULL
     )
     SELECT d.name FROM Disease d JOIN anc a ON d.id = a.id
     """),

    ("Q5", "Running example: rare variants in Diabetes subtypes",
     """
     WITH RECURSIVE sub(id) AS (
         SELECT id FROM Disease WHERE name = 'Diabetes Mellitus'
         UNION
         SELECT d.id FROM Disease d JOIN sub s ON d.parent_id = s.id
     )
     SELECT g.symbol, v.rsid, a.allele_freq
     FROM   Gene          g
     JOIN   Gene_Disease  gd ON gd.gene_id    = g.id
     JOIN   sub           s  ON gd.disease_id = s.id
     JOIN   Variant       v  ON v.gene_id     = g.id
     JOIN   gnomAD        a  ON a.variant_id  = v.id
     WHERE  a.allele_freq < 0.01
     ORDER  BY g.symbol
     """),

    ("Q6", "Genes in any metabolic disease with any rare variant",
     """
     WITH RECURSIVE sub(id) AS (
         SELECT id FROM Disease WHERE name = 'Metabolic Disease'
         UNION
         SELECT d.id FROM Disease d JOIN sub s ON d.parent_id = s.id
     )
     SELECT DISTINCT g.symbol
     FROM   Gene          g
     JOIN   Gene_Disease  gd ON gd.gene_id    = g.id
     JOIN   sub           s  ON gd.disease_id = s.id
     JOIN   Variant       v  ON v.gene_id     = g.id
     JOIN   gnomAD        a  ON a.variant_id  = v.id
     WHERE  a.allele_freq < 0.01
     ORDER  BY g.symbol
     """),

    ("Q7", "MODY genes with any rare variant, depth <= 3 from 'Diabetes Mellitus'",
     """
     WITH RECURSIVE sub(id, depth) AS (
         SELECT id, 0 FROM Disease WHERE name = 'Diabetes Mellitus'
         UNION ALL
         SELECT d.id, s.depth + 1
         FROM   Disease d
         JOIN   sub     s ON d.parent_id = s.id
         WHERE  s.depth < 3
     )
     SELECT DISTINCT g.symbol, v.rsid, a.allele_freq
     FROM   Gene          g
     JOIN   Gene_Disease  gd ON gd.gene_id    = g.id
     JOIN   sub           s  ON gd.disease_id = s.id
     JOIN   Variant       v  ON v.gene_id     = g.id
     JOIN   gnomAD        a  ON a.variant_id  = v.id
     WHERE  a.allele_freq < 0.01
     ORDER  BY g.symbol
     """),
]

HIERQL_FORMS = {
    "Q1": "FIND variant.rsid, af.allele_freq FROM Gene AS gene JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gene.symbol = \"INS\" AND af.allele_freq < 0.001",
    "Q2": """FIND gene.symbol, COUNT(*) AS n_variants
FROM Gene AS gene
  JOIN Variant AS variant ON variant.gene_id = gene.id
GROUP BY gene.symbol
HAVING COUNT(*) > 0
ORDER BY n_variants DESC""",
    "Q3": "FIND d.name FROM Disease AS d WHERE d.id IN DISEASES UNDER \"Diabetes Mellitus\"",
    "Q4": "FIND d.name FROM Disease AS d WHERE d.id IN DISEASES ANCESTORS \"Type 2 Diabetes\"",
    "Q5": "FIND gene.symbol, variant.rsid, af.allele_freq FROM Gene AS gene JOIN Gene_Disease AS gd ON gd.gene_id = gene.id JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gd.disease_id IN DISEASES UNDER \"Diabetes Mellitus\" AND af.allele_freq < 0.01",
    "Q6": """FIND DISTINCT gene.symbol
FROM Gene AS gene
  JOIN Gene_Disease AS gd ON gd.gene_id = gene.id
  JOIN Variant AS variant ON variant.gene_id = gene.id
  JOIN gnomAD AS af ON af.variant_id = variant.id
WHERE gd.disease_id IN DISEASES UNDER "Metabolic Disease"
  AND af.allele_freq < 0.01""",
    "Q7": "FIND gene.symbol, variant.rsid, af.allele_freq FROM Gene AS gene JOIN Gene_Disease AS gd ON gd.gene_id = gene.id JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gd.disease_id IN DISEASES DEPTH \"Diabetes Mellitus\" 3 AND af.allele_freq < 0.01",
}

DATALOG_FORMS = {
    "Q1": """ins_gene(gid) :- gene(gid, "INS").
rare_variant(rsid, af) :-
  ins_gene(gid),
  variant(vid, gid, rsid),
  gnomad(_, vid, af),
  af < 0.001.
.output rare_variant""",
    "Q2": """gene_variant_count(sym, count : count vid) :-
  gene(gid, sym),
  variant(vid, gid, _).
.output gene_variant_count""",
    "Q3": """descendant(id, id) :- disease(id, "Diabetes Mellitus", _).
descendant(child, root) :- disease(child, _, parent), descendant(parent, root).
result(name) :- descendant(id, _), disease(id, name, _).
.output result""",
    "Q4": """ancestor(parent, child) :- disease(child, _, parent), parent != null.
ancestor(ancestor_id, child) :- disease(mid, _, ancestor_id), ancestor(mid, child).
result(name) :- disease(id, "Type 2 Diabetes", _), ancestor(anc, id), disease(anc, name, _).
.output result""",
    "Q5": """descendant(id, id) :- disease(id, "Diabetes Mellitus", _).
descendant(child, root) :- disease(child, _, parent), descendant(parent, root).
rare_gene(sym, rsid, af) :-
  descendant(did, _),
  gene_disease(gid, did),
  gene(gid, sym),
  variant(vid, gid, rsid),
  gnomad(_, vid, af),
  af < 0.01.
.output rare_gene""",
    "Q6": """descendant(id, id) :- disease(id, "Metabolic Disease", _).
descendant(child, root) :- disease(child, _, parent), descendant(parent, root).
rare_gene(sym) :-
  descendant(did, _),
  gene_disease(gid, did),
  gene(gid, sym),
  variant(vid, gid, _),
  gnomad(_, vid, af),
  af < 0.01.
.output rare_gene""",
    "Q7": """desc_depth(id, 0) :- disease(id, "Diabetes Mellitus", _).
desc_depth(child, depth + 1) :- disease(child, _, parent), desc_depth(parent, depth), depth < 3.
rare_gene(sym, rsid, af) :-
  desc_depth(did, _),
  gene_disease(gid, did),
  gene(gid, sym),
  variant(vid, gid, rsid),
  gnomad(_, vid, af),
  af < 0.01.
.output rare_gene""",
}

CYPHER_FORMS = {
    "Q1": """MATCH (g:Gene {symbol:'INS'})<-[:BELONGS_TO]-(v:Variant)-[:ANN]->(a:Annotation)
WHERE a.allele_freq < 0.001
RETURN v.rsid, a.allele_freq""",
    "Q2": """MATCH (g:Gene)<-[:BELONGS_TO]-(v:Variant)
RETURN g.symbol, count(v) AS n_variants
ORDER BY n_variants DESC""",
    "Q3": """MATCH (dm:Disease {name:'Diabetes Mellitus'})<-[:IS_A*0..]-(d:Disease)
RETURN d.name""",
    "Q4": """MATCH (d:Disease {name:'Type 2 Diabetes'})-[:IS_A*1..]->(a:Disease)
RETURN a.name""",
    "Q5": """MATCH (dm:Disease {name:'Diabetes Mellitus'})<-[:IS_A*0..]-(d:Disease)<-[:ASSOC]-(g:Gene)
MATCH (g)<-[:BELONGS_TO]-(v:Variant)-[:ANN]->(a:Annotation)
WHERE a.allele_freq < 0.01
RETURN g.symbol, v.rsid, a.allele_freq
ORDER BY g.symbol""",
    "Q6": """MATCH (m:Disease {name:'Metabolic Disease'})<-[:IS_A*0..]-(d:Disease)<-[:ASSOC]-(g:Gene)
MATCH (g)<-[:BELONGS_TO]-(v:Variant)-[:ANN]->(a:Annotation)
WHERE a.allele_freq < 0.01
RETURN DISTINCT g.symbol
ORDER BY g.symbol""",
    "Q7": """MATCH (dm:Disease {name:'Diabetes Mellitus'})<-[:IS_A*0..3]-(d:Disease)<-[:ASSOC]-(g:Gene)
MATCH (g)<-[:BELONGS_TO]-(v:Variant)-[:ANN]->(a:Annotation)
WHERE a.allele_freq < 0.01
RETURN DISTINCT g.symbol, v.rsid, a.allele_freq
ORDER BY g.symbol""",
}

SPARQL_FORMS = {
    "Q1": """PREFIX bio: <http://example.org/bio#>
SELECT ?rsid ?alleleFreq
WHERE {
  ?gene bio:symbol "INS" .
  ?variant bio:gene ?gene ;
           bio:rsid ?rsid .
  ?annot bio:variant ?variant ;
         bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.001)
}""",
    "Q2": """PREFIX bio: <http://example.org/bio#>
SELECT ?symbol (COUNT(?variant) AS ?nVariants)
WHERE {
  ?gene bio:symbol ?symbol .
  ?variant bio:gene ?gene .
}
GROUP BY ?symbol
ORDER BY DESC(?nVariants)""",
    "Q3": """PREFIX bio: <http://example.org/bio#>
SELECT ?name
WHERE {
  ?root bio:name "Diabetes Mellitus" .
  ?disease bio:subClassOf* ?root ;
           bio:name ?name .
}""",
    "Q4": """PREFIX bio: <http://example.org/bio#>
SELECT ?name
WHERE {
  ?disease bio:name "Type 2 Diabetes" .
  ?disease bio:subClassOf+ ?ancestor .
  ?ancestor bio:name ?name .
}""",
    "Q5": """PREFIX bio: <http://example.org/bio#>
SELECT ?symbol ?rsid ?alleleFreq
WHERE {
  ?root bio:name "Diabetes Mellitus" .
  ?disease bio:subClassOf* ?root .
  ?assoc bio:disease ?disease ;
         bio:gene ?gene .
  ?gene bio:symbol ?symbol .
  ?variant bio:gene ?gene ;
           bio:rsid ?rsid .
  ?annot bio:variant ?variant ;
         bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.01)
}
ORDER BY ?symbol""",
    "Q6": """PREFIX bio: <http://example.org/bio#>
SELECT DISTINCT ?symbol
WHERE {
  ?root bio:name "Metabolic Disease" .
  ?disease bio:subClassOf* ?root .
  ?assoc bio:disease ?disease ;
         bio:gene ?gene .
  ?gene bio:symbol ?symbol .
  ?variant bio:gene ?gene .
  ?annot bio:variant ?variant ;
         bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.01)
}
ORDER BY ?symbol""",
    "Q7": """PREFIX bio: <http://example.org/bio#>
SELECT DISTINCT ?symbol ?rsid ?alleleFreq
WHERE {
  ?root bio:name "Diabetes Mellitus" .
  ?disease bio:subClassOf ?x1 .
  ?x1 bio:subClassOf ?x2 .
  FILTER (?x1 = ?root || ?x2 = ?root || ?disease = ?root)
  ?assoc bio:disease ?disease ;
         bio:gene ?gene .
  ?gene bio:symbol ?symbol .
  ?variant bio:gene ?gene ;
           bio:rsid ?rsid .
  ?annot bio:variant ?variant ;
         bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.01)
}
ORDER BY ?symbol""",
}

# ── Real-data benchmark queries ──────────────────────────────────────────────

REAL_QUERIES = [
    ("Q1", "NLGN4Y chrY variants with AF < 0.001",
     """
     SELECT v.rsid, a.allele_freq
     FROM   Variant v
     JOIN   gnomAD a ON a.variant_id = v.id
     JOIN   Gene   g ON g.id = v.gene_id
     WHERE  g.symbol = 'NLGN4Y'
       AND  a.allele_freq < 0.001
     ORDER  BY a.allele_freq DESC, v.rsid
     """),
    ("Q2", "chrY genes with > 0 variants in the loaded subset",
     """
     SELECT g.symbol, COUNT(*) AS n_variants
     FROM   Gene g
     JOIN   Variant v ON v.gene_id = g.id
     GROUP  BY g.symbol
     HAVING COUNT(*) > 0
     ORDER  BY n_variants DESC, g.symbol
     """),
    ("Q3", "All MONDO descendants of 'gonadal disorder' via is_a",
     """
     WITH RECURSIVE sub(id, path) AS (
         SELECT id, ARRAY[id] AS path
         FROM   Disease
         WHERE  name = 'gonadal disorder'
         UNION ALL
         SELECT e.child_id, s.path || e.child_id
         FROM   DiseaseEdge e
         JOIN   sub s ON e.parent_id = s.id
         WHERE  e.rel = 'is_a'
           AND  NOT e.child_id = ANY(s.path)
     )
     SELECT d.name
     FROM   Disease d
     JOIN   (SELECT DISTINCT id FROM sub) s ON d.id = s.id
     ORDER  BY d.name
     """),
    ("Q4", "Ancestors of '46,XY complete gonadal dysgenesis' via is_a",
     """
     WITH RECURSIVE anc(id, path) AS (
         SELECT e.parent_id AS id, ARRAY[d.id, e.parent_id] AS path
         FROM   Disease d
         JOIN   DiseaseEdge e ON e.child_id = d.id
         WHERE  d.name = '46,XY complete gonadal dysgenesis'
           AND  e.rel = 'is_a'
         UNION ALL
         SELECT e.parent_id, a.path || e.parent_id
         FROM   DiseaseEdge e
         JOIN   anc a ON e.child_id = a.id
         WHERE  e.rel = 'is_a'
           AND  NOT e.parent_id = ANY(a.path)
     )
     SELECT d.name
     FROM   Disease d
     JOIN   (SELECT DISTINCT id FROM anc) a ON d.id = a.id
     ORDER  BY d.name
     """),
    ("Q5", "Rare chrY variants in genes linked to descendants of 'gonadal disorder'",
     """
     WITH RECURSIVE sub(id, path) AS (
         SELECT id, ARRAY[id] AS path
         FROM   Disease
         WHERE  name = 'gonadal disorder'
         UNION ALL
         SELECT e.child_id, s.path || e.child_id
         FROM   DiseaseEdge e
         JOIN   sub s ON e.parent_id = s.id
         WHERE  e.rel = 'is_a'
           AND  NOT e.child_id = ANY(s.path)
     )
     SELECT DISTINCT g.symbol, v.rsid, a.allele_freq
     FROM   Gene         g
     JOIN   Gene_Disease gd ON gd.gene_id = g.id
     JOIN   Variant      v  ON v.gene_id  = g.id
     JOIN   gnomAD       a  ON a.variant_id = v.id
     WHERE  gd.disease_id IN (SELECT DISTINCT id FROM sub)
       AND  a.allele_freq < 0.001
     ORDER  BY g.symbol, a.allele_freq DESC, v.rsid
     """),
    ("Q6", "chrY genes linked to descendants of 'reproductive system disorder' with AF < 0.001",
     """
     WITH RECURSIVE sub(id, path) AS (
         SELECT id, ARRAY[id] AS path
         FROM   Disease
         WHERE  name = 'reproductive system disorder'
         UNION ALL
         SELECT e.child_id, s.path || e.child_id
         FROM   DiseaseEdge e
         JOIN   sub s ON e.parent_id = s.id
         WHERE  e.rel = 'is_a'
           AND  NOT e.child_id = ANY(s.path)
     )
     SELECT DISTINCT g.symbol
     FROM   Gene         g
     JOIN   Gene_Disease gd ON gd.gene_id = g.id
     JOIN   Variant      v  ON v.gene_id  = g.id
     JOIN   gnomAD       a  ON a.variant_id = v.id
     WHERE  gd.disease_id IN (SELECT DISTINCT id FROM sub)
       AND  a.allele_freq < 0.001
     ORDER  BY g.symbol
     """),
    ("Q7", "chrY genes within DEPTH 3 of 'disorder of sexual differentiation' with AF < 0.001",
     """
     WITH RECURSIVE sub(id, depth, path) AS (
         SELECT id, 0 AS depth, ARRAY[id] AS path
         FROM   Disease
         WHERE  name = 'disorder of sexual differentiation'
         UNION ALL
         SELECT e.child_id, s.depth + 1, s.path || e.child_id
         FROM   DiseaseEdge e
         JOIN   sub s ON e.parent_id = s.id
         WHERE  e.rel = 'is_a'
           AND  s.depth < 3
           AND  NOT e.child_id = ANY(s.path)
     )
     SELECT DISTINCT g.symbol, v.rsid, a.allele_freq
     FROM   Gene         g
     JOIN   Gene_Disease gd ON gd.gene_id = g.id
     JOIN   Variant      v  ON v.gene_id  = g.id
     JOIN   gnomAD       a  ON a.variant_id = v.id
     WHERE  gd.disease_id IN (SELECT DISTINCT id FROM sub)
       AND  a.allele_freq < 0.001
     ORDER  BY g.symbol, a.allele_freq DESC, v.rsid
     """),
]

REAL_HIERQL_FORMS = {
    "Q1": "FIND variant.rsid, af.allele_freq FROM Gene AS gene JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gene.symbol = \"NLGN4Y\" AND af.allele_freq < 0.001 ORDER BY af.allele_freq DESC, variant.rsid",
    "Q2": """FIND gene.symbol, COUNT(*) AS n_variants
FROM Gene AS gene
  JOIN Variant AS variant ON variant.gene_id = gene.id
GROUP BY gene.symbol
HAVING COUNT(*) > 0
ORDER BY n_variants DESC, gene.symbol""",
    "Q3": "FIND d.name FROM Disease AS d WHERE d.id IN DISEASES UNDER \"gonadal disorder\" REL = \"is_a\" ORDER BY d.name",
    "Q4": "FIND d.name FROM Disease AS d WHERE d.id IN DISEASES ANCESTORS \"46,XY complete gonadal dysgenesis\" REL = \"is_a\" ORDER BY d.name",
    "Q5": "FIND DISTINCT gene.symbol, variant.rsid, af.allele_freq FROM Gene AS gene JOIN Gene_Disease AS gd ON gd.gene_id = gene.id JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gd.disease_id IN DISEASES UNDER \"gonadal disorder\" REL = \"is_a\" AND af.allele_freq < 0.001 ORDER BY gene.symbol, af.allele_freq DESC, variant.rsid",
    "Q6": "FIND DISTINCT gene.symbol FROM Gene AS gene JOIN Gene_Disease AS gd ON gd.gene_id = gene.id JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gd.disease_id IN DISEASES UNDER \"reproductive system disorder\" REL = \"is_a\" AND af.allele_freq < 0.001 ORDER BY gene.symbol",
    "Q7": "FIND DISTINCT gene.symbol, variant.rsid, af.allele_freq FROM Gene AS gene JOIN Gene_Disease AS gd ON gd.gene_id = gene.id JOIN Variant AS variant ON variant.gene_id = gene.id JOIN gnomAD AS af ON af.variant_id = variant.id WHERE gd.disease_id IN DISEASES DEPTH \"disorder of sexual differentiation\" 3 REL = \"is_a\" AND af.allele_freq < 0.001 ORDER BY gene.symbol, af.allele_freq DESC, variant.rsid",
}

# Real-data queries for alternative languages.
# Souffle fact schema:
#   disease(id:symbol, name:symbol)
#   disease_edge(parent:symbol, child:symbol, rel:symbol)
#   gene(id:number, symbol:symbol)
#   gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
#   variant(id:number, gene_id:number, rsid:symbol)
#   gnomad(id:number, variant_id:number, af:float)

REAL_DATALOG_FORMS = {
    "Q1": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl result(rsid:symbol, af:float)
result(rsid, af) :-
  gene(gid, "NLGN4Y"),
  variant(vid, gid, rsid),
  gnomad(_, vid, af),
  af < 0.001.
.output result(IO=stdout)""",

    "Q2": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl result(sym:symbol, n:number)
result(sym, n) :-
  gene(gid, sym),
  n = count : { variant(_, gid, _) },
  n > 0.
.output result""",

    "Q3": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl root(id:symbol)
root(id) :- disease(id, "gonadal disorder").
.decl descendant(id:symbol)
descendant(id) :- root(id).
descendant(c)  :- descendant(p), disease_edge(p, c, "is_a").
.decl result(name:symbol)
result(name) :- descendant(id), disease(id, name).
.output result(IO=stdout)""",

    "Q4": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl target(id:symbol)
target(id) :- disease(id, "46,XY complete gonadal dysgenesis").
.decl ancestor(id:symbol)
ancestor(p) :- target(c), disease_edge(p, c, "is_a").
ancestor(p) :- ancestor(c), disease_edge(p, c, "is_a").
.decl result(name:symbol)
result(name) :- ancestor(id), disease(id, name).
.output result(IO=stdout)""",

    "Q5": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl root(id:symbol)
root(id) :- disease(id, "gonadal disorder").
.decl descendant(id:symbol)
descendant(id) :- root(id).
descendant(c)  :- descendant(p), disease_edge(p, c, "is_a").
.decl result(sym:symbol, rsid:symbol, af:float)
result(sym, rsid, af) :-
  descendant(did),
  gene_disease(gid, did, _),
  gene(gid, sym),
  variant(vid, gid, rsid),
  gnomad(_, vid, af),
  af < 0.001.
.output result(IO=stdout)""",

    "Q6": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl root(id:symbol)
root(id) :- disease(id, "reproductive system disorder").
.decl descendant(id:symbol)
descendant(id) :- root(id).
descendant(c)  :- descendant(p), disease_edge(p, c, "is_a").
.decl result(sym:symbol)
result(sym) :-
  descendant(did),
  gene_disease(gid, did, _),
  gene(gid, sym),
  variant(vid, gid, _),
  gnomad(_, vid, af),
  af < 0.001.
.output result(IO=stdout)""",

    "Q7": """\
.decl disease(id:symbol, name:symbol)
.decl disease_edge(parent:symbol, child:symbol, rel:symbol)
.decl gene(id:number, symbol:symbol)
.decl gene_disease(gene_id:number, disease_id:symbol, classification:symbol)
.decl variant(id:number, gene_id:number, rsid:symbol)
.decl gnomad(id:number, variant_id:number, af:float)
.input disease .input disease_edge .input gene .input gene_disease .input variant .input gnomad

.decl desc_depth(id:symbol, d:number)
desc_depth(id, 0) :- disease(id, "disorder of sexual differentiation").
desc_depth(c, d+1) :- desc_depth(p, d), disease_edge(p, c, "is_a"), d < 3.
.decl result(sym:symbol, rsid:symbol, af:float)
result(sym, rsid, af) :-
  desc_depth(did, _),
  gene_disease(gid, did, _),
  gene(gid, sym),
  variant(vid, gid, rsid),
  gnomad(_, vid, af),
  af < 0.001.
.output result(IO=stdout)""",
}

# Neo4j property-graph schema:
#   (:Disease {id, name})-[:IS_A]->(:Disease)
#   (:Gene {id, symbol})-[:ASSOC {classification}]->(:Disease)
#   (:Variant {id, rsid})-[:BELONGS_TO]->(:Gene)
#   (:Variant)-[:HAS_AF {allele_freq}]->(:AF)   (AF node carries allele_freq)
# Simpler: store allele_freq as property on a (:gnomAD {allele_freq}) node
#   (:gnomAD {allele_freq})<-[:IN_GNOMAD]-(:Variant)

REAL_CYPHER_FORMS = {
    "Q1": """\
MATCH (g:Gene {symbol:'NLGN4Y'})<-[:BELONGS_TO]-(v:Variant)-[:IN_GNOMAD]->(a:gnomAD)
WHERE a.allele_freq < 0.001
RETURN v.rsid AS rsid, a.allele_freq AS allele_freq
ORDER BY a.allele_freq DESC, v.rsid""",

    "Q2": """\
MATCH (g:Gene)<-[:BELONGS_TO]-(v:Variant)
RETURN g.symbol AS symbol, count(v) AS n_variants
ORDER BY n_variants DESC, g.symbol""",

    "Q3": """\
MATCH (:Disease {name:'gonadal disorder'})<-[:IS_A*0..]-(d:Disease)
RETURN d.name AS name
ORDER BY name""",

    "Q4": """\
MATCH (:Disease {name:'46,XY complete gonadal dysgenesis'})-[:IS_A*1..]->(a:Disease)
RETURN a.name AS name
ORDER BY a.name""",

    "Q5": """\
MATCH (:Disease {name:'gonadal disorder'})<-[:IS_A*0..]-(:Disease)<-[:ASSOC]-(g:Gene)
MATCH (g)<-[:BELONGS_TO]-(v:Variant)-[:IN_GNOMAD]->(a:gnomAD)
WHERE a.allele_freq < 0.001
RETURN DISTINCT g.symbol AS symbol, v.rsid AS rsid, a.allele_freq AS allele_freq
ORDER BY g.symbol, a.allele_freq DESC, v.rsid""",

    "Q6": """\
MATCH (:Disease {name:'reproductive system disorder'})<-[:IS_A*0..]-(:Disease)<-[:ASSOC]-(g:Gene)
MATCH (g)<-[:BELONGS_TO]-(v:Variant)-[:IN_GNOMAD]->(a:gnomAD)
WHERE a.allele_freq < 0.001
RETURN DISTINCT g.symbol AS symbol
ORDER BY g.symbol""",

    "Q7": """\
MATCH (:Disease {name:'disorder of sexual differentiation'})<-[:IS_A*0..3]-(:Disease)<-[:ASSOC]-(g:Gene)
MATCH (g)<-[:BELONGS_TO]-(v:Variant)-[:IN_GNOMAD]->(a:gnomAD)
WHERE a.allele_freq < 0.001
RETURN DISTINCT g.symbol AS symbol, v.rsid AS rsid, a.allele_freq AS allele_freq
ORDER BY g.symbol, a.allele_freq DESC, v.rsid""",
}

# SPARQL queries run against an in-memory rdflib graph.
# Prefix bio: <http://example.org/bio#>
# Triples: disease uri bio:name "...", disease_edge: child bio:isA parent,
#          gene bio:symbol "...", gene_disease: gene bio:associatedWith disease,
#          variant bio:belongsTo gene, variant bio:rsid "...",
#          variant bio:inGnomad gnomad_node, gnomad_node bio:alleleFreq "..."^^xsd:double

REAL_SPARQL_FORMS = {
    "Q1": """\
PREFIX bio: <http://example.org/bio#>
SELECT ?rsid ?alleleFreq
WHERE {
  ?gene bio:symbol "NLGN4Y" .
  ?variant bio:belongsTo ?gene ;
           bio:rsid ?rsid ;
           bio:inGnomad ?af .
  ?af bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.001)
}
ORDER BY DESC(?alleleFreq) ?rsid""",

    "Q2": """\
PREFIX bio: <http://example.org/bio#>
SELECT ?symbol (COUNT(?variant) AS ?n_variants)
WHERE {
  ?gene bio:symbol ?symbol .
  ?variant bio:belongsTo ?gene .
}
GROUP BY ?symbol
HAVING (COUNT(?variant) > 0)
ORDER BY DESC(?n_variants) ?symbol""",

    "Q3": """\
PREFIX bio: <http://example.org/bio#>
SELECT ?name
WHERE {
  ?root bio:name "gonadal disorder" .
  ?disease bio:isA* ?root ;
           bio:name ?name .
}
ORDER BY ?name""",

    "Q4": """\
PREFIX bio: <http://example.org/bio#>
SELECT ?name
WHERE {
  ?target bio:name "46,XY complete gonadal dysgenesis" .
  ?target bio:isA+ ?ancestor .
  ?ancestor bio:name ?name .
}
ORDER BY ?name""",

    "Q5": """\
PREFIX bio: <http://example.org/bio#>
SELECT DISTINCT ?symbol ?rsid ?alleleFreq
WHERE {
  ?root bio:name "gonadal disorder" .
  ?disease bio:isA* ?root .
  ?gene bio:associatedWith ?disease ;
        bio:symbol ?symbol .
  ?variant bio:belongsTo ?gene ;
           bio:rsid ?rsid ;
           bio:inGnomad ?af .
  ?af bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.001)
}
ORDER BY ?symbol DESC(?alleleFreq) ?rsid""",

    "Q6": """\
PREFIX bio: <http://example.org/bio#>
SELECT DISTINCT ?symbol
WHERE {
  ?root bio:name "reproductive system disorder" .
  ?disease bio:isA* ?root .
  ?gene bio:associatedWith ?disease ;
        bio:symbol ?symbol .
  ?variant bio:belongsTo ?gene ;
           bio:inGnomad ?af .
  ?af bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.001)
}
ORDER BY ?symbol""",

    # rdflib does not support bounded property paths (e.g. bio:isA{0,3}).
    # All descendants of 'disorder of sexual differentiation' in this dataset
    # fall within depth 3, so bio:isA* gives the correct result here.
    "Q7": """\
PREFIX bio: <http://example.org/bio#>
SELECT DISTINCT ?symbol ?rsid ?alleleFreq
WHERE {
  ?root bio:name "disorder of sexual differentiation" .
  ?disease bio:isA* ?root .
  ?gene bio:associatedWith ?disease ;
        bio:symbol ?symbol .
  ?variant bio:belongsTo ?gene ;
           bio:rsid ?rsid ;
           bio:inGnomad ?af .
  ?af bio:alleleFreq ?alleleFreq .
  FILTER (?alleleFreq < 0.001)
}
ORDER BY ?symbol DESC(?alleleFreq) ?rsid""",
}

# ── Helpers ───────────────────────────────────────────────────────────────────

CLASSIFICATION_RANK = {
    "Definitive": 5,
    "Strong": 4,
    "Moderate": 3,
    "Limited": 2,
    "Supportive": 1,
    "Disputed": 0,
    "Refuted": 0,
    "No Known Disease Relationship": 0,
    "Animal Model Only": 0,
}


def count_tokens(src: str) -> int:
    """Whitespace-split fallback (not used directly; use language-specific counters)."""
    return len(src.split())


def _tokenize(src: str, pattern: str) -> list:
    """Return all non-whitespace tokens matched by `pattern`."""
    return re.findall(pattern, src)


def count_tokens_sql(src: str) -> int:
    # Keywords, identifiers, dotted names, string literals, numbers,
    # multi-char operators (<=, >=, <>, !=, ||), single-char punct.
    # Does not split on whitespace alone — each lexical unit counts once.
    pat = (
        r"'[^']*'"                            # string literal
        r'|"[^"]*"'                           # double-quoted identifier
        r'|\d+\.\d+|\d+'                      # number
        r'|<>|!=|<=|>=|\|\|'                  # multi-char operators
        r'|[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*'  # identifier or dotted name
        r'|[^\s\w]'                           # any remaining single non-whitespace char
    )
    return len(_tokenize(src, pat))


def count_tokens_hierql(src: str) -> int:
    # Same as SQL but HierQL has no double-quoted identifiers (uses " for strings)
    # and adds -> operator. Reuse SQL tokenizer — structure is identical.
    return count_tokens_sql(src)


def count_tokens_datalog(src: str) -> int:
    # Strip .decl and .input lines — these are shared schema boilerplate that
    # appears identically in every query and does not represent query logic.
    # .output is kept as it is part of the query specification.
    query_lines = [
        l for l in src.splitlines()
        if not l.lstrip().startswith('.decl') and not l.lstrip().startswith('.input')
    ]
    src = '\n'.join(query_lines)
    # Atoms: identifiers, quoted strings, numbers, :- operator, _ wildcard,
    # punctuation (,  .  (  )  =  <  >  !=).
    pat = (
        r'"[^"]*"'             # string literal
        r'|\d+\.\d+|\d+'       # number
        r'|:-'                 # rule neck
        r'|!='                 # not-equal
        r'|[A-Za-z_]\w*'       # identifier / variable
        r'|[^\s\w]'            # punct: ( ) , . : = < > ! _
    )
    return len(_tokenize(src, pat))


def count_tokens_cypher(src: str) -> int:
    # Cypher packs a lot into path patterns with almost no whitespace, e.g.:
    #   (g:Gene {symbol:'INS'})<-[:BELONGS_TO]-(v:Variant)
    # Tokens: identifiers/labels, quoted strings, numbers, relationship-type
    # brackets ([ ]), node parens (( )), arrows (<-, ->, --), * range specs
    # (e.g. *0.., *1..3), property operators ({  }  :  ,  =  .  <  >  !=).
    pat = (
        r"'[^']*'"             # single-quoted string
        r'|"[^"]*"'            # double-quoted string
        r'|\*\d*\.\.\d*'       # path quantifier: *0.., *1..3, *..
        r'|\d+\.\d+|\d+'       # number
        r'|<-|->|--'           # relationship arrows
        r'|<>|!=|<=|>='        # comparison operators
        r'|[A-Za-z_]\w*'       # keyword or identifier
        r'|[^\s\w]'            # remaining punct: ( ) [ ] { } : , . = < > ! |
    )
    return len(_tokenize(src, pat))


def count_tokens_sparql(src: str) -> int:
    # SPARQL has prefixed names (bio:symbol), ?variables, path operators (* + |),
    # URI literals (<...>), string literals, numbers, and standard keywords.
    pat = (
        r'<[^>]*>'                      # URI literal <http://...>
        r'|"[^"]*"'                     # string literal
        r'|\d+\.\d+|\d+'                # number
        r'|\?[A-Za-z_]\w*'              # variable ?foo
        r'|[A-Za-z_]\w*:[A-Za-z_]\w*'  # prefixed name bio:symbol
        r'|[A-Za-z_]\w*'                # keyword or bare identifier
        r'|[^\s\w]'                     # punct: { } ( ) . ; , = < > ! * + | ^
    )
    return len(_tokenize(src, pat))


def count_lines(sql: str) -> int:
    return sum(1 for l in sql.strip().splitlines() if l.strip())


def chunked(rows, size: int):
    for i in range(0, len(rows), size):
        yield rows[i:i + size]


def run_query(cur, sql: str, repeats: int = 3) -> tuple[list, float]:
    rows = None
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        cur.execute(sql)
        rows = cur.fetchall()
        times.append((time.perf_counter() - t0) * 1000)
    return rows, sum(times) / len(times)


def print_rows_preview(rows: list, limit: int) -> None:
    for row in rows[:limit]:
        print("  ", row)
    if len(rows) > limit:
        print(f"  ... {len(rows) - limit} more row(s)")


def mondo_id_from_uri(uri: str) -> str | None:
    if "MONDO_" not in uri:
        return None
    return "MONDO:" + uri.rsplit("MONDO_", 1)[1]


def load_real_dataset(cur, execute_values):
    print("Loading MONDO graph...", flush=True)
    with open(REAL_MONDO_PATH) as f:
        graph = json.load(f)["graphs"][0]

    disease_rows = []
    for node in graph["nodes"]:
        disease_id = mondo_id_from_uri(node.get("id", ""))
        if disease_id:
            disease_rows.append((disease_id, node.get("lbl") or disease_id))

    edge_rows = []
    for edge in graph["edges"]:
        sub = mondo_id_from_uri(edge.get("sub", ""))
        obj = mondo_id_from_uri(edge.get("obj", ""))
        if sub and obj:
            pred = edge.get("pred")
            rel = pred if pred == "is_a" else pred.rsplit("/", 1)[-1]
            edge_rows.append((obj, sub, rel))

    cur.execute(SETUP_SQL_REAL)
    for batch in chunked(disease_rows, 2000):
        execute_values(cur, "INSERT INTO Disease (id, name) VALUES %s", batch)
    for batch in chunked(edge_rows, 4000):
        execute_values(cur, "INSERT INTO DiseaseEdge (parent_id, child_id, rel) VALUES %s", batch)

    print("Scanning GenCC + chrY gnomAD overlap...", flush=True)
    gencc_genes = set()
    with open(REAL_GENCC_PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["gene_symbol"].strip()
            disease = row["disease_curie"].strip()
            if gene and disease.startswith("MONDO:"):
                gencc_genes.add(gene)

    variants_by_gene = defaultdict(list)
    with gzip.open(REAL_GNOMAD_PATH, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, pos, vid, ref, alt, qual, filt, info = line.rstrip("\n").split("\t")[:8]
            af = None
            vep_value = None
            for field in info.split(";"):
                if field.startswith("AF="):
                    try:
                        af = float(field[3:].split(",")[0])
                    except ValueError:
                        af = None
                elif field.startswith("vep="):
                    vep_value = field[4:]
            if af is None or vep_value is None:
                continue

            genes = set()
            for ann in vep_value.split(","):
                parts = ann.split("|")
                if len(parts) > 3 and parts[3] in gencc_genes:
                    genes.add(parts[3])
            if not genes:
                continue

            rsid = vid if vid != "." else f"{chrom}:{pos}:{ref}>{alt}"
            pos_i = int(pos)
            for gene in genes:
                variants_by_gene[gene].append((rsid, chrom, pos_i, ref, alt, af, filt))

    kept_genes = sorted(variants_by_gene)
    print(f"Found {len(kept_genes)} overlapping chrY gene(s): {', '.join(kept_genes)}", flush=True)

    print("Loading genes and deduplicated GenCC assertions...", flush=True)
    execute_values(cur, "INSERT INTO Gene (symbol) VALUES %s", [(g,) for g in kept_genes])
    cur.execute("SELECT id, symbol FROM Gene")
    gene_id_by_symbol = {symbol: gene_id for gene_id, symbol in cur.fetchall()}

    best_gene_disease = {}
    with open(REAL_GENCC_PATH, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene = row["gene_symbol"].strip()
            disease = row["disease_curie"].strip()
            if gene not in gene_id_by_symbol or not disease.startswith("MONDO:"):
                continue
            key = (gene, disease)
            classification = row["classification_title"].strip()
            rank = CLASSIFICATION_RANK.get(classification, -1)
            prev = best_gene_disease.get(key)
            if prev is None or rank > prev[0]:
                best_gene_disease[key] = (rank, classification)

    gene_disease_rows = [
        (gene_id_by_symbol[gene], disease, classification)
        for (gene, disease), (_, classification) in sorted(best_gene_disease.items())
    ]
    if gene_disease_rows:
        execute_values(
            cur,
            "INSERT INTO Gene_Disease (gene_id, disease_id, classification) VALUES %s",
            gene_disease_rows,
        )

    print("Loading chrY variants for overlapping genes...", flush=True)
    variant_rows = []
    for gene in kept_genes:
        gene_id = gene_id_by_symbol[gene]
        for rsid, chrom, pos, ref, alt, af, filt in variants_by_gene[gene]:
            variant_rows.append((gene_id, rsid, chrom, pos, ref, alt, af, filt))

    variant_id_rows = []
    insert_sql = (
        "INSERT INTO Variant (gene_id, rsid, chrom, pos, ref, alt) VALUES %s "
        "RETURNING id"
    )
    for batch in chunked([(r[0], r[1], r[2], r[3], r[4], r[5]) for r in variant_rows], 1000):
        returned = execute_values(cur, insert_sql, batch, fetch=True)
        variant_id_rows.extend([row[0] for row in returned])

    gnomad_rows = [
        (variant_id, variant_rows[i][6], variant_rows[i][7])
        for i, variant_id in enumerate(variant_id_rows)
    ]
    for batch in chunked(gnomad_rows, 2000):
        execute_values(
            cur,
            "INSERT INTO gnomAD (variant_id, allele_freq, filt) VALUES %s",
            batch,
        )

    print(
        f"Loaded {len(disease_rows)} diseases, {len(edge_rows)} edges, "
        f"{len(kept_genes)} genes, {len(gene_disease_rows)} gene-disease pairs, "
        f"and {len(variant_rows)} chrY gene-variant rows.\n",
        flush=True,
    )


# ── Result normalization ──────────────────────────────────────────────────────

def _coerce(v):
    """Normalize floats to 4 significant figures for cross-language comparison.
    Souffle uses 32-bit floats (float32 has ~7 decimal sig-figs), which can
    differ from PostgreSQL's float64 at the 5th–7th digit. 4 sig-figs is safely
    below the float32 precision boundary while keeping distinct AF values distinct."""
    if isinstance(v, float):
        return float(f"{v:.4g}")
    if isinstance(v, str):
        try:
            return float(f"{float(v):.4g}")
        except ValueError:
            return v
    return v


def normalize_rows(rows) -> frozenset:
    """Convert any iterable of rows into a frozenset of tuples for comparison."""
    return frozenset(tuple(_coerce(c) for c in row) for row in rows)


# ── Souffle (Datalog) backend ─────────────────────────────────────────────────

def write_souffle_facts(cur, tmpdir: str):
    """Dump the PostgreSQL tables into Souffle tab-separated fact files."""

    # disease(id, name)
    cur.execute("SELECT id, name FROM Disease")
    with open(os.path.join(tmpdir, "disease.facts"), "w") as f:
        for row in cur.fetchall():
            f.write(f"{row[0]}\t{row[1]}\n")

    # disease_edge(parent_id, child_id, rel)
    cur.execute("SELECT parent_id, child_id, rel FROM DiseaseEdge")
    with open(os.path.join(tmpdir, "disease_edge.facts"), "w") as f:
        for row in cur.fetchall():
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")

    # gene(id, symbol)
    cur.execute("SELECT id, symbol FROM Gene")
    with open(os.path.join(tmpdir, "gene.facts"), "w") as f:
        for row in cur.fetchall():
            f.write(f"{row[0]}\t{row[1]}\n")

    # gene_disease(gene_id, disease_id, classification)
    cur.execute("SELECT gene_id, disease_id, COALESCE(classification, '') FROM Gene_Disease")
    with open(os.path.join(tmpdir, "gene_disease.facts"), "w") as f:
        for row in cur.fetchall():
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")

    # variant(id, gene_id, rsid)
    cur.execute("SELECT id, gene_id, COALESCE(rsid, '') FROM Variant")
    with open(os.path.join(tmpdir, "variant.facts"), "w") as f:
        for row in cur.fetchall():
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")

    # gnomad(id, variant_id, allele_freq)
    # Write allele_freq pre-rounded to 4 sig figs so that Souffle's float32
    # conversion of the already-rounded value is stable and matches _coerce().
    cur.execute("SELECT id, variant_id, allele_freq FROM gnomAD")
    with open(os.path.join(tmpdir, "gnomad.facts"), "w") as f:
        for row in cur.fetchall():
            af = float(f"{row[2]:.4g}")
            f.write(f"{row[0]}\t{row[1]}\t{af}\n")


def run_souffle(program: str, tmpdir: str, repeats: int) -> tuple[list, float]:
    """Write program to tmpdir, run souffle with file output, parse result.csv."""
    dl_path = os.path.join(tmpdir, "query.dl")
    out_csv  = os.path.join(tmpdir, "result.csv")

    # Use file output (no headers) instead of stdout (which adds 3 header lines).
    # Strip any residual (IO=stdout) directives from the program.
    program = program.replace("(IO=stdout)", "")

    with open(dl_path, "w") as f:
        f.write(program)

    times = []
    rows  = None
    for _ in range(repeats):
        # Remove previous output file so a failed run is not silently reused.
        if os.path.exists(out_csv):
            os.remove(out_csv)
        t0 = time.perf_counter()
        result = subprocess.run(
            ["souffle", "-F", tmpdir, "-D", tmpdir, dl_path],
            capture_output=True, text=True
        )
        elapsed = (time.perf_counter() - t0) * 1000
        times.append(elapsed)
        if result.returncode != 0:
            raise RuntimeError(f"Souffle error:\n{result.stderr}")
        # Parse the tab-separated output file (no headers).
        out_rows = []
        if os.path.exists(out_csv):
            with open(out_csv) as fh:
                for line in fh:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    parts = line.split("\t")
                    coerced = []
                    for p in parts:
                        try:
                            coerced.append(float(p))
                        except ValueError:
                            coerced.append(p)
                    out_rows.append(tuple(coerced))
        rows = out_rows

    return rows, sum(times) / len(times)


# ── Neo4j (Cypher) backend ────────────────────────────────────────────────────

def _neo4j_batch(session, query: str, rows: list, batch_size: int = 500) -> None:
    """Send rows to Neo4j in batches using UNWIND to avoid per-row round-trips."""
    for i in range(0, len(rows), batch_size):
        session.run(query, rows=rows[i:i + batch_size])


def setup_neo4j(cur) -> None:
    """Load the relational data from PostgreSQL into Neo4j using batched UNWIND."""
    try:
        from neo4j import GraphDatabase
    except ImportError as e:
        raise SystemExit("Install neo4j: pip install neo4j") from e

    driver = GraphDatabase.driver(NEO4J_URL, auth=NEO4J_AUTH)
    with driver.session() as s:
        # Wipe existing data
        s.run("MATCH (n) DETACH DELETE n")

        # Indexes first (before data so they're ready when edges are created)
        s.run("CREATE INDEX disease_id   IF NOT EXISTS FOR (d:Disease) ON (d.id)")
        s.run("CREATE INDEX disease_name IF NOT EXISTS FOR (d:Disease) ON (d.name)")
        s.run("CREATE INDEX gene_id      IF NOT EXISTS FOR (g:Gene)    ON (g.id)")
        s.run("CREATE INDEX gene_symbol  IF NOT EXISTS FOR (g:Gene)    ON (g.symbol)")
        s.run("CREATE INDEX variant_id   IF NOT EXISTS FOR (v:Variant) ON (v.id)")

        # Disease nodes
        cur.execute("SELECT id, name FROM Disease")
        rows = [{"id": str(r[0]), "name": r[1]} for r in cur.fetchall()]
        _neo4j_batch(s,
            "UNWIND $rows AS r CREATE (:Disease {id: r.id, name: r.name})",
            rows)
        print(f"    {len(rows)} disease nodes", flush=True)

        # IS_A edges
        cur.execute("SELECT parent_id, child_id FROM DiseaseEdge WHERE rel = 'is_a'")
        rows = [{"p": str(r[0]), "c": str(r[1])} for r in cur.fetchall()]
        _neo4j_batch(s,
            "UNWIND $rows AS r "
            "MATCH (p:Disease {id: r.p}), (c:Disease {id: r.c}) "
            "CREATE (c)-[:IS_A]->(p)",
            rows)
        print(f"    {len(rows)} IS_A edges", flush=True)

        # Gene nodes
        cur.execute("SELECT id, symbol FROM Gene")
        rows = [{"id": r[0], "sym": r[1]} for r in cur.fetchall()]
        _neo4j_batch(s,
            "UNWIND $rows AS r CREATE (:Gene {id: r.id, symbol: r.sym})",
            rows)

        # ASSOC edges (Gene)-[:ASSOC]->(Disease)
        cur.execute("SELECT gene_id, disease_id FROM Gene_Disease")
        rows = [{"g": r[0], "d": str(r[1])} for r in cur.fetchall()]
        _neo4j_batch(s,
            "UNWIND $rows AS r "
            "MATCH (g:Gene {id: r.g}), (d:Disease {id: r.d}) "
            "CREATE (g)-[:ASSOC]->(d)",
            rows)

        # Variant nodes + BELONGS_TO edges
        cur.execute("SELECT id, gene_id, COALESCE(rsid,'') FROM Variant")
        rows = [{"id": r[0], "g": r[1], "rsid": r[2]} for r in cur.fetchall()]
        _neo4j_batch(s,
            "UNWIND $rows AS r "
            "MATCH (g:Gene {id: r.g}) "
            "CREATE (:Variant {id: r.id, rsid: r.rsid})-[:BELONGS_TO]->(g)",
            rows)
        print(f"    {len(rows)} variant nodes", flush=True)

        # gnomAD nodes + IN_GNOMAD edges
        cur.execute("SELECT id, variant_id, allele_freq FROM gnomAD")
        rows = [{"id": r[0], "v": r[1], "af": r[2]} for r in cur.fetchall()]
        _neo4j_batch(s,
            "UNWIND $rows AS r "
            "MATCH (v:Variant {id: r.v}) "
            "CREATE (:gnomAD {id: r.id, allele_freq: r.af})<-[:IN_GNOMAD]-(v)",
            rows)
        print(f"    {len(rows)} gnomAD nodes", flush=True)

    driver.close()
    print("  Neo4j loaded.", flush=True)


def run_cypher(query: str, repeats: int) -> tuple[list, float]:
    from neo4j import GraphDatabase
    driver = GraphDatabase.driver(NEO4J_URL, auth=NEO4J_AUTH)
    times = []
    rows = None
    with driver.session() as s:
        for _ in range(repeats):
            t0 = time.perf_counter()
            result = s.run(query)
            rows = [tuple(r.values()) for r in result]
            times.append((time.perf_counter() - t0) * 1000)
    driver.close()
    return rows, sum(times) / len(times)


# ── rdflib (SPARQL) backend ───────────────────────────────────────────────────

_BIO = "http://example.org/bio#"


def setup_rdflib(cur):
    """Build an rdflib Graph from the PostgreSQL data."""
    from rdflib import Graph, URIRef, Literal, XSD
    from rdflib.namespace import RDF

    g = Graph()
    BIO = _BIO

    def uri(kind, pk):
        return URIRef(f"{BIO}{kind}/{pk}")

    # Disease nodes + name triples
    cur.execute("SELECT id, name FROM Disease")
    for did, name in cur.fetchall():
        d_uri = uri("disease", did)
        g.add((d_uri, RDF.type, URIRef(f"{BIO}Disease")))
        g.add((d_uri, URIRef(f"{BIO}name"), Literal(name)))

    # isA edges (child bio:isA parent)
    cur.execute("SELECT parent_id, child_id FROM DiseaseEdge WHERE rel = 'is_a'")
    for parent, child in cur.fetchall():
        g.add((uri("disease", child), URIRef(f"{BIO}isA"), uri("disease", parent)))

    # Gene nodes
    cur.execute("SELECT id, symbol FROM Gene")
    for gid, sym in cur.fetchall():
        g_uri = uri("gene", gid)
        g.add((g_uri, RDF.type, URIRef(f"{BIO}Gene")))
        g.add((g_uri, URIRef(f"{BIO}symbol"), Literal(sym)))

    # Gene-Disease associations (gene bio:associatedWith disease)
    cur.execute("SELECT gene_id, disease_id FROM Gene_Disease")
    for gid, did in cur.fetchall():
        g.add((uri("gene", gid), URIRef(f"{BIO}associatedWith"), uri("disease", did)))

    # Variant nodes
    cur.execute("SELECT id, gene_id, COALESCE(rsid,'') FROM Variant")
    for vid, gid, rsid in cur.fetchall():
        v_uri = uri("variant", vid)
        g.add((v_uri, RDF.type, URIRef(f"{BIO}Variant")))
        g.add((v_uri, URIRef(f"{BIO}rsid"), Literal(rsid)))
        g.add((v_uri, URIRef(f"{BIO}belongsTo"), uri("gene", gid)))

    # gnomAD nodes
    cur.execute("SELECT id, variant_id, allele_freq FROM gnomAD")
    for aid, vid, af in cur.fetchall():
        a_uri = uri("gnomad", aid)
        g.add((a_uri, RDF.type, URIRef(f"{BIO}gnomAD")))
        g.add((a_uri, URIRef(f"{BIO}alleleFreq"), Literal(af, datatype=XSD.double)))
        g.add((uri("variant", vid), URIRef(f"{BIO}inGnomad"), a_uri))

    return g


def run_sparql(graph, query: str, repeats: int) -> tuple[list, float]:
    times = []
    rows = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = graph.query(query)
        rows = [tuple(str(v) for v in r) for r in result]
        times.append((time.perf_counter() - t0) * 1000)
    return rows, sum(times) / len(times)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", choices=["sample", "real"], default="sample")
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--preview-rows", type=int, default=10)
    return p.parse_args()


def get_dataset_spec(dataset: str):
    if dataset == "sample":
        return SAMPLE_QUERIES, HIERQL_FORMS, DATALOG_FORMS, CYPHER_FORMS, SPARQL_FORMS
    return REAL_QUERIES, REAL_HIERQL_FORMS, REAL_DATALOG_FORMS, REAL_CYPHER_FORMS, REAL_SPARQL_FORMS


def _ms(v: float) -> str:
    return f"{v:.2f}" if not math.isnan(v) else "—"


def main():
    args = parse_args()

    try:
        import psycopg2
        from psycopg2.extras import execute_values
    except ImportError as e:
        raise SystemExit("Install psycopg2-binary to run the benchmark script.") from e

    conn = psycopg2.connect(DB_URL)
    conn.autocommit = True
    cur = conn.cursor()

    if args.dataset == "sample":
        print("Setting up sample schema and data...", flush=True)
        cur.execute(SETUP_SQL_SAMPLE)
        print("Done.\n")
    else:
        load_real_dataset(cur, execute_values)

    queries, hierql_forms, datalog_forms, cypher_forms, sparql_forms = get_dataset_spec(args.dataset)

    # ── Set up alternative backends (real dataset only) ───────────────────────
    neo4j_loaded = False
    rdf_graph = None
    souffle_tmpdir = None

    if args.dataset == "real" and (cypher_forms or sparql_forms or datalog_forms):
        # Neo4j
        if cypher_forms:
            try:
                print("Loading Neo4j...", flush=True)
                setup_neo4j(cur)
                neo4j_loaded = True
            except Exception as e:
                print(f"  Neo4j setup failed ({e}); Cypher timings will be skipped.", flush=True)

        # rdflib in-memory graph for SPARQL
        if sparql_forms:
            print("Building rdflib graph for SPARQL...", flush=True)
            rdf_graph = setup_rdflib(cur)
            print(f"  {len(rdf_graph)} triples loaded.", flush=True)

        # Souffle fact files
        if datalog_forms:
            souffle_tmpdir = tempfile.mkdtemp(prefix="hierql_souffle_")
            print(f"Writing Souffle fact files to {souffle_tmpdir}...", flush=True)
            write_souffle_facts(cur, souffle_tmpdir)

    timing_rows = []
    concise_rows = []
    correctness_rows = []

    for qid, desc, sql in queries:
        print(f"── {qid}: {desc}", flush=True)

        # SQL (reference)
        sql_rows, sql_ms = run_query(cur, sql, repeats=args.repeats)
        sql_set = normalize_rows(sql_rows)

        # HierQL (transpiled → SQL)
        hierql = hierql_forms.get(qid, "")
        if hierql:
            t_sql = transpile(hierql)
            hql_rows, hql_ms = run_query(cur, t_sql, repeats=args.repeats)
            hql_set = normalize_rows(hql_rows)
        else:
            hql_ms = math.nan
            hql_set = None

        # Datalog / Souffle
        dl_ms = math.nan
        dl_set = None
        if souffle_tmpdir and qid in datalog_forms:
            try:
                dl_rows, dl_ms = run_souffle(datalog_forms[qid], souffle_tmpdir, args.repeats)
                dl_set = normalize_rows(dl_rows)
            except Exception as e:
                print(f"  Souffle {qid} error: {e}", flush=True)

        # Cypher / Neo4j
        cy_ms = math.nan
        cy_set = None
        if neo4j_loaded and qid in cypher_forms:
            try:
                cy_rows, cy_ms = run_cypher(cypher_forms[qid], args.repeats)
                cy_set = normalize_rows(cy_rows)
            except Exception as e:
                print(f"  Cypher {qid} error: {e}", flush=True)

        # SPARQL / rdflib
        sp_ms = math.nan
        sp_set = None
        if rdf_graph is not None and qid in sparql_forms:
            try:
                sp_rows, sp_ms = run_sparql(rdf_graph, sparql_forms[qid], args.repeats)
                sp_set = normalize_rows(sp_rows)
            except Exception as e:
                print(f"  SPARQL {qid} error: {e}", flush=True)

        # Correctness: compare every available result set against SQL reference
        def match(s, label):
            if s is None:
                return "skip"
            return "OK" if s == sql_set else f"MISMATCH(sql={len(sql_set)},{label}={len(s)})"

        correctness_rows.append((
            qid,
            match(hql_set,  "hql"),
            match(dl_set,   "dl"),
            match(cy_set,   "cy"),
            match(sp_set,   "sp"),
        ))

        timing_rows.append((
            qid, desc,
            _ms(sql_ms), _ms(hql_ms), _ms(dl_ms), _ms(cy_ms), _ms(sp_ms),
            len(sql_rows),
        ))
        concise_rows.append((
            qid,
            count_lines(sql),
            count_tokens_sql(sql),
            count_tokens_hierql(hierql) if hierql else "—",
            count_tokens_datalog(datalog_forms[qid]) if qid in datalog_forms else "—",
            count_tokens_cypher(cypher_forms[qid]) if qid in cypher_forms else "—",
            count_tokens_sparql(sparql_forms[qid]) if qid in sparql_forms else "—",
        ))

        print(f"  SQL {len(sql_rows)} rows  {_ms(sql_ms)} ms"
              f"  HierQL {_ms(hql_ms)} ms"
              f"  Datalog {_ms(dl_ms)} ms"
              f"  Cypher {_ms(cy_ms)} ms"
              f"  SPARQL {_ms(sp_ms)} ms")
        print_rows_preview(sql_rows, args.preview_rows)
        print()

    w = csv.writer(sys.stdout)

    print("\n=== CORRECTNESS TABLE ===")
    w.writerow(["ID", "HierQL", "Datalog", "Cypher", "SPARQL"])
    w.writerows(correctness_rows)

    print(f"\n=== TIMING TABLE ({args.repeats} runs each) ===")
    w.writerow(["ID", "Description", "SQL ms", "HierQL ms", "Datalog ms", "Cypher ms", "SPARQL ms", "Row count"])
    w.writerows(timing_rows)

    print("\n=== CONCISENESS TABLE ===")
    w.writerow(["ID", "SQL lines", "SQL tokens", "HierQL tokens", "Datalog tokens", "Cypher tokens", "SPARQL tokens"])
    w.writerows(concise_rows)

    cur.close()
    conn.close()

    if souffle_tmpdir:
        import shutil
        shutil.rmtree(souffle_tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
