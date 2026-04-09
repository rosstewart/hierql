"""
HierQL preliminary experiments.
Prerequisites: pip install psycopg2-binary
Set DB_URL to your PostgreSQL 16 connection string.
Running this script populates the database, runs all benchmark queries
in hand-written SQL form, times them, counts tokens for multiple query
languages, and prints CSV blocks you can paste directly into the paper.
"""

import time
import csv
import sys
from hierql_transpiler import transpile

DB_URL = "postgresql://postgres:password@localhost:5434/bioql"

# ── Schema + sample data ──────────────────────────────────────────────────────

SETUP_SQL = """
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

# ── Benchmark queries (hand-written SQL) ─────────────────────────────────────
# Each entry: (id, description, sql)

QUERIES = [
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

# ── HierQL surface forms (for token counting / transpilation) ────────────────

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
WHERE gd.disease_id IN DISEASES UNDER \"Metabolic Disease\"
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

# ── Helpers ───────────────────────────────────────────────────────────────────

def count_tokens(sql: str) -> int:
    """Count non-whitespace tokens (words) in a query string."""
    return len(sql.split())

def count_lines(sql: str) -> int:
    return sum(1 for l in sql.strip().splitlines() if l.strip())

def run_query(cur, sql: str, repeats: int = 3) -> tuple[list, float]:
    """Run query, return (rows, mean_ms)."""
    rows = None
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        cur.execute(sql)
        rows = cur.fetchall()
        times.append((time.perf_counter() - t0) * 1000)
    return rows, sum(times) / len(times)

# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    try:
        import psycopg2
    except ImportError as e:
        raise SystemExit("Install psycopg2-binary to run the benchmark script.") from e

    conn = psycopg2.connect(DB_URL)
    conn.autocommit = True
    cur = conn.cursor()

    print("Setting up schema and data...", flush=True)
    cur.execute(SETUP_SQL)
    print("Done.\n")

    timing_rows  = []   # (qid, desc, mean_ms, n_rows)
    concise_rows = []   # (qid, sql_lines, sql_tokens, hierql_tokens, ...)

    for qid, desc, sql in QUERIES:
        rows, mean_ms = run_query(cur, sql)               # hand-written
        hierql = HIERQL_FORMS.get(qid, "")
        if hierql:
            t_sql = transpile(hierql)
            _, t_ms = run_query(cur, t_sql)               # transpiled
        else:
            t_ms = float('nan')
        timing_rows.append((qid, desc, f"{mean_ms:.2f}", f"{t_ms:.2f}", len(rows)))

        concise_rows.append((
            qid,
            count_lines(sql),
            count_tokens(sql),
            count_tokens(hierql) if hierql else "—",
            count_tokens(DATALOG_FORMS[qid]) if qid in DATALOG_FORMS else "—",
            count_tokens(CYPHER_FORMS[qid]) if qid in CYPHER_FORMS else "—",
            count_tokens(SPARQL_FORMS[qid]) if qid in SPARQL_FORMS else "—",
        ))

        print(f"{qid}: {len(rows)} row(s), {mean_ms:.2f} ms avg")
        for r in rows:
            print("  ", r)
        print()

    # ── Print CSV tables for the paper ───────────────────────────────────────

    print("\n=== TIMING TABLE (paste into paper) ===")
    w = csv.writer(sys.stdout)
    w.writerow(["ID", "Description", "Raw mean ms (3 runs)", "Transpiled mean ms (3 runs)", "Row count"])
    w.writerows(timing_rows)

    print("\n=== CONCISENESS TABLE (paste into paper) ===")
    w.writerow(["ID", "SQL lines", "SQL tokens", "HierQL tokens", "Datalog tokens", "Cypher tokens", "SPARQL tokens"])
    w.writerows(concise_rows)

    cur.close()
    conn.close()

if __name__ == "__main__":
    main()
