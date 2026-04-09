"""
hierql_transpiler.py  –  HierQL → PostgreSQL 16 SQL
Operators: UNDER, ANCESTORS, DEPTH, SIBLINGS (SIBLINGS is non-recursive: two-step lookup)
"""
import re
from dataclasses import dataclass
from typing import Optional

# ── Tokenizer ─────────────────────────────────────────────────────────────────

KEYWORDS = {
    'FIND', 'FROM', 'JOIN', 'ON', 'WHERE', 'GROUP', 'HAVING', 'ORDER', 'BY',
    'AS', 'AND', 'OR', 'NOT', 'IN', 'DISTINCT',
    'DISEASES', 'UNDER', 'ANCESTORS', 'DEPTH', 'SIBLINGS', 'REL',
}

@dataclass
class Tok:
    kind: str   # KW | ID | STR | NUM | OP
    val:  str
    def __repr__(self): return f"Tok({self.kind}, {self.val!r})"

_LEX = re.compile(
    r'"([^"]*)"'                                                # g1: "..."
    r"|'([^']*)'"                                               # g2: '...'
    r'|(\d+\.\d+|\d+)'                                         # g3: number
    r'|([<>!]=?|[=(),.*])'                                     # g4: op / punct
    r'|([A-Za-z_][A-Za-z0-9_]*(?:\.[A-Za-z_][A-Za-z0-9_]*)*)' # g5: word / dotted id
    r'|(\s+)',                                                   # g6: whitespace
)

def tokenize(src: str) -> list[Tok]:
    toks = []
    for m in _LEX.finditer(src):
        g1, g2, g3, g4, g5, g6 = m.groups()
        if   g6:             continue
        elif g1 is not None: toks.append(Tok('STR', g1))
        elif g2 is not None: toks.append(Tok('STR', g2))
        elif g3:             toks.append(Tok('NUM', g3))
        elif g4:             toks.append(Tok('OP',  g4))
        elif g5:             toks.append(Tok('KW' if g5.upper() in KEYWORDS else 'ID', g5))
    return toks

# ── AST ───────────────────────────────────────────────────────────────────────

@dataclass
class TableRef:
    table: str
    alias: str

@dataclass
class JoinClause:
    ref:  TableRef
    cond: str

@dataclass
class HierExpr:
    col:   str
    op:    str     # UNDER | ANCESTORS | DEPTH | SIBLINGS
    term:  str
    depth: Optional[int] = None
    rel:   Optional[str] = None

@dataclass
class Query:
    distinct:  bool
    select:    str
    from_ref:  TableRef
    joins:     list   # list[JoinClause]
    where:     list   # list[str | HierExpr]
    group_by:  str
    having:    str
    order_by:  str

# ── Parser ────────────────────────────────────────────────────────────────────

class Parser:
    def __init__(self, toks: list[Tok]):
        self.t   = toks
        self.pos = 0

    def peek(self, n: int = 0) -> Optional[Tok]:
        i = self.pos + n
        return self.t[i] if 0 <= i < len(self.t) else None

    def is_kw(self, *kws, n: int = 0) -> bool:
        p = self.peek(n)
        return p is not None and p.kind == 'KW' and p.val in kws

    def consume(self) -> Tok:
        tok = self.t[self.pos]; self.pos += 1; return tok

    def expect_kw(self, kw: str) -> Tok:
        if not self.is_kw(kw):
            raise SyntaxError(f"Expected keyword {kw!r}, got {self.peek()!r} at pos {self.pos}")
        return self.consume()

    def collect_raw(self, stop: set[str]) -> str:
        """Collect tokens into a SQL fragment, stopping before any keyword in `stop`."""
        parts = []
        while self.peek() and not (self.peek().kind == 'KW' and self.peek().val in stop):
            tok = self.consume()
            parts.append(f"'{tok.val}'" if tok.kind == 'STR' else tok.val)
        return ' '.join(parts)

    def parse_table_ref(self) -> TableRef:
        table = self.consume().val
        self.expect_kw('AS')
        return TableRef(table, self.consume().val)

    def parse_condition(self) -> 'str | HierExpr':
        """
        Detects:
          <col> IN DISEASES (UNDER|ANCESTORS|DEPTH|SIBLINGS) "term" [k] [,] [REL = "edge"]
        Falls back to raw string otherwise.
        """
        if self.is_kw('IN', n=1) and self.is_kw('DISEASES', n=2):
            col  = self.consume().val
            self.consume()              # IN
            self.consume()              # DISEASES
            op   = self.consume().val   # UNDER / ANCESTORS / DEPTH / SIBLINGS
            if op not in ('UNDER', 'ANCESTORS', 'DEPTH', 'SIBLINGS'):
                raise SyntaxError(f"Unknown hierarchy operator: {op!r}")
            term  = self.consume().val  # STR token (quotes already stripped by tokenizer)
            depth = int(self.consume().val) if op == 'DEPTH' and self.peek() and self.peek().kind == 'NUM' else None

            rel = None
            while self.peek():
                if self.peek().kind == 'OP' and self.peek().val == ',':
                    self.consume()
                    continue
                if self.is_kw('REL') or (self.peek().kind == 'ID' and self.peek().val.lower() == 'rel'):
                    self.consume()
                    if not (self.peek() and self.peek().kind == 'OP' and self.peek().val == '='):
                        raise SyntaxError("Expected '=' after rel")
                    self.consume()
                    if not self.peek() or self.peek().kind not in ('STR', 'ID'):
                        raise SyntaxError("Expected relation name after rel=")
                    rel = self.consume().val
                    continue
                if op == 'DEPTH' and depth is None and self.peek().kind == 'NUM':
                    depth = int(self.consume().val)
                    continue
                break

            return HierExpr(col, op, term, depth, rel)
        return self.collect_raw(stop={'AND', 'OR', 'GROUP', 'HAVING', 'ORDER'})

    def parse_cond_list(self) -> list:
        conds = [self.parse_condition()]
        while self.is_kw('AND'):
            self.consume()
            conds.append(self.parse_condition())
        return conds

    def parse(self) -> Query:
        self.expect_kw('FIND')

        distinct = False
        if self.is_kw('DISTINCT'):
            self.consume(); distinct = True

        select = self.collect_raw(stop={'FROM'})
        self.expect_kw('FROM')
        from_ref = self.parse_table_ref()

        joins = []
        while self.is_kw('JOIN'):
            self.consume()
            ref  = self.parse_table_ref()
            self.expect_kw('ON')
            cond = self.collect_raw(stop={'JOIN', 'WHERE', 'GROUP', 'HAVING', 'ORDER'})
            joins.append(JoinClause(ref, cond))

        where = []
        if self.is_kw('WHERE'):
            self.consume()
            where = self.parse_cond_list()

        group_by = ''
        if self.is_kw('GROUP'):
            self.consume(); self.expect_kw('BY')
            group_by = self.collect_raw(stop={'HAVING', 'ORDER'})

        having = ''
        if self.is_kw('HAVING'):
            self.consume()
            having = self.collect_raw(stop={'ORDER'})

        order_by = ''
        if self.is_kw('ORDER'):
            self.consume(); self.expect_kw('BY')
            order_by = self.collect_raw(stop=set())

        return Query(distinct, select, from_ref, joins, where, group_by, having, order_by)

# ── SQL emitter ───────────────────────────────────────────────────────────────

def _slug(s: str) -> str:
    return re.sub(r'[^a-z0-9]+', '_', s.lower()).strip('_')

def _cte_name(op: str, term: str, depth: Optional[int] = None, rel: Optional[str] = None) -> str:
    suffix = ""
    if depth is not None:
        suffix += f"_{depth}"
    if rel is not None:
        suffix += f"_{_slug(rel)}"
    return f"_{op.lower()}_{_slug(term)}{suffix}"

def _sql_str(s: str) -> str:
    return s.replace("'", "''")

def _uses_edge_table(h: HierExpr) -> bool:
    return h.rel is not None

def _make_cte(h: HierExpr) -> tuple[str, str]:
    name = _cte_name(h.op, h.term, h.depth, h.rel)
    q    = _sql_str(h.term)
    rel  = _sql_str(h.rel) if h.rel is not None else None

    if h.op == 'UNDER':
        if _uses_edge_table(h):
            # General edge table path expansion may revisit DAG nodes via multiple parents.
            body = (f"  SELECT id FROM Disease WHERE name = '{q}'\n"
                    f"  UNION\n"
                    f"  SELECT e.child_id\n"
                    f"  FROM   DiseaseEdge e\n"
                    f"  JOIN   {name} s ON e.parent_id = s.id\n"
                    f"  WHERE  e.rel = '{rel}'")
        else:
            # The sample adjacency-list schema is a tree, so UNION ALL is sufficient here.
            body = (f"  SELECT id FROM Disease WHERE name = '{q}'\n"
                    f"  UNION ALL\n"
                    f"  SELECT d.id FROM Disease d JOIN {name} s ON d.parent_id = s.id")
        cols = "id"

    elif h.op == 'ANCESTORS':
        if _uses_edge_table(h):
            body = (f"  SELECT e.parent_id AS id\n"
                    f"  FROM   Disease d\n"
                    f"  JOIN   DiseaseEdge e ON e.child_id = d.id\n"
                    f"  WHERE  d.name = '{q}' AND e.rel = '{rel}'\n"
                    f"  UNION ALL\n"
                    f"  SELECT e.parent_id\n"
                    f"  FROM   DiseaseEdge e\n"
                    f"  JOIN   {name} a ON e.child_id = a.id\n"
                    f"  WHERE  e.rel = '{rel}'")
        else:
            body = (f"  SELECT parent_id AS id FROM Disease WHERE name = '{q}' AND parent_id IS NOT NULL\n"
                    f"  UNION ALL\n"
                    f"  SELECT d.parent_id\n"
                    f"  FROM   Disease d JOIN {name} a ON d.id = a.id\n"
                    f"  WHERE  d.parent_id IS NOT NULL")
        cols = "id"

    elif h.op == 'DEPTH':
        if h.depth is None:
            raise ValueError("DEPTH requires an integer depth bound")
        if _uses_edge_table(h):
            body = (f"  SELECT id, 0 AS depth FROM Disease WHERE name = '{q}'\n"
                    f"  UNION ALL\n"
                    f"  SELECT e.child_id, s.depth + 1\n"
                    f"  FROM   DiseaseEdge e\n"
                    f"  JOIN   {name} s ON e.parent_id = s.id\n"
                    f"  WHERE  e.rel = '{rel}' AND s.depth < {h.depth}")
        else:
            body = (f"  SELECT id, 0 AS depth FROM Disease WHERE name = '{q}'\n"
                    f"  UNION ALL\n"
                    f"  SELECT d.id, s.depth + 1\n"
                    f"  FROM   Disease d JOIN {name} s ON d.parent_id = s.id\n"
                    f"  WHERE  s.depth < {h.depth}")
        cols = "id, depth"

    elif h.op == 'SIBLINGS':
        if _uses_edge_table(h):
            body = (f"  SELECT DISTINCT e2.child_id AS id\n"
                    f"  FROM   Disease d1\n"
                    f"  JOIN   DiseaseEdge e1 ON e1.child_id = d1.id\n"
                    f"  JOIN   DiseaseEdge e2 ON e2.parent_id = e1.parent_id\n"
                    f"  WHERE  d1.name = '{q}'\n"
                    f"    AND  e1.rel = '{rel}'\n"
                    f"    AND  e2.rel = '{rel}'\n"
                    f"    AND  e2.child_id != d1.id")
        else:
            # Non-recursive: parent lookup then sibling expansion in one query
            body = (f"  SELECT d2.id\n"
                    f"  FROM   Disease d1\n"
                    f"  JOIN   Disease d2 ON d2.parent_id = d1.parent_id\n"
                    f"  WHERE  d1.name = '{q}' AND d2.id != d1.id")
        cols = "id"

    else:
        raise ValueError(f"Unknown op: {h.op!r}")

    return name, f"{name}({cols}) AS (\n{body}\n)"

def emit_sql(q: Query) -> str:
    hiers = [c for c in q.where if isinstance(c, HierExpr)]
    seen, cte_frags = set(), []
    for h in hiers:
        name, frag = _make_cte(h)
        if name not in seen:
            seen.add(name); cte_frags.append(frag)

    where_parts = []
    for c in q.where:
        if isinstance(c, HierExpr):
            where_parts.append(
                f"{c.col} IN (SELECT id FROM {_cte_name(c.op, c.term, c.depth, c.rel)})"
            )
        else:
            where_parts.append(c)

    lines = []
    if cte_frags:
        lines.append("WITH RECURSIVE\n" + ",\n".join(cte_frags))

    lines.append(f"SELECT {'DISTINCT ' if q.distinct else ''}{q.select}")
    lines.append(f"FROM   {q.from_ref.table} AS {q.from_ref.alias}")
    for j in q.joins:
        lines.append(f"  JOIN {j.ref.table} AS {j.ref.alias} ON {j.cond}")
    if where_parts:
        lines.append("WHERE  " + "\n  AND  ".join(where_parts))
    if q.group_by:
        lines.append(f"GROUP BY {q.group_by}")
    if q.having:
        lines.append(f"HAVING {q.having}")
    if q.order_by:
        lines.append(f"ORDER BY {q.order_by}")

    return "\n".join(lines)


def transpile(hierql: str) -> str:
    """Main entry point: HierQL string → PostgreSQL 16 SQL string."""
    return emit_sql(Parser(tokenize(hierql)).parse())


# ── Self-test (run as: python hierql_transpiler.py) ───────────────────────────

if __name__ == '__main__':
    # Q2 needs GROUP BY in its HierQL form — update HIERQL_FORMS accordingly:
    #   "Q2": """FIND gene.symbol, COUNT(*) AS n_variants
    # FROM Gene AS gene
    #   JOIN Variant AS variant ON variant.gene_id = gene.id
    # GROUP BY gene.symbol
    # HAVING n_variants > 0
    # ORDER BY n_variants DESC"""

    import sys
    sys.path.insert(0, '.')
    try:
        from experiments import HIERQL_FORMS
    except ImportError:
        print("Place this file alongside experiments.py to self-test.")
        sys.exit(1)

    smoke_tests = dict(HIERQL_FORMS)
    smoke_tests["REL_DEMO"] = (
        'FIND d.name FROM Disease AS d '
        'WHERE d.id IN DISEASES UNDER "Metabolic Disease" REL = "part_of"'
    )

    ok = fail = 0
    for qid, src in sorted(smoke_tests.items()):
        print(f"\n{'─'*60}\n{qid}\n{'─'*60}")
        print(f"HierQL:\n{src}\n")
        try:
            sql = transpile(src)
            print(f"SQL:\n{sql}")
            ok += 1
        except Exception as e:
            print(f"ERROR: {e}")
            fail += 1
    print(f"\n{'='*60}\n{ok} ok, {fail} failed")
