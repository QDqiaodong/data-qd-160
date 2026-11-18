# Exact reaction batch material accounting v1
教学数据整理者输入reactants/products为各1..20个全局唯一化学式，available_moles恰含所有反应物与非负十进制字符串（整数最多9位、小数最多9位）。摩尔量用Fraction精确计算，不转换质量。
化学式支持[A-Z][a-z]?元素token、正整数下标<=1000和圆括号嵌套<=10层，长度<=100、每元素原子数<=10^12；不验证周期表符号真实性，不接受同位素、离子电荷、水合点或反应箭头。
构建原子守恒有理矩阵求零空间；必须恰一维并有全部严格正系数，约为最小互素正整数。欠定/不可正配平整批拒绝，不任意选解；不平衡电荷或判断实际反应可行性。
以available/coefficient最小值确定反应进度，可有并列限量试剂和零进度；输出消耗、剩余、产物理论摩尔上限及每进度原子守恒证据。样例2H2+O2→2H2O，5mol H2/2mol O2时进度2、H2余1、O2余0、H2O理论4mol。不是实验操作指南，不包含温压、动力学、安全性或真实产率。


## Lot mode (`reaction-lots`)

自本版本起提供独立的 `reaction-lots` 入口，供有序批次投料与质量隔离核算；
原有 `reaction-batch` 的输入、输出、退出码与报告结构保持不变。样例为
`fixtures/lots-example.json` 与 `fixtures/lots-diamond.json`。

`reaction-lots --source INPUT --output REPORT` adds ordered batch-charge and
quality-isolation accounting beside the unchanged `reaction-batch` command.
The two commands read different top-level contracts:

- `reaction-batch`: `reactants`, `products`, `available_moles` (v1, unchanged).
- `reaction-lots`: exactly `{"operations": [...]}`.

Both still use `common.read`/`common.save`: at most 8 MB of input, duplicate
JSON keys, non-finite numbers and undefined fields are rejected; the source
SHA-256 is reported; output is atomically replaced and symlinks are refused.
Invalid input or calculation exits 2 and leaves an existing report untouched.

## Operations

Every operation is an object with `op_id` (unique per file, `common.ident`)
and `type` exactly one of `stock`, `run`, `hold`, `release`. No other fields
are allowed.

| type | fields |
| --- | --- |
| `stock` | `op_id`, `lot`, `formula`, `moles` |
| `run` | `op_id`, `inputs`, `products` |
| `hold` | `op_id`, `lot` |
| `release` | `op_id`, `lot` |

- `lot` and `op_id` use `common.ident`.
- `moles` is a `quantity` decimal string (0..9 integer digits, 1..9 fraction
  digits) and must be greater than zero.
- `formula` follows the existing bounded formula grammar.
- `inputs` has 1..40 `{lot, moles}` objects with distinct lots.
- `products` has 1..20 `{lot, formula}` objects with distinct lots.
- Arrays must be JSON arrays; objects must have exactly the listed fields.
  JSON booleans are not accepted where a number or string is required.
- Limits: at most 30,000 operations, 30,000 successfully created lots,
  60,000 direct parent/child edges and 500 hold/release records (rejected
  ones count). Exceeding a limit or any structural/literal error aborts the
  whole batch with exit 2.

Global input errors include duplicate `op_id`, a duplicate input lot, a
duplicate product lot, an array-length breach, an unknown field or type, an
invalid identifier, formula or decimal. Balance failure and a non-common
extraction progress are business rejections, not input errors.

## Ledger semantics

- `stock` creates a new lot with exact `Fraction` inventory. A lot identity is
  never reusable, even after its inventory reaches zero or an operation reuses
  it. Reusing any existing lot identity (stock or product) is rejected.
- `run` resolves every input lot, merges input amounts by chemical formula,
  and uses the existing one-dimensional strictly positive integer balancing
  rule. All merged reactants and all products must number at most 20 species.
  Products must have distinct formulas and share no formula with a reactant.
  Every merged reactant amount must equal the same positive progress times its
  coefficient. Same-formula lots are distinct stock objects but their
  contributions are summed for the reaction.
- On success a `run` deducts each input amount, creates one new lot per
  product with inventory `progress * product coefficient`, and records every
  input lot as a direct parent of every product lot. It never re-solves for a
  greedy maximum extent.
- A failed `run` deducts nothing, creates no product, adds no lineage edge and
  does not reserve proposed product-lot identities. Processing continues.
- `hold` on an existing lot makes it an active isolation root. The isolated
  set is every descendant of every active root (on the accepted-run lineage)
  plus the roots themselves. Holding a lot that is already blocked through
  another root is allowed and makes it an additional independent root. Holding
  an already-active root, or holding an unknown lot, is rejected.
- `release` removes only that lot's own active-root status. Descendants release
  only when no remaining root covers them by any path. Releasing a lot that is
  not an active root, or an unknown lot, is rejected.
- A `run` whose input lot is in the current isolated set is rejected with no
  material or lineage change. Blocked state is derived from lineage, not from
  remaining inventory, so holding a depleted ancestor blocks existing products.
- All checks for one operation are gathered before any commit; no partial
  parent deduction or partial product creation is possible.

## Report

A successful report is exactly:

```json
{
  "schema_version": 1,
  "source_sha256": "<sha256 of source bytes>",
  "result": {
    "operations": [{"op_id": "...", "status": "accepted|rejected"}],
    "lots": [{"lot": "...", "formula": "...", "remaining_moles": "<Fraction>", "blocked": true}],
    "active_holds": ["..."]
  }
}
```

`operations` follows input order. `lots` and `active_holds` are sorted by lot
identifier. `remaining_moles` is the reduced `Fraction` string form.
`blocked` is a JSON boolean. Only final state is emitted: no ancestor lists,
per-operation snapshots or intermediate lot history. Identical source bytes
produce the same report.

Empty `operations` is valid and yields empty arrays.

## Example

`fixtures/lots-example.json` stocks two `H2` lots (0.1 and 0.2) and 0.15
`O2`, runs them into an `H2O` lot, then holds the first `H2` ancestor. The
water lot has 3/10 mol and is blocked even though every parent is depleted.

```sh
reaction-lots --source fixtures/lots-example.json --output run/lots-example.json
```

`fixtures/lots-diamond.json` builds a diamond with two active roots, then
releases one; the shared descendant remains blocked through the other root.

## Scale characteristics

Lineage is stored once as child adjacency; blocked state is recomputed
iteratively from active roots by traversal. No all-pairs reachability matrix
and no per-descendant ancestor-set copies are built, and no traversal uses
recursion. The acceptance-scale envelope (formulas <= 12 chars, <= 2 parent
lots and <= 2 products per run, 18-digit numerator/denominator intermediates)
fits in 192 MiB peak RSS on a 20,000-lot single chain and a convergence graph
within 60,000 edges. Other valid inputs are still processed exactly; the RSS
target is not claimed for them, and values are never truncated or approximated.

## Assay mode (`reaction-assay`)

`reaction-assay --source INPUT --output REPORT` adds a third, independent
command for teaching-data measurement consistency.  It asks whether one sample's
element molar measurements can be explained by a single nonnegative mixture of
the candidate species, and reports the exact tight range of each species amount.
The original `reaction-batch` and `reaction-lots` inputs, outputs and exit codes
are unchanged.  This is an algebraic consistency check of the supplied numbers
only; it makes no claim about whether a real reaction would proceed.

### Input

The input is exactly

```json
{
  "species": [
    {"species_id": "h", "formula": "H2"},
    {"species_id": "w", "formula": "H2O"}
  ],
  "samples": [
    {
      "sample_id": "s1",
      "total_moles": {"lower": "2", "upper": "2"},
      "measurements": [
        {"measurement_id": "mH", "element": "H", "lower": "3", "upper": "3"}
      ]
    }
  ]
}
```

- `species` has 1..60 entries.  `species_id` uses `common.ident` and is unique;
  `formula` uses the existing `chemistry.formula` grammar in full (brackets,
  multipliers, every existing limit) and its *text* is unique.  Two different
  formula texts with identical element counts (for example `H2O` and `OH2`)
  are kept as separate species and are never merged.
- `samples` has 1..8 entries with unique `sample_id`.  Samples are independent
  and share no variables.
- `total_moles` is `{"lower", "upper"}` and `measurements` has 0..120 entries.
  Every `measurement` is exactly `{"measurement_id", "element", "lower",
  "upper"}`, with `measurement_id` unique inside its sample and `element`
  matching ASCII `[A-Z][a-z]?` exactly.  The same element may appear several
  times, and an element may appear in no candidate formula.
- All quantities are `common.quantity` nonnegative decimal strings (JSON
  numbers and booleans are rejected); every interval must have
  `lower <= upper`.  Equal bounds mean an exact equality.
- Arrays and nested objects accept only the fields listed above.  Inputs that
  exceed a limit, repeat an identifier, use an unknown field, an invalid
  formula or element, an inverted interval or a non-string quantity are
  whole-batch input errors: exit 2, no partial report, and any existing report
  is left untouched.  A legal but contradictory sample is *not* an input error;
  it is reported as `inconsistent` in a successful report.

`common.read`/`common.save` still apply: at most 8 MB, duplicate JSON keys and
non-finite numbers rejected, `source_sha256` bound to the exact source bytes,
refusal of an output that is a symlink or the same path as the input, and an
atomic successful replacement.  No database, external solver, network call or
new runtime dependency is used.

### Constraints

For one sample let `n` be the number of species, `x_j >= 0` the moles of
species `j`, and `a_ej` the number of atoms of element `e` in species `j`.  The
constraints are

```
sum_j x_j                 within total_moles         (both bounds)
sum_j a_ej x_j            within every measurement    (all of them)
```

Measurement values are element *atom* moles, not species moles.  There are
`n + 2 + 2m` canonical rows with `m` the number of measurements, in this fixed
order:

1. `n` rows `-x_j <= 0`;
2. `sum_j x_j <= total_moles.upper`;
3. `-sum_j x_j <= -total_moles.lower`;
4. for each measurement in input order, `a_e x <= upper` then `-a_e x <= -lower`.

Zero, repeated and linearly dependent rows keep their positions.  Because
`total_moles.upper` is finite, every feasible sample has a finite tight bound
for every species.  Midpoints are never used and no per-element answer is
assembled separately.

### Report

A successful report is exactly

```json
{
  "schema_version": 1,
  "source_sha256": "<sha256 of source bytes>",
  "result": {
    "species_order": ["..."],
    "samples": [
      {"sample_id": "...", "status": "consistent", "ranges": [
        {"species_id": "...",
         "minimum": {"value": "...", "composition": ["..."], "dual": ["..."]},
         "maximum": {"value": "...", "composition": ["..."], "dual": ["..."]}}
      ]},
      {"sample_id": "...", "status": "inconsistent", "contradiction": ["..."]}
    ]
  }
}
```

`species_order` follows the input species order; `samples` follows input order;
`ranges` follows input species order.  A consistent sample always has one
`ranges` entry per species.  An inconsistent sample has no `ranges` key.

All numbers are reduced `Fraction` strings (`"0"`, `"2"`, `"1/3"`), never JSON
numbers, never negative zero.

### Certificates

Write the canonical system as `A x <= b` with `R = n + 2 + 2m` rows.  For a
**maximum** of species `j`, `c = e_j`; for a **minimum**, `c = -e_j`.  Each
block carries

- `composition`, a feasible point of length `n`: `composition >= 0`,
  `A composition <= b`, and `composition[j] == value`;
- `dual`, a nonnegative vector of length `R` with `A^T dual == c` and
  `b . dual == value` for a maximum, or `A^T dual == -e_j` and
  `b . dual == -value` for a minimum.

The value together with these two vectors proves both attainability and
optimality: any feasible `x` satisfies `x_j <= b . dual` (or `x_j >= -b . dual`).
Different bounds may use different witnesses; the report is a pure function of
the source bytes, so re-running the same source produces identical output.

An inconsistent sample carries `contradiction`, a nonnegative vector of length
`R` with `A^T contradiction == 0` and `b . contradiction == -1`.  This is an
exact Farkas refutation; no tolerance is applied, so a difference of one
nanomole is detected.

For solving, parallel measurement rows of one element are tightened to one
lower/upper pair and the obtained multipliers are spread back onto the exact
canonical rows named by the input, so every reported certificate is checkable
against `A` and `b` exactly as defined above.  The simplex works on the dual
with only `n` rows, uses exact `Fraction`/integer arithmetic throughout, and
selects entering/leaving variables with Dantzig pricing and a Bland fallback so
it cannot cycle.

### Examples

`fixtures/assay-example1.json`: `H2`, `O2`, `H2O` with total 2, H 3, O 1.5.
The unique composition is `["1", "1/2", "1/2"]` and every minimum equals its
maximum.

`fixtures/assay-example2.json`: `H2O` and `OH2` with total 1, H 2, O 1.  Both
species range over `["0", "1"]`; the minimum and maximum use different
witness compositions.

`fixtures/assay-example3.json`: only `H2O` with total 1, H 2 and O
`1.000000001`.  It is `inconsistent` with an exact contradiction vector.

`fixtures/assay-mixed-48.json`: 48 species with overlapping elements and
several degrees of freedom, plus a redundant/duplicate-measurement sample, a
jointly contradictory sample, an empty-measurements sample with
`total_moles = 0`, and a `total_moles = 0` sample that still requires hydrogen.

```sh
reaction-assay --source fixtures/assay-example1.json --output run/assay1.json
reaction-assay --source fixtures/assay-mixed-48.json --output run/assay48.json
```

## Peaks mode (`reaction-peaks`)

`reaction-peaks --source INPUT --output REPORT` adds a fourth, independent
command for shared-drift chromatographic peak identity reconciliation.  A batch
of observed peaks is checked against a reference catalogue under one common
instrument drift per run; the command reports the largest explainable peak
count, every drift that attains it, and the identity each peak takes across all
optimal explanations.  It never promotes one arbitrary matching into a unique
identity, and it says nothing about real chemical structure or reaction
feasibility.  The original `reaction-batch`, `reaction-lots` and
`reaction-assay` inputs, outputs and exit codes are unchanged.

### Input

The input is exactly

```json
{
  "references": [
    {"reference_id": "a", "formula": "A",
     "window": {"lower": "0", "upper": "1"}}
  ],
  "runs": [
    {"run_id": "r1", "drift": {"lower": "-1", "upper": "3"},
     "peaks": [{"peak_id": "p", "formula": "A", "time": "3"}]}
  ]
}
```

- `references` has 1..40 entries.  `reference_id` uses `common.ident` and is
  globally unique; `formula` uses the existing `chemistry.formula` grammar in
  full; `window` is exactly `{"lower", "upper"}` with `lower <= upper`.
- `runs` has 1..8 entries with unique `run_id`.  `drift` is exactly
  `{"lower", "upper"}` with `lower <= upper`; `peaks` has 0..40 entries with
  `peak_id` unique inside its run (the same `peak_id` may appear in other
  runs).
- Times and window bounds are `common.quantity` nonnegative decimal strings;
  drift bounds use the same syntax with one optional leading minus sign (the
  absolute value must still pass `quantity`; `-0` is accepted and normalised
  to `0` on output).  JSON numbers and booleans are rejected everywhere; no
  whitespace, plus signs, exponents or fraction notation.
- Formula texts may repeat freely across references and across peaks; two
  different texts with identical element counts (for example `H2O` and `OH2`)
  are merely compatible and never merge identities.
- Runs are analysed independently against the shared catalogue; a reference
  occupied in one run is still available in the others.
- Any unknown field, duplicate identifier, invalid formula or number, inverted
  interval or breached array limit — in any run — is a whole-batch input
  error: exit 2, no partial report, and any existing report is left untouched.
  It never degrades into `K = 0` or a partial publication of earlier runs.

`common.read`/`common.save` still apply: at most 8 MB, duplicate JSON keys and
non-finite numbers rejected, `source_sha256` bound to the exact source bytes,
refusal of an output that is a symlink or the same path as the input, and an
atomic successful replacement.  No database, external solver, network call or
new runtime dependency is used.

### Semantics

Inside one run every peak is explained by at most one reference and every
reference explains at most one peak.  A peak with time `t` may pair with a
reference whose window is `[lo, hi]` exactly when the formula element counts
are equal and `lo + d <= t <= hi + d`, where `d` is the single shared drift of
the run, drawn from the closed input drift domain.  The pair is therefore
feasible precisely on the closed drift interval `[t - hi, t - lo]` clipped to
the domain.

For a drift `d`, an *explanation* is any set of legal peak/reference pairs
with no shared endpoint; unmatched peaks and references are allowed.  `K` is
the largest explanation size over all drifts and all explanations — there is
no secondary objective such as distance minimisation.  Two peaks that each can
be matched at some drift need not admit one common drift matching both.

The **best drift set** is every `d` in the closed domain that admits a
`K`-sized explanation.  It is reported as merged closed intervals sorted by
`lower`; it may be disconnected and may contain single points, because a
closed endpoint can beat every open neighbourhood around it.  Overlapping or
touching closed intervals are merged; each reported interval has
`lower <= upper` and no open ends or floating-point tolerances are used.

For every peak the report lists the union of `reference_id` values it takes in
at least one `K`-sized explanation over the whole optimal drift set, plus
`can_be_unmatched`, which is true iff some `K`-sized explanation leaves the
peak completely unmatched.  Support can only come from a complete `K`-matching
at one legal drift: edges that merely satisfy the window in isolation, or that
exist only at drifts where no `K`-matching exists, are excluded.  An identity
is certain only when the support is a single reference and
`can_be_unmatched` is false.  `K = 0` is a normal result: the whole input
drift domain is optimal, every support is empty and every
`can_be_unmatched` is true.

The **witness** is a reproducible optimal explanation: its drift is the
minimum of the best drift set, and among the `K`-matchings at that single
drift it chooses the lexicographically smallest assignment vector — peaks
ordered by `peak_id` (Unicode code points), `reference_id` compared by Unicode
code points, null after every `reference_id`.  It is never assembled by
stitching supports across different drifts.  With empty `peaks`, `K` is 0, the
witness drift is the input lower bound and `assignments` is empty.

### Report

A successful report is exactly

```json
{
  "schema_version": 1,
  "source_sha256": "<sha256 of source bytes>",
  "result": {
    "runs": [
      {"run_id": "...", "maximum_matches": 2,
       "best_drifts": [{"lower": "1", "upper": "1"}],
       "peaks": [{"peak_id": "...", "reference_ids": ["..."],
                  "can_be_unmatched": false}],
       "witness": {"drift": "1",
                   "assignments": [{"peak_id": "...", "reference_id": "..."}]}}
    ]
  }
}
```

`runs` follows the input run order.  `peaks` and `assignments` are sorted by
`peak_id` (Unicode code points); `assignments` lists every peak, with JSON
null for an unmatched one.  `reference_ids` is sorted, duplicate-free.
`maximum_matches` is a JSON integer and `can_be_unmatched` a JSON boolean.
All drift bounds in the output are reduced `Fraction` strings (`"0"`, `"-1"`,
`"1/2"`), never JSON numbers, never non-reduced fractions, never negative
zero.  Identical source bytes produce identical reports.

### Algorithm and correctness organisation

All comparisons use exact `Fraction` arithmetic.  Per run, each compatible
peak/reference pair contributes one closed drift interval, so the feasibility
graph `E(d)` changes only at interval endpoints.  The drift domain is swept
left to right over *regions* — each endpoint point and each open segment
between consecutive endpoints — instead of solving every region independently
from scratch: a single maximum matching is maintained incrementally, adding
the edges that open at a point and dropping the ones that close after it,
repaired each time by breadth-first augmenting paths.  No matching or peak
permutation is ever enumerated.

Closed endpoints stay exact because `E(x)` at a point is a superset of the
adjacent segment edge sets: the matching size at a point dominates its
neighbours, so an optimal open segment always has optimal endpoints and the
optimal drift set is a union of closed regions whose closures merge into the
reported closed intervals; isolated optimal points survive as `[x, x]`.

Identity supports are computed per region from one maximum matching via the
standard alternating-path characterisation (orient unmatched edges
left-to-right and matched edges right-to-left): an edge is in some maximum
matching of the region iff it is matched, lies on a directed cycle (same
strongly connected component), is reachable from an unmatched peak, or reaches
an unmatched reference; a peak can be unmatched iff it is reachable from an
unmatched peak (itself included).  Supports and `can_be_unmatched` accumulate
only over regions whose matching size equals the current best, resetting
whenever a larger size is found, so the final values are exactly the union
over the whole optimal drift set.  The witness reuses the same
characterisation greedily at the minimum optimal drift: peaks are fixed in
`peak_id` order, each taking the smallest still-feasible `reference_id` (null
only when none remains), which yields the documented lexicographic minimum.

### Examples

`fixtures/peaks-example1.json`: references `a(A,[0,1])`, `b(B,[0,1])` and
peaks `p(A,1)`, `q(B,2)` with drift domain `[0,2]`.  Only the single drift
`1` aligns both pairs, so `K = 2` and `best_drifts` is exactly `[[1, 1]]`;
every open neighbourhood of `1` explains only one peak.

`fixtures/peaks-example2.json`: references `a(A,[0,1])`, `b(A,[3,4])` and one
peak `p(A,3)` with drift domain `[-1,3]`.  `K = 1`, the optimal drifts are
`[-1, 0]` and `[2, 3]`, `p` supports both `a` and `b` but cannot be unmatched,
and the canonical witness uses drift `-1` with `p` assigned to `b`.

`fixtures/peaks-linkage.json`: references `a(A,[0,0])`, `b(A,[1,1])`,
`c(B,[0,1])` and peaks `p(A,1)`, `q(B,1)` with drift domain `[0,1]`.  `K = 2`
with optimal drifts exactly the two single points `0` and `1`; `p` supports
both `a` and `b`.  Narrowing `c` back to `[0,0]` would leave only drift `1`
optimal and `p` supporting only `a` — a window change in one element group
reshapes the identity support of another.

`fixtures/peaks-scale-40.json`: 40 references and 40 peaks, all formula `A`,
windows `[0,1]`, times `1`, drift domain `[0,1]`.  `K = 40`, the whole domain
is optimal, every peak supports all 40 references and none can be unmatched;
the canonical witness pairs peaks and references in identifier order.  The
support sets are derived from the alternating structure, not by enumerating
the `40!` perfect matchings.

```sh
reaction-peaks --source fixtures/peaks-example1.json --output run/peaks1.json
reaction-peaks --source fixtures/peaks-scale-40.json --output run/peaks40.json
```

## Isotopes mode (`reaction-isotopes`)

`reaction-isotopes --source INPUT --output REPORT` adds a fifth, independent
command for the joint isotope distribution of a molecule and the exclusive
detection probability of user-defined channels.  It is a forward probability
calculation over teaching data: every isotope mass number and abundance, and
every `tagged` flag, comes from the input; no natural-abundance or real
periodic table is built in.  The original `reaction-batch`, `reaction-lots`,
`reaction-assay` and `reaction-peaks` inputs, outputs and exit codes are
unchanged.

### Input

The input is exactly

```json
{
  "elements": [
    {"element": "A",
     "isotopes": [
       {"mass_number": 10, "abundance": "0.5", "tagged": false},
       {"mass_number": 11, "abundance": "0.5", "tagged": true}]}
  ],
  "queries": [
    {"query_id": "q", "formula": "A2", "cutoff": 2,
     "channels": [
       {"channel_id": "c1", "mass": {"lower": 20, "upper": 22},
        "tag_count": {"lower": 0, "upper": 2}}]}
  ]
}
```

- `elements` has 1..20 entries.  `element` is an ASCII `[A-Z][a-z]?` token,
  unique across the array; it is *not* checked against the periodic table.
  `isotopes` has 1..6 entries with distinct `mass_number`.  `mass_number` is a
  JSON integer 1..300 (booleans are rejected).  `abundance` uses the existing
  `quantity` nonnegative decimal-string grammar and range, plus `<= 1`; the
  abundances of one element sum to exactly 1.  A zero-abundance entry is legal
  and still participates in the minimum-`mass_number` base definition.
  `tagged` is a JSON boolean.
- `queries` has 1..4 entries with unique `query_id` (`common.ident`).
  `formula` uses the existing `chemistry.formula` grammar and its limits, and
  additionally the parsed total atom count must be `<= 80`; every element that
  actually appears must have an isotope table (unused tables are legal).
  `cutoff` is a JSON integer 0..128 (booleans rejected).  `channels` has 0..16
  entries, in priority order, with `channel_id` unique inside the query.
- `mass` and `tag_count` are exactly `{"lower", "upper"}`.  Mass bounds are
  JSON integers 0..10^9 with `lower <= upper` and `upper <= B + cutoff`;
  tag bounds are JSON integers 0..80 with `lower <= upper` and may exceed the
  atom count, producing an empty gate.  Every integer rejects booleans.
- All arrays and nested objects accept only the fields listed above.  Unknown
  fields, duplicate identifiers, a missing element table, an over-limit atom
  count, a bad abundance sum, an inverted or over-limit interval, or any
  literal error in any query is a whole-batch input error: exit 2, no partial
  report, and any existing report is left untouched.

`common.read`/`common.save` still apply: at most 8 MB, duplicate JSON keys and
non-finite numbers rejected, `source_sha256` bound to the exact source bytes,
refusal of an output that is a symlink or the same path as the input, and an
atomic successful replacement.  No database, external solver, network call or
new runtime dependency is used.

### Semantics

For element `e` with `n_e` atoms, `b_e` is the minimum `mass_number` among
*all* of its input isotope entries (including zero-abundance ones).  The base
mass is `B = sum_e n_e * b_e`.  Each atom is sampled independently from its
element's table; the nominal mass `M` is the sum of the sampled mass numbers
and `T` counts the samples whose `tagged` flag is true.  The offset
`S = M - B` is always nonnegative.  Results with `0 <= S <= cutoff` are
retained; those with `S > cutoff` are `overflow`.

The probability of a `(M, T)` pair is the sum over all atom assignments: two
atoms taking one light and one heavy isotope are two distinct assignments, so
their weights add.  The joint law is never assembled from the two marginals.

Channels are checked in input order; each retained `(M, T)` is assigned to the
first channel whose closed `mass` interval contains `M` **and** whose closed
`tag_count` interval contains `T`.  An outcome assigned to no channel is
`missed`, and includes every truncated tail.  `overflow` is the `S > cutoff`
part, which is therefore a subset of `missed`; it is never added to `missed`
again.  Because assignment depends on the joint event, shrinking an earlier
channel changes the exclusive probability and conditional tag moments of later
channels even when their own bounds are unchanged.  Queries are independent:
one query never consumes probability or isotope tables from another.

### Statistics object

Every event block is exactly

```json
{"probability": "p", "first_moment": "u", "second_factorial_moment": "v",
 "mean_tags": "mu", "variance_tags": "var"}
```

with `p = P(E)`, `u = E[T * 1_E]`, `v = E[T*(T-1) * 1_E]`.  When `p > 0`,
`mu = u/p` and `var = (v + u)/p - mu^2`; when `p = 0`, `u = v = "0"` and
`mu` and `var` are JSON null.  All numbers are reduced `Fraction` strings
(`"0"`, `"1"`, `"1/3"`), never JSON decimals, never non-reduced fractions,
never negative zero.  A tiny negative variance is not rounded to zero by
tolerance; the exact arithmetic cannot produce one.

### Report

A successful report is exactly

```json
{
  "schema_version": 1,
  "source_sha256": "<sha256>",
  "result": {"queries": [
    {"query_id": "...", "base_mass": 20,
     "envelope": [{"mass": 20, "probability": "1/4"}],
     "overflow": {"...": "..."},
     "channels": [{"channel_id": "...", "stats": {"...": "..."}}],
     "missed": {"...": "..."},
     "total": {"...": "..."}}
  ]}
}
```

`queries` follows input order.  `base_mass` is the JSON integer `B`.
`envelope` lists every retained mass with positive probability in ascending
mass order; zero-probability masses are omitted (sparse gaps are not filled).
`channels` follows input order and always keeps zero-probability channels.
`total` is the untruncated molecule statistic with `probability` exactly `"1"`.

The envelope probabilities plus `overflow.probability` sum to 1.  The
`channels` and `missed` blocks partition every outcome, so their three moment
components (`p`, `u`, `v`) add term by term to `total`; conditional means and
variances are not additive.  Identical source bytes produce identical reports.

### Algorithm

The implementation keeps an explicit joint state `(S, T)` for retained offsets
and multiplies the per-atom choices cycle by cycle; a state whose offset has
already passed `cutoff` is dropped because offsets never decrease, and its
probability and tag moments remain exactly recoverable as `total - kept`.  The
alternative exact representation — the full joint state without the cutoff
pruning plus a monotone-tag tail moment summary — is equivalent: both preserve
the same atom-position multiplicity and the same first-come channel boundary,
and both conserve the truncated tail's probability and moments.  Small
molecules are checked in the tests by an independent per-position enumeration
that does not call this solver or any recurrence.

### Examples

`fixtures/isotopes-example1.json`: element `A` with masses 10 (untagged) and
11 (tagged), both 1/2.  `A2`, `cutoff = 2`, no channels: the joint law is
`(20,0):1/4`, `(21,1):1/2`, `(22,2):1/4`, so `total` has `u = 1`,
`v = 1/2`, `mean_tags = 1`, `variance_tags = 1/2`.  The `gate` query's channel
`mass=[21,21]`, `tag_count=[0,0]` has probability 0 even though each marginal
is positive.  The `linkage` query shows a prefix change reshaping the next
channel's conditional moments.

`fixtures/isotopes-example2.json`: elements `A`, `B` and `C` with several
isotopes, zero-abundance entries and `tagged` flags whose mass offset differs
from the tag increment; it exercises a nested formula, a sparse envelope, a
tagged overflow tail and empty gates.

```sh
reaction-isotopes --source fixtures/isotopes-example1.json --output run/isotopes1.json
reaction-isotopes --source fixtures/isotopes-example2.json --output run/isotopes2.json
```
