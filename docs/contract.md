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

## Split mode (`reaction-split`)

`reaction-split --source INPUT --output REPORT` adds a sixth, independent
command for teaching-data hygiene: it turns classified chemistry samples into
a train/test listing in which records that could leak into the other side stay
together.  It does not train a model, infer a chemical structure or claim that
two formulas with equal element counts denote the same real molecule.  The
original five commands keep their inputs, reports and exit behaviour.

### Input

The input is exactly

```json
{
  "records": [
    {"record_id": "a", "formula": "H2O", "experiment_id": "e1",
     "label": "positive", "pin": "free"}
  ],
  "test_quotas": {
    "positive": {"lower": 1, "upper": 1, "target": 1},
    "negative": {"lower": 0, "upper": 0, "target": 0}
  }
}
```

- `records` has 0..600 entries.  Each entry is exactly
  `{record_id, formula, experiment_id, label, pin}`.
  `record_id` uses `common.ident` and is globally unique.
  `formula` uses the full `chemistry.formula` grammar and every existing
  limit; the periodic table is not consulted and the isotope command's
  80-atom bound does not apply.  `experiment_id` uses `common.ident` and may
  repeat.  `label` is exactly `positive` or `negative`; `pin` is exactly
  `free`, `train` or `test`.
- At most 300 records may carry one label.
- `test_quotas` is exactly `{positive, negative}` and each block is exactly
  `{lower, upper, target}`, all JSON integers 0..300 with booleans rejected
  and `lower <= upper`.  `target` need not lie inside the interval or inside
  the actual label count.  `lower` above the actual total is legal but may be
  infeasible.  Quotas constrain the raw test record counts per label, not
  group counts; there is no train quota and no total-test quota.
- Arrays and nested objects accept only the listed fields.  Any unknown field,
  duplicate `record_id`, invalid formula, label or pin, over-limit array or
  per-label count, or literal error (including a boolean where an integer is
  required) is a whole-batch input error: exit 2, no partial report, and an
  existing report is left untouched.  A legal input that cannot be split is
  *not* an error; it produces an `infeasible` report.

`common.read`/`common.save` still apply: at most 8 MB, duplicate JSON keys and
non-finite numbers rejected, `source_sha256` bound to the exact source bytes,
refusal of an output that is a symlink or the same path as the input, and an
atomic successful replacement.  No database, external solver, network call or
new runtime dependency is used.

### Leakage groups

Two records are linked when their parsed formulas have exactly the same
element-key to integer-count mapping, or when they share an `experiment_id`.
Both relations are undirected, and the reported groups are their joint
transitive closure: a chain such as equal formula, then shared experiment,
then equal formula again lands in one group.  `H2O` and `OH2` are linked
because their element counts are identical; `H4O2` is not equal to `H2O`
(no empirical-formula reduction is performed).  Labels never cut a group, so
a group can hold both labels, and grouping is decided before quota planning.
Each group keeps every original record: nothing is dropped or merged into a
single sample.

### Planning

`pin = free` imposes nothing.  A group containing at least one `train` pin is
forced to train, one containing at least one `test` pin is forced to test,
and a group containing both is a conflict.  Compatible same-direction pins
are harmless.  All groups are checked for pin conflicts before any quota
arithmetic; if any conflict exists the result is `pin_conflict` even when
quotas would also fail.

Otherwise let `P` and `N` be the test-side positive and negative record
counts and `D = |P - positive.target| + |N - negative.target|`.  A complete
listing is valid when every group is on one side, every pin is honoured and
both test counts lie in their closed quota intervals.  Among all valid
listings the canonical one minimises `D`; ties are broken by the 0/1
assignment vector over all record ids in Unicode code-point order with
`train = 0` and `test = 1`, taking the lexicographically smallest.  Because a
group cannot be split, this is the same as comparing group assignments in
increasing minimum-`record_id` order, and the tie-break is resolved on the
original record identity vector rather than by choosing a count cell first.
Label order, source array order, disjoint-set roots and hash iteration order
never break a tie.  Exactly one canonical listing is reported.

### Report

A successful report is exactly

```json
{
  "schema_version": 1,
  "source_sha256": "<sha256>",
  "result": {"groups": [...], "status": "feasible",
             "split": {"train": ["..."], "test": ["..."],
                       "test_counts": {"positive": 1, "negative": 0},
                       "deviation": 0}}
}
```

`result` always has `groups` and `status`.  Each group is exactly
`{"group_id", "record_ids", "counts", "forced"}`: `group_id` is its smallest
`record_id`, `record_ids` is sorted by Unicode code point (so `group_ids` are
increasing down the array), `counts` adds the original records per label, and
`forced` is `free`, `train`, `test` or `conflict`.  Counts are record counts,
never formula or experiment counts.

Feasible results add `split`, whose `train` and `test` are sorted disjoint
`record_id` arrays covering every record, plus the integer `test_counts` and
`deviation`.  An infeasible result with a pin conflict is exactly
`{"groups": [...], "status": "infeasible", "reason": "pin_conflict",
"conflict_groups": [...]}` with every `forced = conflict` id sorted.  An
infeasible result without a pin conflict is exactly
`{"groups": [...], "status": "infeasible", "reason": "quotas"}`.  No
infeasible report carries `split`, an approximate listing or a fabricated
deviation; groups are never split, records never dropped and quotas never
relaxed.  All counts and `deviation` are JSON integers.

Empty `records` with both `lower` bounds 0 is feasible with empty arrays,
`test_counts` of 0 and `deviation` equal to the two targets; otherwise it is
`quotas`-infeasible.

### Algorithm and state

Grouping is a disjoint-set union over a first-seen index per element-count
tuple and per experiment id, so each record triggers at most two unions and
the closure is exact by construction.

Planning uses only free groups.  For every suffix of free groups the planner
stores a reachability table: one Python integer per possible positive count,
whose bits are the reachable negative counts.  The table is at most 301
integers of at most 301 bits whatever the number of groups, and every
snapshot step costs at most 301 shifted ORs; memory is bounded by
`O(free_groups * 301 * 301 bits)` of shared rows, far below the 192 MiB
ceiling at 600 records.  Targets and quotas give at most a few hundred
`(P, N)` pairs at the optimal deviation, found by scanning the reachable
table rather than assigning groups.

The canonical assignment is then rebuilt over all groups sorted by
`group_id`: forced groups are fixed, and each free group is placed on train
if the remaining suffix can still complete an optimal in-quota assignment and
on test otherwise.  Preferring train at the earliest code point is exactly
the documented lexicographic minimum, and every step consults the same exact
reachability table, so the global deviation, quota membership and
canonical identity are preserved together.  No subset of group assignments
is enumerated, and no random seed or greedy count choice stands in for the
optimum.

### Examples

`fixtures/split-example1.json`: `a(H2O,e1)`, `b(OH2,e2)`, `c(CO2,e2)`,
`d(O2,e3)`.  `a` and `b` share element counts, `b` and `c` share `e2`, so
`a, b, c` form one multi-hop group and `d` is alone.  Quotas
`positive [1,2] target 1`, `negative [0,2] target 1` yield
`test = [a, b, c]`, `deviation = 1`.

`fixtures/split-example2.json`: six free records in two formula groups
(`a, b, c` with `A`; `d, e, f` with `B`), both quotas `[1,2] target 1`.  The
two count choices `(1,2)` and `(2,1)` both have `D = 1`; the canonical test
side is `[d, e, f]`, not the smaller-positive count cell.

`fixtures/split-example3.json`: the bridge case.  As written, `a` and `b`
share `e1` while pinned `test` and `train`, so the result is
`pin_conflict` with `conflict_groups = [a]`.  Removing the shared experiment
(or changing `b.experiment_id` back to `e2`) separates the groups and the same
pins become feasible.

`fixtures/split-scale-600.json`: 600 independent free components
`r000..r599`, formulas `A1..A600`, unique experiments, 300 records per label;
both quotas are exact `150` with target `150`.  The canonical test side is
`r150..r299` plus `r450..r599`, `deviation = 0`, 600 groups; the CLI run stays
well under 192 MiB peak RSS.

```sh
reaction-split --source fixtures/split-example1.json --output run/split1.json
reaction-split --source fixtures/split-example2.json --output run/split2.json
reaction-split --source fixtures/split-example3.json --output run/split3.json
reaction-split --source fixtures/split-scale-600.json --output run/split600.json
```

## Certify mode (`reaction-certify`)

`reaction-certify --source INPUT --output REPORT` adds a seventh, independent
command for offline teaching concentration records: it maintains traceable
piecewise instrument calibrations, corrected observations and protected
certificates over an ordered batch stream, and reports the calibration values
and certificate evidence knowable at each requested batch prefix.  It is a
deterministic accounting of supplied teaching data: one file is computed from
the empty state, no state persists across files, and nothing is inferred about
real instrument performance or experimental safety.  The original six commands
keep their entries, inputs and reports.

### Input

The input is exactly `{"batches": [...], "queries": [...]}` with no other
top-level keys.  Repeated JSON keys, unknown keys, non-finite values, invalid
identifiers/formulas/quantities, out-of-bound values, duplicate batch or query
ids, a repeated `(entity, id)` inside one batch, illegal intervals or curve
nodes, duplicate/non-positive certificate members, inverted thresholds and an
illegal `as_of` are all whole-file input errors (exit 2, no partial report,
an existing report left untouched).  The complete file is structurally
validated before any batch runs, so a malformed trailing record cannot be
hidden by earlier valid batches.

- `batches` has 0..300 entries processed in array order.  Each is exactly
  `{"batch_id", "changes"}`; `batch_id` uses `common.ident`, is case
  sensitive and unique across the file.  A batch may carry an empty `changes`
  array.  All changes across the file total at most 6000.
- A change is exactly
  `{"entity", "action": "put", "record"}` or
  `{"entity", "action": "delete", "id"}`, where `entity` is exactly
  `calibration`, `observation` or `certificate`.  The three entities have
  independent id spaces.  At most one change per `(entity, id)` may appear in
  one batch; a deleted id may be recreated in a later batch.
- Calibration: exactly
  `{"id", "instrument", "start", "end", "points"}`.  `start`/`end` are JSON
  integers 0..10^9 with `start < end` (booleans rejected).  `points` has
  2..64 `{"raw", "value"}` nodes with strictly increasing `raw`; `value`
  need not be monotone.  All put calibration nodes across the file total at
  most 20000.  Every `raw`/`value` is a nonnegative `chemistry.quantity`
  decimal string (JSON numbers and fraction literals are rejected).
- Observation: exactly
  `{"id", "instrument", "time", "formula", "raw"}`.  `time` is an integer
  0..10^9, `formula` uses the full existing formula grammar and `raw` is a
  quantity string.  An observation may legally have no covering calibration
  or a raw value outside the curve domain; that only matters when a
  certificate asks it for support.  `instrument` is only a matching key.
- Certificate: exactly
  `{"id", "formula", "members", "lower", "upper", "max_range"}`.  `members`
  has 1..64 entries `{"observation_id", "weight"}` with distinct
  `observation_id`s and positive weights; `lower <= upper`.  Member entries
  across all put certificate records total at most 2000.  One observation may
  support several certificates.  `weight`, `lower`, `upper`, `max_range` are
  quantity strings.
- `queries` has 0..60 entries, each exactly `{"query_id", "as_of"}` with a
  unique `common.ident` id and an integer `as_of` in `0..len(batches)`.
  Order need not increase and the same prefix may repeat.  `as_of = k` is
  strictly the view knowable after the first `k` batches; rejected batches
  still occupy a sequence number but have no state effect.

`common.read`/`common.save` still apply: at most 8 MB, duplicate JSON keys and
non-finite numbers rejected, `source_sha256` bound to the exact source bytes,
refusal of a symlink or same-path output, and atomic successful replacement.

### Calibration and observation semantics

Calibrations of one instrument take effect on the half-open interval
`[start, end)`.  The curve interpolates exactly with rational arithmetic
between strictly increasing raw nodes:

```
y = y0 + (x - x0) * (y1 - y0) / (x1 - x0)
```

End nodes belong to the domain; a raw value exactly equal to an internal node
returns that node's value.  Adjacent valid intervals may meet with the first
interval's `end` equal to the next one's `start`, and that time belongs only
to the later interval.  Within one committed view, valid intervals of the same
instrument must not overlap; different instruments are unconstrained.  An
observation whose time no interval covers is `unavailable` with reason
`no_calibration` (`calibration_id` and `value` null).  If a covering
calibration exists but the raw value lies outside its node domain, the
observation is `unavailable` with reason `outside_curve`, the selected
`calibration_id` still reported and `value` null.  A retroactive correction may
split, merge or reshape a past calibration; no calibration receives priority.

### Batch acceptance and protected certificates

Each batch applies all its puts and deletes simultaneously to a candidate view
on top of the last accepted state — `put` replaces the whole record of an
existing id, otherwise creates it — and only the final candidate is validated,
so the result never depends on change ordering.  The business rejection
reasons, in fixed output priority, are:

1. `missing_delete` — a delete names an absent key of that entity;
2. `overlap` — two same-instrument calibration intervals overlap;
3. `missing_member` — a certificate references an absent observation;
4. `formula_mismatch` — a member's parsed element-count dictionary differs
   from the certificate target's (`H2O` and `OH2` match; `A` and `A2` never
   reduce to each other);
5. `unavailable_member` — a member observation has no value in this view;
6. `mean_out_of_range` — the exact weighted mean falls outside
   `[lower, upper]`;
7. `range_exceeded` — `max(value) - min(value)` exceeds `max_range`.

With member values `v_i` and positive weights `w_i`,
`mean = sum(w_i * v_i) / sum(w_i)` and `range = max(v_i) - min(v_i)`.
Members, weights and thresholds are never auto-adjusted.  If any reason
applies, the whole batch — all three entity kinds — rolls back to the prior
state, processing continues with later batches, and no partial modification is
published.  Certificates untouched by the batch are still rechecked, so a
retroactive calibration that would invalidate an active certificate is
rejected unless the same batch also withdraws or replaces that certificate and
the final view is fully legal.

### Report

The report keeps the three top-level keys `schema_version` (fixed at 1),
`source_sha256` and `result`.  `result` is exactly `{"batches", "queries"}`.

- `batches` follows source order; each entry is exactly
  `{"batch_id", "status", "reason"}` with `status` `accepted` or `rejected`
  and `reason` null on acceptance or the first applicable reason code.
- `queries` follows input order; each entry is exactly
  `{"query_id", "as_of", "observations", "certificates"}`.  Both entity
  arrays sort by id in Python (Unicode code point) order and list only
  entities existing at that prefix; the empty state yields empty arrays.
  No calibration catalogue snapshot is emitted — each observation's source
  calibration is the evidence.
- An observation item has exactly the keys `id`, `instrument`, `time`,
  `formula`, `raw`, `status`, `calibration_id`, `value`, `reason`.  The first
  five mirror the view's record with `raw` canonicalised; `time` stays a JSON
  integer.  On availability `calibration_id` is the selected calibration,
  `value` is the exact interpolated value and `reason` is null; otherwise the
  two unavailable forms above apply.
- A certificate item has exactly the keys `id`, `formula`, `lower`, `upper`,
  `max_range`, `mean`, `range`, `members`.  Members sort by
  `observation_id`; each is exactly
  `{"observation_id", "weight", "value", "calibration_id"}` with values taken
  from the query's own view.  Every reported certificate passes in that view;
  rejected candidates never appear.
- Every output `raw`, `value`, `weight`, `lower`, `upper`, `max_range`,
  `mean` and `range` is a simplest `Fraction` string (`"0"`, `"4"`,
  `"1/3"`): no JSON decimals, unreduced fractions or negative zero.
  `time`/`start`/`end`/`as_of` remain JSON integers.

### Algorithm and state

The parser validates every batch and query before accounting and converts all
quantities to immutable `Fraction` records and all formulas to parsed
element-count tuples.  Processing keeps one mutable working state of three
entity dictionaries; each batch builds shallow candidate dictionaries, applies
the simultaneous changes, sorts same-instrument intervals to detect any
overlap, then re-evaluates every active certificate against an indexed
calibration view.  On rejection the candidate is discarded and the prior
dictionaries remain the state, so rejected and accepted batches, and every
historical query, always refer to one consistent committed view.

Only the prefixes actually named by queries are retained as snapshots, each a
tuple of `(id, record)` references over the immutable records.  A later
accepted retroactive correction therefore cannot alter an earlier prefix's
values, availability, calibration source or certificate membership, while
memory stays bounded by the at most 61 retained views rather than by 300 full
state copies.  An independent reference implementation used by the tests
(`tests/reference_certify.py`) re-derives every prefix directly from the empty
state with its own interpolation and certificate checks, never calling the
production evaluation functions, and the tests additionally permute the
changes inside multi-change batches to confirm order independence.

### Examples

`fixtures/certify-example.json` is the documented calibration/withdrawal
walkthrough.  Batch 1 installs the identity curve `c` on `[0, 20)`,
observation `o` (time 10, raw 4, formula `H2O`) and certificate `s` targeting
the same element counts via `OH2` with mean pinned to 4; it accepts.  Batch 2
replaces `c` by a two-times curve, which would move `o` to 8 and invalidate
`s`, so the batch is rejected (`mean_out_of_range`) and the identity state
survives.  Batch 3 applies the same replacement together with deleting `s` and
accepts.  Queries at prefixes 1 and 2 still report value 4 and the original
certificate with source `c`; prefix 3 reports value 8 and no certificate;
prefix 0 is empty.

```sh
reaction-certify --source fixtures/certify-example.json --output run/certify-example.json
```

