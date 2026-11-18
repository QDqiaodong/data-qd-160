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
