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
