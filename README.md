# data-qd-160 · 化学反应配料核算

Exact formula balancing and limiting-reagent batch accounting。业务输入、计算和当前边界详见 [数据契约](docs/contract.md)。

## 安装与演示

Python >=3.11，运行和测试仅用标准库；先确认 `python3 --version`。
每仓建立独立环境，构建工具锁版本。首次安装无公共缓存时可能联网；业务运行不访问网络。

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-build.lock
.venv/bin/python -m pip install --no-build-isolation --no-deps .
.venv/bin/python -m pip check
.venv/bin/reaction-batch --source fixtures/input.json --output run/report.json
.venv/bin/reaction-lots --source fixtures/lots-example.json --output run/lots-example.json
.venv/bin/python -m unittest discover -s tests -v
```

本地已有wheel时，安装构建工具可加 `--no-index --find-links <公共wheel缓存目录>`。
修改源码后重新安装。本仓不共享其他仓环境，不依赖Docker、外部账号或持续服务。

输入JSON最多8MB，拒绝重复键、非有限数值和未定义字段；源字节SHA256随结果保留。
输出是可独立读取的JSON；成功后原子替换，输入/计算失败退出2并保留旧报告。
输入和输出不得相同，输出符号链接拒绝。单操作者批处理，不承诺并发发布或断电恢复。
同输入重跑不累加。测试使用临时目录，运行产物在忽略的run/。

src/reaction_batch/model.py 为原批次业务模型，lot_model.py 为投料/隔离账本；common.py 提供受限JSON及报告发布，
cli.py 与 lots.py 衔接两个命令的数据链路。fixtures 为明确自制样例，代码与样例采用MIT许可。

## 批次投料与隔离（`reaction-lots`）

`reaction-lots --source INPUT --output REPORT` 处理 `{"operations":[...]}` 有序事实流：
`stock` 建 lot，`run` 按同式合并、既有严格正配平规则扣料并生成新产品 lot 与直接父子边，
`hold`/`release` 管理活动隔离根。被隔离集合为所有活动根沿成功 run 谱系的后代并含根自身；
`blocked` 由谱系而非当前库存推导。报告仅含最终 `operations` 状态、按标识排序的 `lots`
（`remaining_moles` 为最简 Fraction 字符串）与 `active_holds`。

业务拒绝（未知 lot、库存不足、被隔离投料、无法配平或不同一进度、身份重用等）记为
`rejected` 并保持事发前状态；结构/字段/字面值/上限错误整批退出 2 并保留旧报告。
无输出祖先快照、无全对可达矩阵、遍历非递归；2 万 lot 单链与 6 万边汇合图可在
192MiB 内完成。完整契约见 [docs/contract.md](docs/contract.md)。

```sh
.venv/bin/reaction-lots --source fixtures/lots-diamond.json --output run/lots-diamond.json
```

## 防泄漏训练/测试切分（`reaction-split`）

`reaction-split --source INPUT --output REPORT` 为化学实验分类样本生成不泄漏的
训练/测试清单。元素计数完全相同的化学式（如 `H2O` 与 `OH2`）或共享
`experiment_id` 的样本取无向传递闭包，每个分量不可拆分；闭包先于配额规划，
含相反 `pin` 的分量报 `pin_conflict`，否则在分量整体归属上联合满足两类测试
配额并最小化两类测试数到目标的绝对偏差之和，同分时按 `record_id` 升序的
0/1 归属向量取字典序最小。规划器按正例计数行保存负例可达位图，只保留有界
计数状态与规范前缀，不枚举分量子集。

```sh
.venv/bin/reaction-split --source fixtures/split-example1.json --output run/split1.json
.venv/bin/reaction-split --source fixtures/split-example2.json --output run/split2.json
.venv/bin/reaction-split --source fixtures/split-example3.json --output run/split3.json
.venv/bin/reaction-split --source fixtures/split-scale-600.json --output run/split600.json
```

输入严格为 `{"records":[...],"test_quotas":{...}}`：0..600 条样本、每标签最多
300 条、配额为 0..300 的 JSON 整数闭区间。原五个命令的输入、输出与失败行为
不变。完整契约见 [docs/contract.md](docs/contract.md)。

## 业务与样例

# Exact reaction batch material accounting v1
教学数据整理者输入reactants/products为各1..20个全局唯一化学式，available_moles恰含所有反应物与非负十进制字符串（整数最多9位、小数最多9位）。摩尔量用Fraction精确计算，不转换质量。
化学式支持[A-Z][a-z]?元素token、正整数下标<=1000和圆括号嵌套<=10层，长度<=100、每元素原子数<=10^12；不验证周期表符号真实性，不接受同位素、离子电荷、水合点或反应箭头。
构建原子守恒有理矩阵求零空间；必须恰一维并有全部严格正系数，约为最小互素正整数。欠定/不可正配平整批拒绝，不任意选解；不平衡电荷或判断实际反应可行性。
以available/coefficient最小值确定反应进度，可有并列限量试剂和零进度；输出消耗、剩余、产物理论摩尔上限及每进度原子守恒证据。样例2H2+O2→2H2O，5mol H2/2mol O2时进度2、H2余1、O2余0、H2O理论4mol。不是实验操作指南，不包含温压、动力学、安全性或真实产率。

## 测量一致性核对（`reaction-assay`）

`reaction-assay --source INPUT --output REPORT` 核对候选物种组成与重叠的元素摩尔测量区间：
把每个样品化为一个精确线性规划，判断是否存在共同组成，并给出每种组分可被数据识别的紧界。
报告为每个界给出可行组成见证与非负对偶证书，无法解释的样品给出精确 Farkas 矛盾向量。
与 `reaction-batch`、`reaction-lots` 相互独立，原输入、输出与退出码不变。

```sh
.venv/bin/reaction-assay --source fixtures/assay-example1.json --output run/assay1.json
.venv/bin/reaction-assay --source fixtures/assay-mixed-48.json --output run/assay48.json
```

输入严格为 `{"species":[...],"samples":[...]}`；候选物种 1..60、样品 1..8、
每样品测量 0..120。全部数值为最简 Fraction 字符串；相等区间视为精确等式，
`1.000000001` 与 `1` 的差异不会被容差吞掉。完整契约、证书定义与样例见
[docs/contract.md](docs/contract.md)。

## 共同时间漂移下的峰身份核对（`reaction-peaks`）

`reaction-peaks --source INPUT --output REPORT` 对一批色谱观测峰与参考条目作共同时间
漂移下的全局一对一身份核对：同一运行共享一个仪器漂移量，峰与参考一一配对。报告最大可
解释峰数 K、能达到 K 的全部漂移区间（闭区间、可含单点、可不相连）、每个峰在所有最优
解释中的身份支持及可否不被匹配，并给出规范见证（最优漂移最小值处字典序最小的
assignment 向量）。支持只来自同一合法漂移下的完整 K 匹配，不从跨漂移点的边拼出身份。

```sh
.venv/bin/reaction-peaks --source fixtures/peaks-example1.json --output run/peaks1.json
.venv/bin/reaction-peaks --source fixtures/peaks-scale-40.json --output run/peaks40.json
```

输入严格为 `{"references":[...],"runs":[...]}`；参考 1..40、运行 1..8、每运行峰 0..40。
漂移边界为可带一个前导负号的十进制字符串，时刻与窗口为非负十进制字符串；全部输出边界
为最简 Fraction 字符串。算法按漂移区间端点扫描并用增量最大匹配与交错路径结构求支持，
不枚举匹配。完整契约见 [docs/contract.md](docs/contract.md)。

## 分子同位素联合分布与选择性检出核算（`reaction-isotopes`）

`reaction-isotopes --source INPUT --output REPORT` 按用户提供的同位素丰度表计算分子的精确名义质量包络，以及按名义质量与标记原子数联合触发的通道检出概率、条件标记数均值/方差和未检出部分。
同位素质量数、丰度与 `tagged` 均来自输入，不内置自然丰度或真实元素表；每个原子独立抽样，同元素多个原子仍按位置区分，因此 `(M,T)` 概率是各原子赋值概率的累加，而不是两个边缘分布相乘。
通道按输入顺序排他归属，`missed` 含全部截断尾部，`overflow` 是 `missed` 的子集，不重复相加。
所有概率与矩为精确最简有理数字符串；查询相互独立；相同源字节重跑报告相同。
原 `reaction-batch`、`reaction-lots`、`reaction-assay`、`reaction-peaks` 接口保持兼容。

```sh
.venv/bin/reaction-isotopes --source fixtures/isotopes-example1.json --output run/isotopes1.json
.venv/bin/reaction-isotopes --source fixtures/isotopes-example2.json --output run/isotopes2.json
```

完整契约、通道排他规则、矩守恒与样例说明见 [docs/contract.md](docs/contract.md)。
