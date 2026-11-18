本目录为原创虚构合成数据，使用MIT许可。没有真实个人、线上实验、生产日志、交易、项目进度或搜索用户信息；不来自外部平台。

- `input.json`：原 `reaction-batch` 样例（2H2 + O2 → 2H2O）。
- `lots-example.json`：`reaction-lots` 文档示例，两个同式 H2 lot 与 O2 合并投料，随后 hold 耗尽祖先。
- `lots-diamond.json`：菱形谱系双活动根与单根释放的隔离示例。
- `assay-example1.json`：`reaction-assay` 文档示例，唯一组成 `["1","1/2","1/2"]`。
- `assay-example2.json`：同元素组成不同化学式（H2O/OH2），两端点各有组成见证。
- `assay-example3.json`：`1.000000001` 与 `1` 的精确矛盾，返回 Farkas 证书。
- `assay-mixed-48.json`：48 物种多自由度、冗余与重复区间、联合矛盾、总量 0 的可行/不可行、空测量与未出现元素样例。
- `peaks-example1.json`：`reaction-peaks` 文档示例，仅单点漂移 1 同时容纳两条边，K=2、最优漂移仅 [1,1]。
- `peaks-example2.json`：两个分离的最优漂移区间 [-1,0] 与 [2,3]，峰 p 支持 a、b 且不能不匹配，规范见证 d=-1 配 b。
- `peaks-linkage.json`：B 组参考窗口变化联动改变 A 组峰 p 的身份支持，最优漂移为两个单点 0 与 1。
- `peaks-scale-40.json`：40 参考 × 40 峰全等价匹配，K=40、整段漂移域最优、每峰支持全部 40 个参考。
