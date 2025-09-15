# 两原子定向最近邻替位脚本 使用教程（pymatgen）

更新时间：2025-09-04（已更新：仅保留 either 模式；新增“同胞内邻居”约束与 15° 硬性角度限制；输出目录命名规则更新）
脚本路径：sub-direction/dir_nn_substitute.py

## 1. 功能简介
本脚本用于在晶体基体中，生成沿指定晶体方向 [u v w] 的两个最近邻替位（Substitutional）缺陷：
- 选择首个替位位点（支持：按元素/按站点索引/按 Wyckoff 字母）
- 在该位点的第一近邻壳层中，沿给定方向筛选“前向最近邻”作为第二替位位点
- 同时替换两个位点元素（可相同或不同）
- 可选：原胞标准化、超胞扩展
- 导出 VASP POSCAR 与运行摘要 summary.json，并对最近邻与方向进行验证

关键依赖：pymatgen（Structure、Lattice、SpacegroupAnalyzer、ReplaceSiteSpeciesTransformation、get_all_neighbors、Poscar 等）

## 2. 环境准备
- Python 3.9+ 推荐
- 安装 pymatgen（二选一）：
  - Conda：`conda install -c conda-forge pymatgen`
  - Pip：`pip install pymatgen`
- 结构文件格式：POSCAR/CONTCAR、CIF、VASP5 POSCAR 等常见晶体文件

## 3. 脚本位置与运行
- 脚本：`sub-direction/dir_nn_substitute.py`
- 运行（Windows PowerShell/终端）：
  - `python sub-direction/dir_nn_substitute.py`
- 全程交互式，根据提示输入即可。

## 4. 交互参数说明
1) 结构文件路径：输入 POSCAR/CIF 等文件的路径。
2) 是否转换为原胞 (y/n)：使用 SpacegroupAnalyzer.find_primitive() 以获得一致原胞。
3) 是否扩展超胞：输入如 `1 1 1`（不变）或 `2 2 2`（扩胞）。
4) 首位点选择方式：
   - 1 按元素：输入元素符号（如 Si），脚本列出候选索引再选择
   - 2 按索引：直接输入站点索引（脚本会打印站点列表）
   - 3 按 Wyckoff：输入字母（如 a/b/c），脚本列出对应索引再选择
5) 晶体方向 [u v w]：例如 `1 0 0`、`1 1 0`、`1 1 1`（基于分数坐标方向，脚本会转换并归一化为笛卡尔方向）。
6) 两个替位元素 D1、D2：可相同或不同（如 Al 与 Mg）。
7) 参数：
   - 距离容差 tol（Å），默认 0.10：用于第一近邻壳层判定区间 [d1−tol, d1+tol]
   - 角度容差 θ_tol（度），默认 10：用于方向筛选的夹角阈值（仅 either 模式）
   - 自动回退 auto_fallback：默认 `y`，当在设定 θ_tol 内找不到匹配时，脚本会将角度放宽至 45° 再尝试（仅 either 模式）
8) 输出目录：默认目录名为 `<基体>-<D1>-<D2>-[uvw]`，例如 `U-Pu-Pu-[010]`，可自定义。

## 4.1 方向选择与自动回退（更新）

- 方向选择：仅保留 either 模式（不区分正/反向）。筛选规则为：
  1) 在角度阈值 θ_tol 内，优先选择与方向单位向量 d̂ 夹角最小的候选；
  2) 若最小角度相同，则选择距离更近的候选（保证是真正的最近邻）。
- 自动回退 auto_fallback：当在给定 θ_tol 内找不到满足条件的候选时，脚本会自动将角度容差放宽至 45° 再尝试一次。

选择建议：
- 为获得更稳健的“第一近邻”，建议保持较小的 tol（如 0.05–0.15 Å），θ_tol 可按 5–15° 视结构各向异性调节。


## 5. 工作流程与判据（要点）
- 使用 `Structure.get_all_neighbors(r)`（r ≥ max(a,b,c)）获取周期性最近像的邻居；仅保留同胞内候选（image == [0,0,0]），避免跨边界选择。
- 第一近邻壳层：计算最小非零距离 d1，保留距离在 [d1−tol, d1+tol] 的邻居。
- 方向筛选（either）：
  - 计算从首位点到候选邻居“最近像”的位移向量 Δr（利用邻居对象的 image）。
  - 要求与 d̂ 的夹角 ≤ θ_tol（不区分正/反向）。
  - 在满足角度阈值的候选中，先选夹角最小者；若夹角相同，选距离更小者。
- 替换：`ReplaceSiteSpeciesTransformation({i: D1, j: D2})` 同时替换两个站点。
- 验证：
  - 再次判定 j 是否在 i 的第一壳层且为同胞内候选（image == [0,0,0]），并检查方向夹角；结果记录于 summary.json。

## 6. 输出文件
- POSCAR：替位后的结构文件。
- summary.json：运行摘要，包括：
  - input_file、space_group、supercell、primary_index（首位点 i）、secondary_index（第二位点 j）
  - elements（D1/D2）、direction_uvw、tol_A、theta_tol_deg
  - neighbor_choice（被选中邻居的投影/角度/距离/使用的 image）
  - verification（是否第一壳层、角度是否满足、具体数值）
  - formula_before / formula_after（替位前/后化学式）

## 7. 常见问题与建议
- 未找到符合方向条件的最近邻：
  - 使用“自动回退”：将角度容差放宽至 45° 后再试一次（仅 either 模式）
  - 手动放宽 θ_tol（如 10° → 15°/20°/30°）
  - 根据晶体对称性与局域配位，调整 [u v w]
  - 微调 tol（0.05–0.15 Å 范围尝试），或扩大超胞以减小边界/镜像效应
- Wyckoff 显示 “-”：
  - 可能对称性解析失败；可先转换为原胞，或改为索引/元素选择
- 第二位点不是第一壳层：
  - 调整 tol 以稳定第一壳层识别
- 输出失败：
  - 确保目标目录可写或更换无特殊字符的路径

## 8. 扩展方向（可选开发）
- 多首位点/多方向的批量枚举与导出
- 对称等价去重（SpacegroupAnalyzer.are_symmetrically_equivalent 或 StructureMatcher）
- 自动化超胞选择（基于最短镜像距离阈值）
- 增加 CLI 非交互参数支持，便于脚本化批处理

## 9. 快速试跑清单
1) 安装 pymatgen，准备 POSCAR/CIF
2) 运行：`python sub-direction/dir_nn_substitute.py`
3) 按提示选择首位点、输入 [u v w]、D1/D2、tol/θ_tol、输出目录
4) 检查输出：POSCAR 与 summary.json
5) 若不满足预期：参考第 7 节调整参数重试

—— 完 ——

