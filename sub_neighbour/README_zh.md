# sub_neighbour.py 中文使用教程

本教程介绍如何在 `sub_neighbour/sub_neighbour.py` 中使用交互式模式与命令行参数模式，对给定结构的目标原子进行 1NN/2NN/3NN 等壳层邻居识别，并导出结果。

## 1. 脚本概览
- 作用：选择“目标原子”（支持索引或元素名的“第一个出现”），按距离识别近邻分壳（1NN、2NN、3NN …）。
- 输入来源：
  - 结构文件：POSCAR/CONTCAR/CIF 等。
- 分壳方式：
  - 指定 `cutoffs`（每一壳的距离上限，单位 Å）。
  - 自动分壳（指定壳层数 `shells` 与最大半径 `max-radius`，按距离突变阈值 `min-dist-gap` 划分）。

## 2. 运行环境与依赖
- Python 3.8+（推荐 Anaconda 环境）
- 必需：`pymatgen`
- 可选：`PyYAML`（当选择 `yaml` 输出时更优；未安装会回退为 JSON）

安装示例：
```bash
pip install pymatgen
# 如需 YAML 输出
pip install pyyaml
```

## 3. 路径位置
- 脚本：`C:/Users/lenovo/Desktop/缺陷产生脚本/sub_neighbour/sub_neighbour.py`
- 教程：`C:/Users/lenovo/Desktop/缺陷产生脚本/sub_neighbour/README_zh.md`

## 4. 交互式模式（推荐上手）
当不带参数执行脚本时，会自动进入交互式模式，逐步引导你输入（仅支持“结构文件”作为输入来源），并且支持“两次替位 + 按 xNN 生成替位结构”。

启动：
```bash
cd C:/Users/lenovo/Desktop/缺陷产生脚本
python sub_neighbour/sub_neighbour.py
```

交互问题与示例回答（示例：用 POSCAR，先替位 Zr，再在 1nn 处替位 Pu，自动分 3 壳）：
1) 结构文件路径（可为目录，自动识别 POSCAR/CONTCAR/CIF）：`POSCAR`
2) 选择索引基：`0` 或 `1`：`0`
3) 请输入要被替位的目标原子索引：`0`
4) 第一次替位元素（如 `Zr`）：`Zr`
5) 第二次替位元素（如 `Pu`）：`Pu`
6) 分壳方式：`1`（cutoffs）或 `2`（自动 shells）：`2`
7) shells（默认 3）：回车
8) max-radius（Å，默认 8.0）：回车
9) min-dist-gap（Å，默认 0.15）：回车
10) 输出格式：`json` / `yaml` / `txt`：回车（默认 json）
11) 输出文件路径（留空打印到终端）：`neighbour_report.json`
12) 对称等价归并？`y/n`：`n`
13) 打印详细信息？`y/n`：`y`
14) 请选择要生成的 xNN（输入 `1`/`2`/`3` 或 `1,2` 或 `all`）：`1`
15) 是否对该壳的全部位置生成？`y/n`：`y`

说明：
- 上述步骤会在当前目录下创建形如 `Zr-Pu-1nn/` 的文件夹，并在其中写出多个 `POSCAR_index-*.vasp`（对应第二个替位原子选择的不同邻居位置）。
- 同时会写出 `mapping.json`，记录生成的文件与对应的结构索引、距离等信息，便于追踪。
 - 邻居搜索已固定为“仅原胞内”（不跨边界、无周期像）。
 - 替位结构的输出目录定位为“输入结构文件所在目录”的子目录，例如：若结构文件为 `D:/proj/POSCAR`，输出目录则为 `D:/proj/Zr-Pu-1nn/`。

## 5. 命令行参数模式（适合脚本化/批处理）
查看帮助：
```bash
python sub_neighbour/sub_neighbour.py --help
```

示例1：元素 + cutoffs 分壳
```bash
python sub_neighbour/sub_neighbour.py \
  --structure POSCAR \
  --target Si \
  --cutoffs 2.5,3.5,4.5 \
  --output neighbour_report.json \
  --verbose
```

示例2：索引 + 自动分壳（3 壳）
```bash
python sub_neighbour/sub_neighbour.py \
  --structure POSCAR \
  --target 0 \
  --shells 3 \
  --max-radius 8.0 \
  --format yaml \
  --verbose
```

（已移除）不再支持在 CLI 中通过 `supercell_info.json` 输入（交互模式亦默认使用结构文件）。

### 关键参数说明（核心）
- `--structure`：结构文件路径（可给目录，自动寻找 POSCAR/CONTCAR/CIF）
- `--target`：
  - 数字：目标原子索引；配合 `--index-base 0|1`
  - 元素名：如 `Si`，选择“第一个出现的 Si 原子”
- 分壳方式二选一：
  - `--cutoffs a,b,c`：明确每壳上限（Å）
  - `--shells N` 与 `--max-radius R`：在 [0, R] 自动分 N 壳；`--min-dist-gap`（默认 0.15 Å）为壳层边界阈值
- `--format`：`json|yaml|txt`（默认 `json`）
- `--output`：输出文件路径；未给则打印到终端
- `--verbose`：打印更多过程信息

## 6. 输入结构文件说明
脚本将自动识别常见文件名（如 `POSCAR`、`CONTCAR`、`*.cif` 等）。若提供目录路径，会在目录内尝试查找这些常见文件名。

## 7. 输出说明
邻居统计输出包含：
- `target_site`：目标原子 index、element、frac_coords、cart_coords
- `shells`（每一壳）：
  - `shell`（壳编号，从 1 开始）
  - `cutoff`（本壳距离上限/代表值）
  - `count`（原子数）
  - `items`：每个邻居的 index、element、distance、image（周期像）、delta_cart（位移矢量）
- `stats`：`total_neighbors`（总邻居数）、`element_counts`（元素分布）

替位结构导出：
- 目录命名：`First-Second-xnn`，例如 `Zr-Pu-1nn/`。
- 文件命名：`POSCAR_index-<j>.vasp`（j 为结构索引，遵循所选的索引基显示）。
- 元数据：`mapping.json`，记录每个导出文件对应的邻居信息（索引、距离、原始元素等）。
 - POSCAR 元素行顺序：基于“原结构元素（去除将替入的元素）+ [第二次替位元素, 第一次替位元素]”的唯一序列，并使用 `sort_structure=False` 保留站点顺序。例如原结构含 `U`，第一次替位 `Pu`、第二次替位 `Zr`，则元素行顺序为 `U Zr Pu`。

## 8. 常见问题与排错
- 中文帮助输出报 `UnicodeEncodeError`：脚本已内置 Windows 控制台编码兜底设置（UTF-8 + errors=replace）。
- 提示缺少输入：
  - 交互模式会引导选择来源；命令行模式需至少提供 `--structure`。
- 目标元素不存在/索引越界：检查结构的元素种类与原子总数，并确认 `--index-base`（0 或 1）。
- 自动分壳效果不理想：调整 `--min-dist-gap`（如 0.2~0.3 Å）、增大 `--max-radius`，或改用明确的 `--cutoffs`。
- 邻居为空：`max-radius` 过小，请增大（如 8.0~10.0 Å）。
 - 不跨边界说明：已强制仅原胞内搜索邻居；若你确实需要考虑周期像，可告知，我们可提供一个可选开关（如 `--allow-pbc`）。
 - POSCAR 元素顺序：若你希望固定为某个全局顺序（例如始终 `U Zr Pu O ...`），可联系补充一个交互/CLI 选项。

## 9. 建议工作流示例
- 与其他流程配合（基于结构文件）：
```bash
python sub_neighbour/sub_neighbour.py --structure POSCAR --target O --shells 3
```
- 快速配位检查（以 Si 金刚石为例）：
```bash
python sub_neighbour/sub_neighbour.py --structure POSCAR --target Si --shells 3 --verbose
```

## 10. 后续可选增强
- 对称等价归并（`--symmetry` 开关当前为占位）：可基于空间群将等价位归并并统计多重度。
- 导出 CSV，便于 Excel 浏览。
- 增加“仅同胞像/合并不同像”的筛选项。
- 增加配位环境识别与价态过滤。
