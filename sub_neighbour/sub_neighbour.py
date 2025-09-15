#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
sub_neighbour.py

功能：
- 在给定结构中，先确定一个目标原子（优先“第一个该元素”），随后基于距离识别 1NN、2NN、3NN ... 壳层邻居。
- 支持两种分壳方式：
  1) --cutoffs a,b,c  通过用户给定的每壳距离上限（Å）直接切分。
  2) --shells N --max-radius R  在 [0, R] 内自动分层，按距离突变（最小间隙阈值 --min-dist-gap）或截取前 N 壳。
- 支持从结构文件（POSCAR/CONTCAR/CIF/等）或（尽量兼容的）supercell_info.json 构建 pymatgen Structure。
- 输出 JSON（默认）/YAML/TXT 三种格式。

依赖：pymatgen（必须）。若选择 YAML 输出且未安装 PyYAML，则自动回退为 JSON 并给出警告。

示例：
- 指定元素，按 cutoffs：
  python sub_neighbour.py --structure POSCAR --target Si --cutoffs 2.5,3.5,4.5 --output neighbour_report.json
- 指定索引，自动分壳：
  python sub_neighbour.py --structure POSCAR --target 0 --shells 3 --max-radius 8.0 --format yaml
- 使用 supercell_info.json：
  python sub_neighbour.py --supercell-info supercell_info.json --target O --shells 3
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

try:
    import yaml as _yaml  # type: ignore
    _HAS_YAML = True
except Exception:
    _yaml = None
    _HAS_YAML = False

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.io.cif import CifParser
from pymatgen.core.lattice import Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer  # 可选，用于对称等价


# ========== Windows 控制台编码兼容 ==========
try:
    # 避免在 GBK 控制台下打印包含非 ASCII 字符时出现 UnicodeEncodeError
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    if hasattr(sys.stderr, "reconfigure"):
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

# ========== 简单日志 ==========
# ========== 简单日志 ==========

def info(msg: str) -> None:
    print(f"[INFO] {msg}")


def warn(msg: str) -> None:
    print(f"[WARN] {msg}")


def err(msg: str) -> None:
    print(f"[ERROR] {msg}")


# ========== 工具函数 ==========

def ask(prompt: str, default: Optional[str] = None) -> str:
    """交互式输入，支持默认值。返回去除首尾空白的字符串。"""
    if default is not None:
        prompt = f"{prompt} [默认: {default}]: "
    else:
        prompt = f"{prompt}: "
    try:
        s = input(prompt)
    except EOFError:
        s = ""
    s = (s or "").strip()
    return s if s else (default or "")

def _to_float_list(x: Sequence[Any]) -> List[float]:
    return [float(v) for v in x]


def _neighbor_distance(n: Any) -> float:
    # 兼容 PeriodicNeighbor/Neighbor 的不同属性名
    if hasattr(n, "nn_distance"):
        return float(getattr(n, "nn_distance"))
    if hasattr(n, "distance"):
        return float(getattr(n, "distance"))
    # 兜底计算
    try:
        return float(np.linalg.norm(np.array(n.vec)))
    except Exception:
        raise ValueError("未知邻居对象，无法获取距离")


def parse_cutoffs(s: str) -> List[float]:
    parts = [p for p in s.replace(";", ",").split(",") if p.strip()]
    vals = [float(p.strip()) for p in parts]
    if any(v <= 0 for v in vals):
        raise ValueError("cutoffs 中的距离必须为正数")
    return vals


# ========== 结构读取 ==========

def load_structure_from_path(path: str) -> Structure:
    # 若为目录，尝试自动识别常见文件名
    if os.path.isdir(path):
        candidates = [
            "POSCAR", "CONTCAR", "POSCAR.vasp", "CONTCAR.vasp",
            "structure.cif", "input.cif",
        ]
        for name in candidates:
            p = os.path.join(path, name)
            if os.path.isfile(p):
                info(f"检测到目录，使用文件: {p}")
                path = p
                break
        else:
            raise FileNotFoundError(f"目录中未找到常见结构文件: {path}")

    if not os.path.isfile(path):
        raise FileNotFoundError(f"未找到结构文件: {path}")

    info(f"读取结构文件: {path}")

    # 1) 优先自动识别
    try:
        return Structure.from_file(path)
    except Exception:
        pass

    # 2) 尝试 POSCAR
    try:
        poscar = Poscar.from_file(path)
        return poscar.structure
    except Exception:
        pass

    # 3) 尝试 CIF
    try:
        parser = CifParser(path)
        structs = parser.get_structures()
        if structs:
            return structs[0]
    except Exception:
        pass

    # 4) 兜底：文本按 POSCAR
    try:
        with open(path, "r", encoding="utf-8") as f:
            txt = f.read()
        return Structure.from_str(txt, fmt="poscar")
    except Exception as e:
        raise ValueError(f"无法解析结构文件: {path}，请提供标准 POSCAR/CONTCAR 或 .cif。原始错误: {e}")


def load_structure_from_supercell_info(path: str) -> Structure:
    """尽力从 supercell_info.json 构建 Structure。
    期望字段（尽量兼容，多路尝试）：
    {
      "lattice": [[a1],[a2],[a3]]  # 3x3，单位 Å
      "sites": [
         {"element": "Si", "frac_coords": [x,y,z]}  或
         {"specie": "Si", "f": [x,y,z]} 或
         {"symbol": "Si", "frac": [x,y,z]}
      ]
    }
    若失败，抛出异常并提示改用 --structure。
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"未找到 supercell_info.json: {path}")
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)
    # lattice
    lat = None
    for key in ("lattice", "latt", "cell", "matrix"):
        if key in data:
            mat = data[key]
            if isinstance(mat, (list, tuple)) and len(mat) == 3:
                lat = Lattice(mat)
                break
    if lat is None:
        raise ValueError("supercell_info.json 缺少 lattice 字段或格式错误")
    # sites
    sites = data.get("sites") or data.get("atoms") or data.get("positions")
    if not sites or not isinstance(sites, list):
        raise ValueError("supercell_info.json 缺少 sites/atoms/positions 列表")
    species: List[str] = []
    fcoords: List[List[float]] = []
    for s in sites:
        if not isinstance(s, dict):
            raise ValueError("sites 条目应为对象")
        elem = s.get("element") or s.get("specie") or s.get("symbol") or s.get("type")
        if elem is None:
            raise ValueError("sites 条目缺少元素字段 element/specie/symbol/type")
        fc = s.get("frac_coords") or s.get("f") or s.get("frac") or s.get("fract_coords")
        if fc is None:
            # 兼容给出笛卡尔坐标的情况
            cart = s.get("cart_coords") or s.get("c") or s.get("cart")
            if cart is None:
                raise ValueError("sites 条目缺少分数坐标 frac_coords/f，或缺少笛卡尔坐标 cart_coords/cart")
            # 转分数
            fc = lat.get_fractional_coords(cart)
        species.append(str(elem))
        fcoords.append(_to_float_list(fc))
    return Structure(lattice=lat, species=species, coords=fcoords, coords_are_cartesian=False)


def load_structure(structure_path: Optional[str], supercell_info_path: Optional[str]) -> Structure:
    if structure_path:
        return load_structure_from_path(structure_path)
    if supercell_info_path:
        warn("未提供 --structure，尝试从 supercell_info.json 构建结构……")
        return load_structure_from_supercell_info(supercell_info_path)
    raise ValueError("必须提供 --structure 或 --supercell-info 之一")


# ========== 站点列表展示 ==========
def list_sites_simple(struct: Structure, show_cart: bool = False) -> None:
    """打印结构中所有站点：索引、元素、分数坐标（可选笛卡尔）。"""
    header = "索引  元素   分数坐标 (f.)                    "
    if show_cart:
        header += "笛卡尔坐标 (Å)"
    print(header)
    print("-" * len(header))
    for i, site in enumerate(struct.sites):
        f = site.frac_coords
        line = f"{i:>4d}  {str(site.specie):>3s}   [ {f[0]:8.5f} {f[1]:8.5f} {f[2]:8.5f} ]"
        if show_cart:
            c = site.coords
            line += f"   ( {c[0]:9.5f} {c[1]:9.5f} {c[2]:9.5f} )"
        print(line)


# ========== 目标原子选择 ==========

def select_target_site(struct: Structure, target: str, index_base: int) -> int:
    """target 可以是整数索引（字符串形式）或元素符号。
    - 若是纯数字：按 index_base 换算。
    - 否则：按元素符号，选择结构中该元素“第一次出现”的原子。
    """
    t = target.strip()
    if t.isdigit() or (t.startswith("-") and t[1:].isdigit()):
        idx = int(t)
        if index_base == 1:
            idx -= 1
        if idx < 0 or idx >= len(struct):
            raise ValueError(f"目标索引越界：{idx}，结构共有 {len(struct)} 个原子")
        return idx
    # 按元素名
    cand = [i for i, s in enumerate(struct.sites) if str(s.specie) == t]
    if not cand:
        # 提示可用元素
        elems = {}
        for s in struct.sites:
            e = str(s.specie)
            elems[e] = elems.get(e, 0) + 1
        raise ValueError(f"结构中不存在元素 {t}。可用元素及计数：{elems}")
    return cand[0]


# ========== 邻居收集与分壳 ==========

@dataclass
class NeighborItem:
    index: int
    element: str
    distance: float
    image: Tuple[int, int, int]
    delta_cart: Tuple[float, float, float]


def collect_neighbors(struct: Structure, target_index: int, rmax: float, use_pbc: bool = True) -> List[NeighborItem]:
    """收集近邻。

    当 use_pbc 为 True 时，使用 pymatgen 的 get_neighbors（包含周期像）。
    当 use_pbc 为 False 时，仅在原胞内（不跨边界）按分数坐标直接差值计算距离，禁止使用周期像。
    """
    site = struct[target_index]
    lat = struct.lattice
    items: List[NeighborItem] = []

    if use_pbc:
        # 标准 PBC 邻居（可能跨边界）
        neighs = struct.get_neighbors(site, rmax, include_index=True)
        for n in neighs:
            j = int(n.index)
            if j == target_index and _neighbor_distance(n) <= 1e-8:
                continue
            pj = struct[j]
            try:
                im = tuple(int(x) for x in n.image)
            except Exception:
                im = tuple(n.image)
            delta_f = pj.frac_coords + np.array(im, dtype=float) - site.frac_coords
            delta = lat.get_cartesian_coords(delta_f)
            dist = float(np.linalg.norm(delta))
            if dist <= 1e-8:
                continue
            items.append(NeighborItem(
                index=j,
                element=str(pj.specie),
                distance=dist,
                image=(int(im[0]), int(im[1]), int(im[2])),
                delta_cart=(float(delta[0]), float(delta[1]), float(delta[2]))
            ))
        items.sort(key=lambda x: x.distance)
        return items

    # 非 PBC：仅原胞内。直接用分数坐标差，不考虑任何周期像。
    pi_f = site.frac_coords
    for j, pj in enumerate(struct.sites):
        if j == target_index:
            continue
        delta_f = np.array(pj.frac_coords, dtype=float) - np.array(pi_f, dtype=float)
        delta = lat.get_cartesian_coords(delta_f)
        dist = float(np.linalg.norm(delta))
        if dist <= 1e-8 or dist > rmax + 1e-8:
            continue
        items.append(NeighborItem(
            index=int(j),
            element=str(pj.specie),
            distance=dist,
            image=(0, 0, 0),
            delta_cart=(float(delta[0]), float(delta[1]), float(delta[2]))
        ))
    items.sort(key=lambda x: x.distance)
    return items


@dataclass
class Shell:
    shell: int
    cutoff: float  # 该壳的上限（或代表性距离）
    items: List[NeighborItem]


def split_shells_by_cutoffs(neighs: List[NeighborItem], cutoffs: List[float]) -> List[Shell]:
    if not cutoffs:
        return []
    shells: List[Shell] = []
    start = 0
    for k, c in enumerate(cutoffs, start=1):
        grp: List[NeighborItem] = []
        while start < len(neighs) and neighs[start].distance <= c + 1e-8:
            grp.append(neighs[start])
            start += 1
        shells.append(Shell(shell=k, cutoff=float(c), items=grp))
    return shells


def split_shells_auto(neighs: List[NeighborItem], n_shells: int, rmax: float, min_gap: float) -> List[Shell]:
    if not neighs:
        return []
    # 限制在 rmax 内
    neighs = [n for n in neighs if n.distance <= rmax + 1e-8]
    if not neighs:
        return []
    dists = [n.distance for n in neighs]
    # 找边界：相邻距离差 >= min_gap 视为壳层分界
    gaps = [dists[i+1] - dists[i] for i in range(len(dists)-1)]
    boundaries: List[int] = [i for i, g in enumerate(gaps) if g >= min_gap]
    # 构造 shells，最多 n_shells
    shells: List[Shell] = []
    start = 0
    shell_idx = 1
    for b in boundaries:
        if shell_idx > n_shells:
            break
        end = b + 1  # 切片上界
        grp = neighs[start:end]
        if grp:
            shells.append(Shell(shell=shell_idx, cutoff=max(x.distance for x in grp), items=grp))
            shell_idx += 1
        start = end
        if shell_idx > n_shells:
            break
    # 末尾残余作为最后一壳
    if shell_idx <= n_shells and start < len(neighs):
        grp = neighs[start:]
        shells.append(Shell(shell=shell_idx, cutoff=max(x.distance for x in grp), items=grp))
    # 仅保留前 n_shells
    shells = shells[:n_shells]
    return shells


# ========== 对称等价（可选） ==========

def symmetry_grouping(struct: Structure, shells: List[Shell]) -> List[Shell]:
    """当前实现为占位：仅原样返回。若后续需要，可基于 SpacegroupAnalyzer 做等价归并。"""
    return shells


# ========== 输出渲染 ==========

def build_summary(struct: Structure, target_index: int, shells: List[Shell]) -> Dict[str, Any]:
    target = struct[target_index]
    target_info = {
        "index": int(target_index),
        "element": str(target.specie),
        "frac_coords": [float(x) for x in target.frac_coords],
        "cart_coords": [float(x) for x in target.coords],
    }
    out_shells: List[Dict[str, Any]] = []
    element_counts: Dict[str, int] = {}
    total_neighbors = 0
    for s in shells:
        items = []
        for it in s.items:
            items.append({
                "index": it.index,
                "element": it.element,
                "distance": it.distance,
                "image": list(it.image),
                "delta_cart": list(it.delta_cart),
            })
            element_counts[it.element] = element_counts.get(it.element, 0) + 1
            total_neighbors += 1
        out_shells.append({
            "shell": s.shell,
            "cutoff": s.cutoff,
            "count": len(s.items),
            "items": items,
        })
    summary = {
        "target_site": target_info,
        "shells": out_shells,
        "stats": {
            "total_neighbors": total_neighbors,
            "element_counts": element_counts,
        }
    }
    return summary


# ========== 替位与导出工具 ==========

def _ensure_dir(path: str) -> None:
    try:
        os.makedirs(path, exist_ok=True)
    except Exception:
        pass


def _replace_species(struct: Structure, index: int, new_element: str) -> Structure:
    """返回一个新结构：将 index 位置替换为 new_element。"""
    new_struct = struct.copy()
    new_struct[index] = new_element
    return new_struct


def _unique_in_order(seq: Sequence[str]) -> List[str]:
    seen = set()
    out: List[str] = []
    for s in seq:
        if s not in seen:
            seen.add(s)
            out.append(s)
    return out


def _reorder_structure_by_element_order(struct: Structure, element_order: List[str]) -> Structure:
    """按给定元素顺序对站点进行分组重排，返回新的 Structure。

    - 仅按符号匹配（如 "U", "Zr", "Pu"）。
    - 若结构中存在 element_order 未包含的元素，会附加在末尾（保持其首次出现顺序）。
    """
    # 目标顺序 + 余下元素的顺序
    existing = [str(s.specie) for s in struct.sites]
    rest = [e for e in _unique_in_order(existing) if e not in element_order]
    final_order = element_order + rest

    species: List[str] = []
    fcoords: List[List[float]] = []
    for el in final_order:
        for site in struct.sites:
            if str(site.specie) == el:
                species.append(el)
                fcoords.append([float(x) for x in site.frac_coords])
    return Structure(lattice=struct.lattice, species=species, coords=fcoords, coords_are_cartesian=False)


def export_substitution_variants(
    struct: Structure,
    first_index: int,
    first_element: str,
    candidates: List[NeighborItem],
    second_element: str,
    out_dir: str,
    index_base: int = 0,
) -> List[str]:
    """为一组候选邻居位置生成替位结构并输出到 POSCAR 文件。

    生成流程：
    1) 先在 first_index 替换 first_element；
    2) 对于每一个候选邻居位置 j，再替换为 second_element；
    3) 输出到 out_dir 下的多个 POSCAR_* 文件。
    返回写出的文件路径列表。
    """
    _ensure_dir(out_dir)
    written: List[str] = []
    # 基础元素顺序：以原始结构的元素首次出现顺序为基底
    base_order = _unique_in_order([str(s.specie) for s in struct.sites])
    # 规则：
    # - 保留原有元素的顺序，但移除将被新引入的元素（若原结构不存在则无影响）；
    # - 新引入元素顺序优先按 [second_element, first_element]，以满足示例期望 [U, Zr, Pu]。
    base_core = [e for e in base_order if e not in {first_element, second_element}]
    desired_order = base_core + _unique_in_order([second_element, first_element])

    for it in candidates:
        j = int(it.index)
        s1 = _replace_species(struct, first_index, first_element)
        s2 = _replace_species(s1, j, second_element)
        # 依据 desired_order 重排站点，确保 POSCAR 元素序列为 [U, Zr, Pu] 这类唯一有序序列
        s2_ordered = _reorder_structure_by_element_order(s2, desired_order)
        # 显示索引按 index_base
        disp_j = j + 1 if index_base == 1 else j
        fname = f"POSCAR_index-{disp_j}.vasp"
        fpath = os.path.join(out_dir, fname)
        # 关闭结构内部自动排序，保留我们重排后的站点顺序
        Poscar(s2_ordered, sort_structure=False).write_file(fpath)
        written.append(fpath)
    # 额外写出一个 mapping.json 以便记录索引映射
    try:
        mapping = {
            "index_base": index_base,
            "first_index_display": first_index + 1 if index_base == 1 else first_index,
            "candidates": [
                {
                    "neighbor_index_display": (int(it.index) + 1 if index_base == 1 else int(it.index)),
                    "neighbor_index_0based": int(it.index),
                    "distance": float(it.distance),
                    "element_before": it.element,
                    "file": f"POSCAR_index-{(int(it.index) + 1 if index_base == 1 else int(it.index))}.vasp",
                }
                for it in candidates
            ],
        }
        with open(os.path.join(out_dir, "mapping.json"), "w", encoding="utf-8") as f:
            json.dump(mapping, f, ensure_ascii=False, indent=2)
    except Exception:
        pass
    return written


def dump_output(data: Dict[str, Any], fmt: str, out_path: Optional[str]) -> None:
    text = None
    if fmt.lower() == "json":
        text = json.dumps(data, ensure_ascii=False, indent=2)
    elif fmt.lower() == "yaml":
        if not _HAS_YAML:
            warn("未安装 PyYAML，已回退为 JSON 输出")
            text = json.dumps(data, ensure_ascii=False, indent=2)
        else:
            text = _yaml.safe_dump(data, allow_unicode=True, sort_keys=False)  # type: ignore
    elif fmt.lower() == "txt":
        # 简要文本汇总
        lines: List[str] = []
        t = data.get("target_site", {})
        lines.append(f"Target index={t.get('index')} element={t.get('element')}")
        for s in data.get("shells", []):
            lines.append(f"Shell {s.get('shell')}: count={s.get('count')} cutoff={s.get('cutoff')}")
            for item in s.get("items", []):
                lines.append(
                    f"  idx={item['index']:>4d} el={item['element']:>2s} dist={item['distance']:.4f} "
                    f"img={item['image']} d=({item['delta_cart'][0]:.4f} {item['delta_cart'][1]:.4f} {item['delta_cart'][2]:.4f})"
                )
        text = "\n".join(lines)
    else:
        raise ValueError("不支持的输出格式，仅支持 json|yaml|txt")

    if out_path:
        # 确保输出目录存在
        try:
            import os as _os
            parent = _os.path.dirname(out_path)
            if parent and not _os.path.isdir(parent):
                _os.makedirs(parent, exist_ok=True)
        except Exception:
            pass
        # 自动追加默认文件名
        if _os.path.isdir(out_path):
            if fmt.lower() == "json":
                out_path = _os.path.join(out_path, "output.json")
            elif fmt.lower() == "yaml":
                out_path = _os.path.join(out_path, "output.yaml")
            elif fmt.lower() == "txt":
                out_path = _os.path.join(out_path, "output.txt")
        elif out_path.endswith(_os.path.sep):
            if fmt.lower() == "json":
                out_path = _os.path.join(out_path, "output.json")
            elif fmt.lower() == "yaml":
                out_path = _os.path.join(out_path, "output.yaml")
            elif fmt.lower() == "txt":
                out_path = _os.path.join(out_path, "output.txt")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(text)
        info(f"已写出: {out_path}")
    else:
        print(text)


# ========== 主流程 ==========

def main(argv: Optional[Sequence[str]] = None) -> int:
    # 若无参数，则进入交互式模式
    if argv is None:
        argv = sys.argv[1:]

    if len(argv) == 0:
        info("进入交互式模式（Interactive Mode）")
        # 1) 输入结构（仅支持结构文件）
        structure_path: Optional[str] = None
        structure_path = ask("请输入结构文件路径（可为目录，自动识别 POSCAR/CONTCAR/CIF）", default="POSCAR")
        if structure_path is not None:
            structure_path = structure_path.strip().strip('"').strip("'")

        # === 读取结构并展示站点 ===
        try:
            struct = load_structure(structure_path, None)
        except Exception as e:
            err(str(e))
            return 2
        # 解析实际使用的结构文件路径，用于后续输出目录定位
        def _resolve_structure_file(path: str) -> Optional[str]:
            if not path:
                return None
            p = path.strip().strip('"').strip("'")
            if os.path.isdir(p):
                for name in [
                    "POSCAR", "CONTCAR", "POSCAR.vasp", "CONTCAR.vasp",
                    "structure.cif", "input.cif",
                ]:
                    cand = os.path.join(p, name)
                    if os.path.isfile(cand):
                        return cand
                return None
            return p if os.path.isfile(p) else None

        resolved_struct_file = _resolve_structure_file(structure_path or "")
        base_output_dir = os.path.dirname(resolved_struct_file) if resolved_struct_file else os.getcwd()
        show_cart_ans = ask("是否同时显示笛卡尔坐标？y/n", default="n").lower().startswith("y")
        list_sites_simple(struct, show_cart=show_cart_ans)

        # 2) 目标选择（仅按索引选择要被替位的元素）
        index_base_str = ask("选择索引基：0 或 1", default="0")
        try:
            index_base = 1 if str(index_base_str).strip() == "1" else 0
        except Exception:
            index_base = 0
        idx_in = ask("请输入要被替位的目标原子索引")
        try:
            if not idx_in.strip():
                raise ValueError("未输入索引")
            idx_val = int(idx_in)
            target_index = idx_val - 1 if index_base == 1 else idx_val
            if target_index < 0 or target_index >= len(struct):
                raise ValueError("索引越界")
        except Exception:
            err("索引无效或越界，请重新运行并输入正确的数值索引")
            return 2

        # 3) 确认替位元素（两次替位）
        first_elem = ask("请输入第一次替位元素（例如 Zr）")
        if not first_elem:
            err("未输入第一次替位元素")
            return 2
        second_elem = ask("请输入第二次替位元素（例如 Pu）")
        if not second_elem:
            err("未输入第二次替位元素")
            return 2

        # 4) 分壳设置（用于确定第二个替位原子相对第一个原子的 xNN）
        sh_mode = ask("选择分壳方式：1) cutoffs  2) 自动(shells)", default="2")
        cutoffs_s: Optional[str] = None
        shells_n: Optional[int] = None
        max_radius: float = 8.0
        min_dist_gap: float = 0.15
        if sh_mode == "1":
            cutoffs_s = ask("请输入 cutoffs（以逗号分隔，如 2.5,3.5,4.5）")
        else:
            shells_n_str = ask("请输入壳层数 shells（留空默认 3）", default="3")
            try:
                shells_n = int(shells_n_str)
            except Exception:
                shells_n = 3
            max_radius_str = ask("请输入最大半径 max-radius（Å）", default="8.0")
            try:
                max_radius = float(max_radius_str)
            except Exception:
                max_radius = 8.0
            min_gap_str = ask("自动分层最小距离间隙 min-dist-gap（Å）", default="0.15")
            try:
                min_dist_gap = float(min_gap_str)
            except Exception:
                min_dist_gap = 0.15

        # 5) 其他
        fmt = ask("输出格式：json / yaml / txt", default="json").lower() or "json"
        out_path = ask("输出文件路径（留空则打印到终端）", default="") or None
        if out_path:
            out_path = out_path.strip().strip('"').strip("'")
        symm = ask("是否尝试对称等价归并？y/n", default="n").lower().startswith("y")
        verb = ask("是否打印详细信息？y/n", default="n").lower().startswith("y")

        # === 执行 ===（target_index 已由用户索引选择得到）

        # rmax & cutoffs
        rmax = float(max_radius)
        cutoffs: Optional[List[float]] = None
        if cutoffs_s is not None and cutoffs_s.strip():
            try:
                cutoffs = parse_cutoffs(cutoffs_s)
                rmax = max(cutoffs) + 1e-6
            except Exception as e:
                err(f"解析 cutoffs 失败: {e}")
                return 2
        elif shells_n is None:
            shells_n = 3

        if verb:
            info(f"目标原子 index={target_index}, 元素={struct[target_index].specie}")
            info(f"搜索半径 rmax={rmax} Å")

        # 仅在原胞内搜索邻居（不跨边界）
        neighs = collect_neighbors(struct, target_index, rmax, use_pbc=False)

        if cutoffs is not None:
            shells = split_shells_by_cutoffs(neighs, cutoffs)
        else:
            shells = split_shells_auto(neighs, n_shells=int(shells_n), rmax=rmax, min_gap=float(min_dist_gap))

        if symm:
            shells = symmetry_grouping(struct, shells)

        # 生成替位结构（第二原子在 xNN）
        # 询问用户生成哪一层壳：如 1、2、3 或者 1,2 或 all
        select_shell_str = ask("请选择要生成的 xNN（输入 1/2/3 或 1,2 或 all）", default="1")
        selected_shell_ids: List[int] = []
        if select_shell_str.lower() == "all":
            selected_shell_ids = [s.shell for s in shells]
        else:
            try:
                selected_shell_ids = [int(x.strip()) for x in select_shell_str.replace(" ", "").split(",") if x.strip()]
            except Exception:
                selected_shell_ids = [1]

        # 遍历选定壳层并生成结构
        for sid in selected_shell_ids:
            cand_shell = next((s for s in shells if s.shell == sid), None)
            if cand_shell is None:
                warn(f"未找到第 {sid} 壳，跳过")
                continue
            # 列出该壳候选位置
            print(f"第 {sid} 壳候选位置（共 {len(cand_shell.items)} 个）：")
            for k, it in enumerate(cand_shell.items, start=1):
                disp_idx = it.index + 1 if index_base == 1 else it.index
                print(f"  [{k:>2d}] 结构索引={disp_idx} 元素={it.element} 距离={it.distance:.4f} Å")

            gen_all = ask("是否对该壳的全部位置生成？y/n", default="y").lower().startswith("y")
            chosen_items: List[NeighborItem]
            if gen_all:
                chosen_items = list(cand_shell.items)
            else:
                pick = ask("请输入要生成的序号列表（如 1,3,5；序号非索引）", default="1").strip()
                ids: List[int] = []
                try:
                    ids = [int(x.strip()) for x in pick.split(",") if x.strip()]
                except Exception:
                    ids = [1]
                chosen_items = []
                for m in ids:
                    if 1 <= m <= len(cand_shell.items):
                        chosen_items.append(cand_shell.items[m-1])
                    else:
                        warn(f"序号 {m} 越界，已跳过")

            # 输出目录：First-Second-xnn，例如 Zr-Pu-1nn
            out_dir = os.path.join(base_output_dir, f"{first_elem}-{second_elem}-{sid}nn")
            written_files = export_substitution_variants(
                struct=struct,
                first_index=target_index,
                first_element=first_elem,
                candidates=chosen_items,
                second_element=second_elem,
                out_dir=out_dir,
                index_base=index_base,
            )
            info(f"第 {sid} 壳已生成 {len(written_files)} 个变体，目录: {out_dir}")

        # 同时保留邻居统计输出（可选）
        out = build_summary(struct, target_index, shells)
        try:
            dump_output(out, fmt=fmt, out_path=out_path)
        except Exception as e:
            err(f"写出失败: {e}")
            return 2

        if verb:
            for s in shells:
                info(f"Shell {s.shell}: count={len(s.items)} cutoff={s.cutoff:.4f} Å")
        return 0

    # ===== 非交互式：CLI 参数解析 =====
    parser = argparse.ArgumentParser(
        description="按距离识别 1NN/2NN/3NN…，支持 cutoffs 或自动分壳",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    g_in = parser.add_argument_group("输入结构")
    g_in.add_argument("--structure", type=str, default=None, help="结构文件路径（POSCAR/CONTCAR/CIF/等）")

    g_sel = parser.add_argument_group("目标选择")
    g_sel.add_argument("--target", required=True, type=str, help="目标原子：索引(按 index-base) 或 元素符号（取第一个该元素）")
    g_sel.add_argument("--index-base", type=int, choices=(0, 1), default=0, help="输入索引的基准（0或1）")

    g_shell = parser.add_argument_group("分壳设置")
    g_shell.add_argument("--cutoffs", type=str, default=None, help="以逗号分隔的壳层上限（Å），如 2.5,3.5,4.5")
    g_shell.add_argument("--shells", type=int, default=None, help="自动模式分配的壳层数")
    g_shell.add_argument("--max-radius", dest="max_radius", type=float, default=8.0, help="搜索的最大半径（Å）")
    g_shell.add_argument("--min-dist-gap", dest="min_dist_gap", type=float, default=0.15, help="自动分层时的最小距离间隙阈值（Å）")

    g_misc = parser.add_argument_group("其他")
    g_misc.add_argument("--pbc", action="store_true", help="是否使用周期性边界（当前计算默认使用 PBC）")
    g_misc.add_argument("--symmetry", action="store_true", help="是否尝试对称等价归并（当前为占位，暂不改变结果）")
    g_misc.add_argument("--output", type=str, default=None, help="输出文件路径（不指定则打印到终端）")
    g_misc.add_argument("--format", type=str, default="json", choices=("json", "yaml", "txt"), help="输出格式")
    g_misc.add_argument("--verbose", action="store_true", help="打印详细信息")

    args = parser.parse_args(argv)

    if (args.cutoffs is None) and (args.shells is None):
        # 若都未给出，默认 shells=3
        args.shells = 3
        info("未指定 cutoffs/shells，已默认 shells=3")

    if (args.cutoffs is not None) and (args.shells is not None):
        err("--cutoffs 与 --shells 不能同时使用，请二选一")
        return 2

    try:
        struct = load_structure(args.structure, None)
    except Exception as e:
        err(str(e))
        return 2

    try:
        target_index = select_target_site(struct, args.target, args.index_base)
    except Exception as e:
        err(str(e))
        return 2

    # rmax 确定：cutoffs 模式取最大 cutoff；自动模式取 max-radius
    rmax = float(args.max_radius)
    cutoffs: Optional[List[float]] = None
    if args.cutoffs is not None:
        try:
            cutoffs = parse_cutoffs(args.cutoffs)
        except Exception as e:
            err(f"解析 cutoffs 失败: {e}")
            return 2
        rmax = max(cutoffs) + 1e-6

    if args.verbose:
        info(f"目标原子 index={target_index}, 元素={struct[target_index].specie}")
        info(f"搜索半径 rmax={rmax} Å")

    # 邻居收集
    # 仅在原胞内搜索邻居（不跨边界）
    neighs = collect_neighbors(struct, target_index, rmax, use_pbc=False)

    # 分壳
    shells: List[Shell]
    if cutoffs is not None:
        shells = split_shells_by_cutoffs(neighs, cutoffs)
    else:
        assert args.shells is not None
        shells = split_shells_auto(neighs, n_shells=int(args.shells), rmax=rmax, min_gap=float(args.min_dist_gap))

    if args.symmetry:
        shells = symmetry_grouping(struct, shells)

    # 汇总与输出
    out = build_summary(struct, target_index, shells)
    try:
        dump_output(out, fmt=args.format, out_path=args.output)
    except Exception as e:
        err(f"写出失败: {e}")
        return 2

    if args.verbose:
        # 打印每壳统计
        for s in shells:
            info(f"Shell {s.shell}: count={len(s.items)} cutoff={s.cutoff:.4f} Å")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        err("用户取消")
        sys.exit(130)
