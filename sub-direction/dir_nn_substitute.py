import os
import json
import math
import sys
from datetime import datetime
from typing import List, Tuple, Optional, Dict

import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.vasp import Poscar
from pymatgen.transformations.site_transformations import ReplaceSiteSpeciesTransformation, RemoveSitesTransformation
from pymatgen.io.cif import CifParser


def info(msg: str):
    print(f"[INFO] {msg}")


def warn(msg: str):
    print(f"[WARN] {msg}")


def err(msg: str):
    print(f"[ERROR] {msg}")


def ask(prompt: str, default: Optional[str] = None) -> str:
    if default is not None:
        prompt = f"{prompt} [默认: {default}]: "
    else:
        prompt = f"{prompt}: "
    s = input(prompt).strip()
    return s if s else (default or "")


def parse_int_triplet(s: str) -> Tuple[int, int, int]:
    parts = s.replace(",", " ").split()
    if len(parts) != 3:
        raise ValueError("需要三个整数，如 '1 0 1'")
    return tuple(int(x) for x in parts)  # type: ignore


def parse_float(s: str) -> float:
    try:
        return float(s)
    except Exception:
        raise ValueError("请输入数字")

# 将包含 numpy 类型的对象递归转换为 JSON 可序列化的内建类型
def to_jsonable(obj):
    try:
        import numpy as _np  # 局部导入，避免全局污染
    except Exception:  # numpy 肯定已装，此处仅兜底
        _np = None

    if _np is not None and isinstance(obj, _np.generic):
        return obj.item()
    if isinstance(obj, (list, tuple)):
        return [to_jsonable(x) for x in obj]
    if isinstance(obj, dict):
        return {to_jsonable(k): to_jsonable(v) for k, v in obj.items()}
    return obj


def load_structure(path: str) -> Structure:
    # 若为目录，尝试自动识别常见文件名
    if os.path.isdir(path):
        candidates = [
            "POSCAR", "CONTCAR", "POSCAR.vasp", "CONTCAR.vasp",
            "structure.cif", "input.cif"
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

    # 1) 优先按 pymatgen 自动识别（依赖扩展名或特征）
    try:
        return Structure.from_file(path)
    except Exception as e:
        last_err = e

    # 2) 无扩展名时，尝试按 VASP POSCAR 解析
    try:
        poscar = Poscar.from_file(path)
        return poscar.structure
    except Exception:
        pass

    # 3) 再尝试 CIF 解析
    try:
        parser = CifParser(path)
        structs = parser.get_structures()
        if structs:
            return structs[0]
    except Exception:
        pass

    # 4) 兜底：读文本后按 POSCAR 格式尝试
    try:
        with open(path, "r", encoding="utf-8") as f:
            txt = f.read()
        return Structure.from_str(txt, fmt="poscar")
    except Exception:
        pass

    raise ValueError(
        f"无法解析结构文件: {path}。请使用标准 POSCAR/CONTCAR 或 .cif，或改为输入具体文件名。原始错误: {last_err}"
    )


def standardize_to_primitive(struct: Structure) -> Structure:
    try:
        sga = SpacegroupAnalyzer(struct, symprec=1e-3)
        prim = sga.find_primitive()
        if prim is not None:
            info("已转换为原胞 (primitive cell)")
            return prim
        else:
            warn("未找到更小的原胞，保持输入结构不变")
            return struct
    except Exception as e:
        warn(f"标准化失败，保持原结构: {e}")
        return struct


def maybe_expand_supercell(struct: Structure, scale: Optional[Tuple[int, int, int]]) -> Structure:
    if scale is None:
        return struct
    if any(n <= 0 for n in scale):
        raise ValueError("超胞放大需为正整数，如 1 1 1")
    if scale == (1, 1, 1):
        return struct
    info(f"扩展超胞: {scale}")
    return struct * scale


def list_sites(struct: Structure, wyckoffs: Optional[List[str]] = None):
    header = f"索引  元素  分数坐标 (f.)                   Wyckoff"
    print(header)
    print("-" * len(header))
    for i, site in enumerate(struct.sites):
        f = site.frac_coords
        wy = wyckoffs[i] if wyckoffs else "-"
        print(f"{i:>4d}  {str(site.specie):>2s}   [{f[0]:7.4f} {f[1]:7.4f} {f[2]:7.4f}]   {wy}")


def get_wyckoffs(struct: Structure) -> List[str]:
    try:
        sga = SpacegroupAnalyzer(struct, symprec=1e-3)
        ds = sga.get_symmetry_dataset()
        if ds and "wyckoffs" in ds:
            return list(ds["wyckoffs"])  # type: ignore
    except Exception:
        pass
    return ["-"] * len(struct)


def select_primary_site(struct: Structure) -> int:
    wy = get_wyckoffs(struct)
    print("请选择首个替位原子选择方式：\n 1) 按元素(element)  2) 按站点索引(index)  3) 按Wyckoff字母(wyckoff)")
    mode = ask("输入 1/2/3", default="1")
    list_sites(struct, wy)

    if mode == "2":
        idx = int(ask("请输入首位点索引 index"))
        if idx < 0 or idx >= len(struct):
            raise ValueError("索引超出范围")
        return idx

    if mode == "3":
        letters = sorted({w for w in wy if w != "-"})
        if not letters:
            warn("无法识别Wyckoff，改用索引选择")
            return int(ask("请输入首位点索引 index"))
        print("可选Wyckoff字母:", ", ".join(letters))
        sel = ask("输入Wyckoff字母，例如 'a'")
        cand = [i for i, w in enumerate(wy) if w == sel]
        if not cand:
            raise ValueError("该Wyckoff无站点")
        print("该Wyckoff对应站点:", cand)
        return int(ask("从上面索引中选择一个"))

    # default: by element
    elem = ask("请输入元素符号，如 Si")
    try:
        Element(elem)
    except Exception:
        raise ValueError("无效元素符号")
    cand = [i for i, s in enumerate(struct.sites) if str(s.specie) == elem]
    if not cand:
        raise ValueError(f"结构中无该元素: {elem}")
    print("候选索引:", cand)
    idx = int(ask("选择一个索引", default=str(cand[0])))
    if idx not in cand:
        raise ValueError("所选索引不在候选中")
    return idx


def frac_direction_to_cart_unit(lattice, uvw: Tuple[int, int, int]) -> np.ndarray:
    v_frac = np.array(uvw, dtype=float)
    if np.allclose(v_frac, 0):
        raise ValueError("方向 [u v w] 不能全为 0")
    v_cart = lattice.get_cartesian_coords(v_frac)
    norm = np.linalg.norm(v_cart)
    if norm < 1e-8:
        raise ValueError("方向向量长度异常")
    return v_cart / norm


def get_first_neighbor_shell(struct: Structure, idx: int, tol: float, same_cell_only: bool = True) -> Tuple[List, float]:
    # 构建足够大的半径 r，确保囊括第一壳层
    r = max(struct.lattice.abc) + 1e-6
    all_neighs = struct.get_all_neighbors(r)
    neighs = all_neighs[idx]

    # 过滤掉自身（距离零）并可选限制为同一超胞（image == [0, 0, 0]）
    def _is_same_cell(n) -> bool:
        try:
            im = tuple(map(int, n.image))
        except Exception:
            im = tuple(n.image)
        return im == (0, 0, 0)

    candidates = []
    for n in neighs:
        if n.nn_distance <= 1e-6:
            continue
        if same_cell_only and not _is_same_cell(n):
            continue
        candidates.append(n)

    if not candidates:
        raise RuntimeError("同一超胞范围内未找到任何邻居（请检查结构、超胞或适当放宽阈值）")

    d1 = min(n.nn_distance for n in candidates)
    shell = [n for n in candidates if abs(n.nn_distance - d1) <= tol]
    return shell, d1


def choose_forward_neighbor(struct: Structure, i: int, shell: List, d_hat: np.ndarray, theta_tol_deg: float,
                            auto_fallback: bool = False,
                            max_theta_fallback: float = 45.0) -> Tuple[int, Dict]:
    """
    仅保留 'either' 模式：在角度阈值内，优先选择角度最小者；若角度相同，则选择距离最小者。
    auto_fallback: 若在 theta_tol_deg 内无解，则自动将角度放宽至 max_theta_fallback。
    """
    lat = struct.lattice
    pi = struct[i]

    def select_with_mode(theta_allow):
        best_primary_score = None
        best_secondary_score = None
        best_info_local = None
        for n in shell:
            j = n.index
            pj = struct[j]
            delta_f = pj.frac_coords + np.array(n.image, dtype=float) - pi.frac_coords
            delta = lat.get_cartesian_coords(delta_f)
            dist = np.linalg.norm(delta)
            if dist < 1e-8:
                continue
            proj = float(np.dot(delta, d_hat))
            cosang = max(-1.0, min(1.0, abs(proj) / dist))
            ang = math.degrees(math.acos(cosang))
            if ang > theta_allow:
                continue
            # 新评分：先角度小，再距离小
            primary_score = -ang
            secondary_score = -dist
            is_better = False
            if best_primary_score is None:
                is_better = True
            elif primary_score > best_primary_score:
                is_better = True
            elif primary_score == best_primary_score and secondary_score > best_secondary_score:
                is_better = True
            if is_better:
                best_primary_score = primary_score
                best_secondary_score = secondary_score
                best_info_local = {
                    "j": j,
                    "proj": proj,
                    "angle_deg": ang,
                    "dist": dist,
                    "image": list(map(int, n.image)),
                    "mode": "either",
                }
        return best_info_local

    # 第一次尝试，使用给定 θ_tol
    cand = select_with_mode(theta_tol_deg)
    if cand is not None:
        return cand["j"], cand

    # 自动回退：仅扩大角度到 max_theta_fallback
    if auto_fallback:
        cand = select_with_mode(max_theta_fallback)
        if cand is not None:
            return cand["j"], cand

    raise RuntimeError("在给定角度与回退策略下仍未找到合适最近邻，请调整参数或方向")


def verify_first_neighbor_direction(struct: Structure, i: int, j: int, d_hat: np.ndarray, tol: float, theta_tol_deg: float) -> Dict:
    shell, d1 = get_first_neighbor_shell(struct, i, tol)
    lat = struct.lattice
    pi = struct[i]
    # 寻找 j 是否在第一壳层（任一像）
    ok_shell = False
    angle_ok = False
    angle_val = None
    proj_val = None
    dist_val = None
    image_val = None
    for n in shell:
        if n.index != j:
            continue
        delta_f = struct[j].frac_coords + np.array(n.image, dtype=float) - pi.frac_coords
        delta = lat.get_cartesian_coords(delta_f)
        dist = np.linalg.norm(delta)
        if dist < 1e-8:
            continue
        proj = float(np.dot(delta, d_hat))
        cosang = max(-1.0, min(1.0, abs(proj) / dist))
        ang = math.degrees(math.acos(cosang))
        angle_ok = ang <= theta_tol_deg
        ok_shell = True
        angle_val = ang
        proj_val = proj
        dist_val = dist
        image_val = list(map(int, n.image))
        break
    return {
        "in_first_shell": ok_shell,
        "angle_ok": angle_ok,
        "angle_deg": angle_val,
        "proj": proj_val,
        "dist": dist_val,
        "image": image_val,
        "d1_min": d1,
    }


# ====== POSCAR 物种顺序与写出控制 ======

def unique_symbols_in_first_appearance(struct: Structure) -> List[str]:
    seen = set()
    order = []
    for site in struct.sites:
        sym = str(site.specie)
        if sym not in seen:
            seen.add(sym)
            order.append(sym)
    return order


def reorder_structure_by_symbol_order(struct: Structure, symbol_order: List[str]) -> Structure:
    symbol_to_sites = {sym: [] for sym in symbol_order}
    others = []
    for site in struct.sites:
        sym = str(site.specie)
        if sym in symbol_to_sites:
            symbol_to_sites[sym].append(site)
        else:
            others.append(site)
    ordered_sites = []
    for sym in symbol_order:
        ordered_sites.extend(symbol_to_sites.get(sym, []))
    # 若存在未包含在 symbol_order 的物种，顺序追加在末尾
    ordered_sites.extend(others)
    return Structure.from_sites(ordered_sites)


def main():
    info("两原子定向最近邻替位 缺陷生成脚本 (pymatgen)")

    in_path = ask("请输入结构文件路径 (支持 POSCAR/CONTCAR/cif 等)")
    try:
        struct = load_structure(in_path)
    except Exception as e:
        err(str(e))
        sys.exit(1)

    # 新增：询问输入是否已是超胞
    is_input_supercell = ask("输入的结构文件是否已经是超胞？(y/n)", default="n").lower().startswith("y")

    # 原胞标准化（保持原有逻辑）
    if ask("是否转换为原胞? (y/n)", default="y").lower().startswith("y"):
        struct = standardize_to_primitive(struct)

    # 超胞扩展（若输入已是超胞，则跳过扩胞提问与处理）
    sc = None
    if not is_input_supercell:
        sc_in = ask("是否扩展超胞? 输入如 '1 1 1'，留空表示不变", default="")
        if sc_in.strip():
            try:
                sc = parse_int_triplet(sc_in)
            except Exception as e:
                err(f"超胞输入错误: {e}")
                sys.exit(1)
            struct = maybe_expand_supercell(struct, sc)
    else:
        info("检测到输入为超胞：已跳过‘扩展超胞’步骤。")

    # 选择首个替位位点
    try:
        i = select_primary_site(struct)
    except Exception as e:
        err(f"首位点选择失败: {e}")
        sys.exit(1)

    # 晶向 [u v w]
    dir_in = ask("请输入晶体方向 [u v w]，如 '1 0 1'", default="1 0 0")
    try:
        uvw = parse_int_triplet(dir_in)
        d_hat = frac_direction_to_cart_unit(struct.lattice, uvw)
    except Exception as e:
        err(f"方向输入错误: {e}")
        sys.exit(1)

    # 缺陷类型选择
    print("请选择缺陷类型：\n 1) 替位-替位 (SS)  2) 替位-空位 (SV)  3) 空位-空位 (VV)")
    mode_in = ask("输入 1/2/3", default="1")
    defect_mode = {"1": "SS", "2": "SV", "3": "VV"}.get(mode_in, "SS")

    # 根据模式采集元素输入与标签
    elem1 = None
    elem2 = None
    label1 = None
    label2 = None

    if defect_mode == "SS":
        # 替位-替位：需要两个元素
        elem1 = ask("请输入第一个替位元素 D1，如 'Al'")
        elem2 = ask("请输入第二个替位元素 D2，如 'Al' (可与D1相同)")
        try:
            Element(elem1); Element(elem2)
        except Exception:
            err("替位元素符号不合法")
            sys.exit(1)
        label1, label2 = str(elem1), str(elem2)
    elif defect_mode == "SV":
        # 替位-空位：仅第一个元素
        elem1 = ask("请输入替位元素 D1，如 'Al'")
        try:
            Element(elem1)
        except Exception:
            err("替位元素符号不合法")
            sys.exit(1)
        label1, label2 = str(elem1), "Va"
    else:  # VV
        # 空位-空位：无元素输入
        label1, label2 = "Va", "Va"

    # 参数与候选查找（循环以支持 15° 硬性角度限制的重试）
    while True:
        # 容差参数
        tol = parse_float(ask("距离容差 tol (Å)", default="0.10"))
        theta_tol = parse_float(ask("角度容差 θ_tol (度)", default="10"))

        # 自动回退
        auto_fb = ask("自动回退策略? (y/n)", default="y").lower().startswith("y")

        # 查找第一壳层 & 方向筛选（仅 either 模式）
        try:
            shell, d1min = get_first_neighbor_shell(struct, i, tol)
            j, cand_info = choose_forward_neighbor(struct, i, shell, d_hat, theta_tol,
                                                   auto_fallback=auto_fb,
                                                   max_theta_fallback=45.0)
        except Exception as e:
            err(f"筛选方向最近邻失败: {e}")
            sys.exit(1)

        # 硬性角度限制检查（15°）
        angle_deg = cand_info.get("angle_deg")
        if angle_deg is None:
            err("错误：内部错误，未获得候选的角度信息。")
            sys.exit(1)
        if angle_deg > 15.0:
            err(f"错误：选中的邻居与指定方向 [u v w] 的夹角为 {angle_deg:.2f}°，超过了 15° 的限制。")
            if ask("是否返回重新设置距离容差和角度容差？(y/n)", default="y").lower().startswith("y"):
                continue  # 重新输入 tol / theta_tol 并重试
            else:
                sys.exit(1)

        # 角度检查通过
        info(f"选中最近邻 j={j} (投影 {cand_info['proj']:.3f} Å, 角度 {angle_deg:.2f}°) 作为第二位点")
        info(f"角度检查通过：{angle_deg:.2f}° ≤ 15°，继续执行缺陷构建。")
        break

    # 执行缺陷构建（按模式分支）
    try:
        if defect_mode == "SS":
            trans = ReplaceSiteSpeciesTransformation({i: elem1, j: elem2})
            struct_defect = trans.apply_transformation(struct)
        elif defect_mode == "SV":
            struct_tmp = ReplaceSiteSpeciesTransformation({i: elem1}).apply_transformation(struct)
            struct_defect = RemoveSitesTransformation([j]).apply_transformation(struct_tmp)
        else:  # VV
            struct_defect = RemoveSitesTransformation([i, j]).apply_transformation(struct)
    except Exception as e:
        err(f"执行缺陷构建失败: {e}")
        sys.exit(1)

    # 验证（仅 either 模式）
    ver = verify_first_neighbor_direction(struct, i, j, d_hat, tol, theta_tol)
    angle_str = f"{ver['angle_deg']:.2f}" if ver['angle_deg'] is not None else "nan"
    info(f"验证: 第一壳层={ver['in_first_shell']}, 角度满足={ver['angle_ok']}, 角度={angle_str}°")

    # 输出目录与文件（重命名为 <Base>-<D1>-<D2>-[uvw]，其中空位用 Va 表示）
    base_host = struct.composition.reduced_formula  # 基体（化简式）
    folder_name = f"{base_host}-{label1}-{label2}-[{uvw[0]}{uvw[1]}{uvw[2]}]"
    out_dir_default = os.path.join(
        os.path.dirname(os.path.abspath(in_path)),
        folder_name
    )
    out_dir = ask("输出目录", default=out_dir_default)
    os.makedirs(out_dir, exist_ok=True)

    # 写 POSCAR（使用标准文件名 POSCAR）
    poscar_path = os.path.join(out_dir, "POSCAR")
    # 初始化顺序，确保在写文件异常时也可用于 summary 和打印
    host_order = unique_symbols_in_first_appearance(struct)
    defect_order_all = unique_symbols_in_first_appearance(struct_defect)
    new_syms = [s for s in defect_order_all if s not in host_order]
    desired_order = host_order + new_syms
    try:
        # 写入主 POSCAR（按 desired_order 分组，避免如 "U Pu U" 的重复顺序）
        struct_defect_sorted = reorder_structure_by_symbol_order(struct_defect, desired_order)
        Poscar(struct_defect_sorted, sort_structure=False).write_file(poscar_path)
        info(f"已保存: {poscar_path}")
    except Exception as e:
        err(f"写入POSCAR失败: {e}")

    # 写日志
    try:
        sga0 = SpacegroupAnalyzer(struct, symprec=1e-3)
        spg = sga0.get_space_group_symbol()
    except Exception:
        spg = "-"

    summary = {
        "input_file": os.path.abspath(in_path),
        "space_group": spg,
        "supercell": sc if sc is not None else [1, 1, 1],
        "mode": defect_mode,
        "primary_index": int(i),
        "secondary_index": int(j),
        "elements": {"D1": str(label1), "D2": str(label2)},
        "direction_uvw": [int(uvw[0]), int(uvw[1]), int(uvw[2])],
        "auto_fallback": auto_fb,
        "tol_A": float(tol),
        "theta_tol_deg": float(theta_tol),
        "neighbor_choice": to_jsonable(cand_info),
        "verification": to_jsonable(ver),
        "formula_before": struct.composition.formula,
        "formula_after": struct_defect.composition.formula,
        "poscar_species_order": desired_order,
    }
    log_path = os.path.join(out_dir, "summary.json")
    with open(log_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    info(f"已保存日志: {log_path}")

    print("\n====== 总结 ======")
    print(f"模式: {defect_mode}  | 首位点 i={i}, 二位点 j={j}, 方向 [u v w] = {uvw}")
    print(f"第一壳层最小距离 d1 = {d1min:.3f} Å, 选择投影 = {cand_info['proj']:.3f} Å, 角度 = {cand_info['angle_deg']:.2f}°")
    print(f"POSCAR 输出: {poscar_path}")
    print(f"POSCAR 物种顺序: {desired_order}")
    print(f"日志输出: {log_path}")


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n用户中断。")

