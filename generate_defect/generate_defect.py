#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Unified CLI for PyDefect defect generation workflows.
Single-file, modular organization with:
- LoggerManager
- ConfigManager
- DefectSetOps
- EntriesOps
- EntryOps
- InterstitialOps

Subcommands:
- defect-set (ds)
- defect-entries (de)
- make-entry (me)
- append-interstitial (ai)
- pop-interstitial (pi)
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Third-party deps
from monty.serialization import loadfn
from tqdm import tqdm
import yaml
import json

# Standalone build: no pydefect dependency required
# This function is kept for backward compatibility and now does nothing.
def _ensure_pydefect_import(pydefect_src: Optional[str] = None) -> None:
    return



# ------------- Logging -------------
import logging


class LoggerManager:
    @staticmethod
    def setup(verbose: bool = False, quiet: bool = False, log_file: Optional[str] = None) -> logging.Logger:
        level = logging.INFO
        if verbose and not quiet:
            level = logging.DEBUG
        if quiet and not verbose:
            level = logging.WARNING
        logger = logging.getLogger("generate_defect")
        logger.setLevel(logging.DEBUG)
        # Clean handlers to avoid duplication when re-running
        logger.handlers = []
        ch = logging.StreamHandler()
        ch.setLevel(level)
        ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(ch)
        if log_file:
            fh = logging.FileHandler(log_file, encoding="utf-8")
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s"))
            logger.addHandler(fh)
        return logger


# ------------- Config -------------
DEFAULT_CONFIG: Dict[str, Any] = {
    "paths": {
        "supercell_info": "supercell_info.json",
        "defect_in": "defect_in.yaml",
        "perfect_dir": "perfect",
        "log_file": None,
    },

    "defect_set": {
        "overwritten_oxi_states": {},  # e.g., {"O": -2}
        "dopants": [],                 # e.g., ["H", "Li"]
        "ele_neg_diff": 2.0,
        "keywords": [],
    },
    "entries": {
        "skip_existing": True,
    },
    "entry": {
        "charge": None,  # int or None -> auto
        "perfect_poscar": None,  # default: <perfect_dir>/POSCAR
    },
    "interstitial": {
        "base_structure": None,  # POSCAR for primitive/unit cell; default: from supercell_info if available
        "info": "",
    },
}


class ConfigError(ValueError):
    pass


class ConfigManager:
    def __init__(self, config_file: Optional[str]):
        self._file = config_file
        self.cfg = DEFAULT_CONFIG.copy()
        if config_file:
            p = Path(config_file)
            if not p.is_file():
                raise ConfigError(f"Config file not found: {config_file}")
            with p.open("r", encoding="utf-8") as f:
                loaded = yaml.safe_load(f) or {}
            self._deep_update(self.cfg, loaded)

    def merge_cli_overrides(self, overrides: Dict[str, Any]) -> None:
        self._deep_update(self.cfg, overrides)

    def validate(self) -> None:
        paths = self.cfg.get("paths", {})
        if not isinstance(paths, dict):
            raise ConfigError("paths must be a mapping")
        # Minimal validation; detailed checks are done per operation

    @staticmethod
    def _deep_update(dst: Dict[str, Any], src: Dict[str, Any]) -> None:
        for k, v in (src or {}).items():
            if isinstance(v, dict) and isinstance(dst.get(k), dict):
                ConfigManager._deep_update(dst[k], v)
            else:
                dst[k] = v


# ------------- PyDefect imports (lazy) -------------

def _imports():
    # Removed pydefect dependency in standalone build
    raise ImportError("PyDefect-dependent operations are unavailable in the standalone build.")


# ------------- Operations -------------
@dataclass
class Context:
    logger: logging.Logger
    cfg: Dict[str, Any]


# ----- Standalone helpers (no pydefect) -----
from pymatgen.core import Element

def _charge_set(ox_state: int) -> List[int]:
    if ox_state >= 0:
        charges = [i for i in range(ox_state + 1)]
        if ox_state % 2 == 1:
            charges.insert(0, -1)
    else:
        charges = [i for i in range(ox_state, 1)]
        if ox_state % 2 == 1:
            charges.append(1)
    return charges


def _oxi_state(element: str, overwritten: Dict[str, int]) -> int:
    if element in overwritten:
        return int(overwritten[element])
    try:
        el = Element(element)
        commons = el.common_oxidation_states
        if commons:
            # choose the first common oxidation state
            return int(commons[0])
    except Exception:
        pass
    return 0


def _electronegativity(element: str) -> Optional[float]:
    try:
        x = Element(element).X
        return float(x) if x is not None else None
    except Exception:
        return None


class SupercellInfoOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    def run(self, poscar_path: str, matrix: Optional[List[int]] = None,
            symprec: Optional[float] = None, angle_tol: Optional[float] = None) -> Path:
        """
        从 POSCAR 读取结构（可为原胞或常规胞），构建 supercell_info 并写入 JSON。
        - 若提供 9 个整数的 matrix，将按 3x3 变换矩阵使用；否则使用单位矩阵 [1,0,0;0,1,0;0,0,1]。
        - 使用 pymatgen 完成结构读取、超胞构建与空间群识别。
        - 路径策略（保持向后兼容）：
          1) 若 CLI 或配置中显式指定了输出路径（paths.supercell_info），优先使用该路径；
          2) 否则，默认将 supercell_info.json 写入到 POSCAR 所在目录。
        """
        from pymatgen.core import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

        paths = self.ctx.cfg["paths"]

        # 1) 基本校验 + 解析 POSCAR 路径
        p = Path(poscar_path)
        if not p.is_file():
            raise FileNotFoundError(f"未找到 POSCAR 文件：{poscar_path}")
        poscar_dir = p.parent  # POSCAR 所在目录

        # 2) 解析输出路径（实现自动路径推导 + 配置优先级）
        #    - 若 paths.supercell_info 显式给出且不等于默认名 "supercell_info.json"，则使用该路径
        #    - 否则默认写入 POSCAR 同目录的 supercell_info.json
        cfg_supercell_info: Optional[str] = paths.get("supercell_info")
        use_auto_out = (
            cfg_supercell_info is None
            or str(cfg_supercell_info).strip().lower() == "supercell_info.json"
        )
        out_path = (poscar_dir / "supercell_info.json").resolve() if use_auto_out else Path(cfg_supercell_info).resolve()

        # 3) 读取 POSCAR
        try:
            structure = Structure.from_file(p)
        except Exception as e:
            raise ValueError(f"读取 POSCAR 失败，请检查文件格式/编码：{poscar_path}\n原因：{e}")

        # 4) 处理 matrix: 将 9 元素展平向量转为 3x3
        if matrix is not None:
            if len(matrix) != 9:
                raise ValueError("--matrix 需要提供 9 个整数，按行优先构成 3x3 矩阵")
            mat3x3: List[List[int]] = [matrix[0:3], matrix[3:6], matrix[6:9]]
        else:
            mat3x3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]

        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(
                f"[dry-run] 将从 {poscar_path} 读取结构，并生成 supercell_info 到 {out_path}"
            )
            return out_path

        # 5) 对称性参数（仅用于空间群识别，不改变结构）
        sp_symprec = symprec if symprec is not None else 1e-3
        sp_angle_tol = angle_tol if angle_tol is not None else 5.0

        # 6) 构建超胞
        supercell = structure.copy()
        try:
            supercell.make_supercell(mat3x3)
        except Exception as e:
            raise ValueError(f"按矩阵生成超胞失败：{mat3x3}，原因：{e}")

        # 7) 空间群识别（基于原始结构）
        try:
            sga = SpacegroupAnalyzer(structure, symprec=sp_symprec, angle_tolerance=sp_angle_tol)
            sg_symbol = sga.get_space_group_symbol()
        except Exception:
            sg_symbol = "Unknown"

        # 8) 元素计数
        elem_counts = {str(el): int(q) for el, q in supercell.composition.get_el_amt_dict().items()}

        # 9) 组织并写出 JSON（确保目录存在）
        out = {
            "structure": supercell.as_dict(),
            "space_group": sg_symbol,
            "transformation_matrix": mat3x3,
            "elements": elem_counts,
            "interstitials": []
        }
        try:
            out_path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            # 目录创建失败（权限/路径非法等）
            raise PermissionError(f"创建输出目录失败：{out_path.parent}\n原因：{e}")
        try:
            with open(out_path, "w", encoding="utf-8") as f:
                json.dump(out, f, ensure_ascii=False, indent=2)
        except PermissionError as e:
            raise PermissionError(f"写入 supercell_info.json 失败，可能无写入权限：{out_path}\n原因：{e}")
        except OSError as e:
            raise OSError(f"写入 supercell_info.json 失败：{out_path}\n原因：{e}")
        self.ctx.logger.info(f"已写出 supercell_info: {out_path}")
        return out_path


class DefectSetOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    def run(self) -> Path:
        paths = self.ctx.cfg["paths"]
        supercell_info_path = Path(paths["supercell_info"]).resolve()
        if not supercell_info_path.is_file():
            raise FileNotFoundError(f"Missing {supercell_info_path}")
        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] Would load supercell_info from {supercell_info_path}")
        # 读取 supercell_info.json（本地自定义结构）
        with open(supercell_info_path, "r", encoding="utf-8") as f:
            sc_info = json.load(f)
        ds_cfg = self.ctx.cfg["defect_set"]
        overwritten = ds_cfg.get("overwritten_oxi_states") or {}
        dopants = ds_cfg.get("dopants") or []
        ele_neg_diff = float(ds_cfg.get("ele_neg_diff", 2.0))
        keywords = ds_cfg.get("keywords") or []
        force_charge = ds_cfg.get("force_charge", None)

        host_elements = list((sc_info.get("elements") or {}).keys())

        # 生成缺陷集合（简化版）：空位 + 掺杂取代 (+ 可选间隙)
        defect_map: Dict[str, List[int]] = {}

        # 空位：Va_Element
        for host in host_elements:
            charges = _charge_set(-_oxi_state(host, overwritten))
            name = f"Va_{host}"
            defect_map[name] = charges

        # 取代：dopant_host（限制电负性差）
        for inp in dopants:
            xin = _electronegativity(inp)
            ox_in = _oxi_state(inp, overwritten)
            for host in host_elements:
                if inp == host:
                    continue
                xh = _electronegativity(host)
                if xin is None or xh is None or abs(xin - xh) > ele_neg_diff:
                    continue
                ox_host = _oxi_state(host, overwritten)
                charges = _charge_set(ox_in - ox_host)
                name = f"{inp}_{host}"
                defect_map[name] = charges

        # 间隙：inp_iN（如 supercell_info 指定 interstitials，按 dopants 生成）
        inters = sc_info.get("interstitials") or []
        if inters and dopants:
            for i, _ in enumerate(inters, start=1):
                for inp in dopants:
                    charges = _charge_set(_oxi_state(inp, overwritten))
                    name = f"{inp}_i{i}"
                    defect_map[name] = charges

        # 关键词筛选（若提供）
        if keywords:
            defect_map = {k: v for k, v in defect_map.items() if any(kw in k for kw in keywords)}

        # 覆盖电荷（若提供 --q）
        if force_charge is not None:
            try:
                fc = int(force_charge)
            except Exception:
                fc = 0
            defect_map = {k: [fc] for k in defect_map.keys()}

        # 写出 YAML
        out_path = Path(paths["defect_in"]).resolve()
        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] 将写出缺陷集合到 {out_path}")
        else:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            with open(out_path, "w", encoding="utf-8") as f:
                yaml.safe_dump(defect_map, f, allow_unicode=True, sort_keys=True)
            self.ctx.logger.info(f"Wrote defect set: {out_path}")
        return out_path


class EntriesOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    def _load_supercell_structure(self, supercell_info_path: Path):
        from pymatgen.core import Structure
        with open(supercell_info_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        if "structure" not in data:
            raise ValueError("supercell_info.json missing 'structure'")
        structure = Structure.from_dict(data["structure"])
        interstitials = data.get("interstitials") or []
        return structure, interstitials

    @staticmethod
    def _write_prior_info_yaml(dir_path: Path, charge: int) -> None:
        content = f"charge: {int(charge)}\n"
        (dir_path / "prior_info.yaml").write_text(content, encoding="utf-8")

    @staticmethod
    def _site_symmetry_symbol(structure) -> str:
        try:
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            sga = SpacegroupAnalyzer(structure, symprec=1e-3, angle_tolerance=5.0)
            return sga.get_point_group_symbol()
        except Exception:
            return "N/A"

    @staticmethod
    def _defect_entry_json(name: str, charge: int, structure) -> Dict[str, Any]:
        # Choose a heuristic defect center: use barycenter of species removed/added later if available.
        # Here we fall back to geometric center (0,0,0) if unknown; callers should override if known.
        defect_center = [0.0, 0.0, 0.0]
        try:
            # Use lattice-centered origin as default; many downstream tools only need a valid list length.
            defect_center = [0.0, 0.0, 0.0]
        except Exception:
            pass
        site_sym = EntriesOps._site_symmetry_symbol(structure)
        return {
            "@module": "pydefect.input_maker.defect_entry",
            "@class": "DefectEntry",
            "@version": "0.2.4",
            "name": name,
            "charge": int(charge),
            "structure": structure.as_dict(),
            "perturbed_structure": None,
            "site_symmetry": site_sym,
            "defect_center": defect_center,
        }

    @staticmethod
    def _choose_first_site_index(structure, element_symbol: str) -> Optional[int]:
        for i, site in enumerate(structure.sites):
            if str(site.specie) == element_symbol:
                return i
        return None

    @staticmethod
    def _parse_defect_name(raw_name: str) -> Tuple[str, Dict[str, Any]]:
        """Return (kind, info) where kind in {vacancy, substitution, interstitial}.
        info contains keys depending on kind.
        Supported patterns:
          - Vacancy:  Va_X or Va_X1 (index optional)
          - Subst.:   A_on_B or A_B or A_on_B1 (index optional)
          - Interst.: A_iN (N = 1-based interstitial index)
        """
        name = raw_name.strip()
        # Interstitial: A_iN
        if "_i" in name:
            try:
                sp, tail = name.split("_i", 1)
                idx = int(tail)
                return "interstitial", {"dopant": sp, "index": idx, "base_name": f"{sp}_i{idx}"}
            except Exception:
                pass
        # Vacancy: Va_X or Va_X1
        if name.startswith("Va_"):
            body = name[3:]
            head = ''.join([c for c in body if c.isalpha()])
            num = ''.join([c for c in body if c.isdigit()])
            site_index = int(num) if num else 1
            return "vacancy", {"element": head, "site_index": site_index, "base_name": f"Va_{head}{site_index}"}
        # Substitution: A_on_B or A_B
        if "_on_" in name:
            a, b = name.split("_on_")
            return "substitution", {"dopant": a, "host": b, "site_index": 1, "base_name": f"{a}_on_{b}1"}
        if "_" in name:
            a, b = name.split("_", 1)
            return "substitution", {"dopant": a, "host": b, "site_index": 1, "base_name": f"{a}_on_{b}1"}
        # Fallback: treat as interstitial label without index (unsupported)
        return "unknown", {"raw": name, "base_name": name}

    @staticmethod
    def _apply_defect(structure, kind: str, info: Dict[str, Any], interstitials: List[Dict[str, Any]]):
        from pymatgen.core import Element
        defect_center_frac: Optional[List[float]] = None
        s = structure.copy()
        if kind == "vacancy":
            elem = info["element"]
            # choose index-th occurrence of elem (1-based)
            occ = 0
            target_idx: Optional[int] = None
            for i, site in enumerate(s.sites):
                if str(site.specie) == elem:
                    occ += 1
                    if occ == info.get("site_index", 1):
                        target_idx = i
                        break
            if target_idx is None:
                raise ValueError(f"No site found for vacancy {elem}{info.get('site_index', 1)}")
            defect_center_frac = s[target_idx].frac_coords.tolist()
            # remove site
            s.remove_sites([target_idx])
        elif kind == "substitution":
            host = info["host"]
            dopant = info["dopant"]
            occ = 0
            target_idx = None
            for i, site in enumerate(s.sites):
                if str(site.specie) == host:
                    occ += 1
                    if occ == info.get("site_index", 1):
                        target_idx = i
                        break
            if target_idx is None:
                raise ValueError(f"No site found for substitution {dopant}_on_{host}{info.get('site_index', 1)}")
            defect_center_frac = s[target_idx].frac_coords.tolist()
            s[target_idx] = Element(dopant)
        elif kind == "interstitial":
            idx = info["index"]
            if idx <= 0 or idx > len(interstitials):
                raise IndexError(f"Interstitial index out of range: {idx}")
            fc = interstitials[idx - 1].get("frac_coords")
            if not fc or len(fc) != 3:
                raise ValueError(f"Invalid interstitial frac_coords at index {idx}")
            dopant = info.get("dopant")
            if dopant is None:
                raise ValueError("Interstitial dopant element is required")
            s.append(Element(dopant), fc)
            defect_center_frac = [float(fc[0]), float(fc[1]), float(fc[2])]
        else:
            raise ValueError(f"Unsupported defect kind for {info}")
        return s, defect_center_frac

    def run(self) -> None:
        paths = self.ctx.cfg["paths"]
        supercell_info_path = Path(paths["supercell_info"]).resolve()
        defect_in_path = Path(paths["defect_in"]).resolve()
        perfect_dir = Path(paths["perfect_dir"]).resolve()
        if not supercell_info_path.is_file():
            raise FileNotFoundError(f"Missing {supercell_info_path}")
        if not defect_in_path.is_file():
            raise FileNotFoundError(f"Missing {defect_in_path}; run defect-set first.")

        # Load supercell and interstitial hints
        structure, interstitials = self._load_supercell_structure(supercell_info_path)

        # Ensure perfect dir and POSCAR
        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] Would ensure {perfect_dir} and write POSCAR")
        else:
            perfect_dir.mkdir(exist_ok=True)
            structure.to(filename=str(perfect_dir / "POSCAR"))

        # Load defect set mapping
        with open(defect_in_path, "r", encoding="utf-8") as f:
            defect_map: Dict[str, List[int]] = yaml.safe_load(f) or {}

        created, skipped = [], []
        from tqdm import tqdm as _tqdm
        for name, charges in _tqdm(sorted(defect_map.items(), key=lambda kv: kv[0]), desc="entries"):
            kind, info = self._parse_defect_name(name)
            if kind == "unknown":
                self.ctx.logger.warning(f"Skip unknown defect pattern: {name}")
                continue
            base_name = info["base_name"]
            # Build defect structure once per name, then reuse per charge
            try:
                defect_structure, defect_center = self._apply_defect(structure, kind, info, interstitials)
            except Exception as e:
                self.ctx.logger.warning(f"Skip {name}: {e}")
                continue

            # site symmetry from the unrelaxed defect structure
            site_sym = self._site_symmetry_symbol(defect_structure)

            for q in sorted(set(int(c) for c in (charges or []))):
                dir_name = f"{base_name}_{q}"
                dir_path = Path(dir_name)
                if dir_path.exists() and self.ctx.cfg["entries"].get("skip_existing", True):
                    skipped.append(dir_name)
                    continue

                if self.ctx.cfg.get("dry_run"):
                    self.ctx.logger.info(f"[dry-run] Would create {dir_path}/POSCAR + prior_info.yaml + defect_entry.json")
                    created.append(dir_name)
                    continue

                dir_path.mkdir(exist_ok=True)
                # POSCAR
                defect_structure.to(filename=str(dir_path / "POSCAR"))
                # prior_info.yaml
                self._write_prior_info_yaml(dir_path, q)
                # defect_entry.json
                entry_dict = {
                    "@module": "pydefect.input_maker.defect_entry",
                    "@class": "DefectEntry",
                    "@version": "0.2.4",
                    "name": base_name.rsplit("_", 1)[0],  # without charge suffix
                    "charge": int(q),
                    "structure": defect_structure.as_dict(),
                    "perturbed_structure": None,
                    "site_symmetry": site_sym,
                    "defect_center": defect_center if defect_center is not None else [0.0, 0.0, 0.0],
                }
                with open(dir_path / "defect_entry.json", "w", encoding="utf-8") as jf:
                    json.dump(entry_dict, jf, ensure_ascii=False, indent=2)
                created.append(dir_name)

        self.ctx.logger.info(f"Created: {len(created)}, Skipped: {len(skipped)}")
        if skipped:
            self.ctx.logger.debug(f"Skipped entries: {skipped}")


class EntryOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    @staticmethod
    def _infer_defect_center(perfect_structure, defect_structure) -> List[float]:
        import numpy as np
        lat = perfect_structure.lattice
        # Map defect sites to perfect sites and find unmatched/mismatched species
        p_coords = perfect_structure.frac_coords
        d_coords = defect_structure.frac_coords
        if len(d_coords) > len(p_coords):
            # Interstitial: find extra site by nearest match threshold
            used = set()
            centers = []
            for j, d in enumerate(d_coords):
                # nearest perfect distance
                dm = np.array([lat.get_distance_and_image(d, p)[0] for p in p_coords])
                if dm.min() > 0.3:  # Angstrom
                    centers.append(d.tolist())
            if centers:
                return [float(x) for x in centers[0]]
        elif len(d_coords) < len(p_coords):
            # Vacancy: find missing perfect site
            used = set()
            centers = []
            for i, p in enumerate(p_coords):
                dm = np.array([lat.get_distance_and_image(p, d)[0] for d in d_coords])
                if dm.min() > 0.3:
                    centers.append(p.tolist())
            if centers:
                return [float(x) for x in centers[0]]
        # Substitution or fallback: choose site with largest species mismatch
        try:
            from pymatgen.analysis.structure_matcher import StructureMatcher
            sm = StructureMatcher(primitive_cell=False, attempt_supercell=False)
            mapping = sm.get_mapping(perfect_structure, defect_structure)
            if mapping:
                # pick first mapped index
                idx = list(mapping.keys())[0]
                return [float(x) for x in perfect_structure[idx].frac_coords]
        except Exception:
            pass
        return [0.0, 0.0, 0.0]

    def run(self, dir_: str, name: str, charge: Optional[int] = None, perfect_poscar: Optional[str] = None) -> Path:
        from pymatgen.core import Structure
        d = Path(dir_).resolve()
        if not (d / "POSCAR").is_file():
            raise FileNotFoundError(f"POSCAR not found in {d}")
        if perfect_poscar is None:
            perfect_dir = Path(self.ctx.cfg["paths"]["perfect_dir"]).resolve()
            perfect_poscar = str(perfect_dir / "POSCAR")
        if not Path(perfect_poscar).is_file():
            raise FileNotFoundError(f"Perfect POSCAR not found: {perfect_poscar}")

        if charge is None:
            # Standalone build cannot infer charge from INCAR/POTCAR; default to 0
            if self.ctx.cfg.get("dry_run"):
                self.ctx.logger.info(f"[dry-run] Would estimate charge state from INCAR/POTCAR in {d}")
                charge = 0
            self.ctx.logger.warning("charge not specified; defaulting to 0 in standalone mode")

        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] Would make defect_entry.json in {d} for name={name}, q={charge}")
            return d / "defect_entry.json"

        defect_structure = Structure.from_file(d / "POSCAR")
        perfect_structure = Structure.from_file(perfect_poscar)

        # site symmetry and defect center
        try:
            from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
            site_sym = SpacegroupAnalyzer(defect_structure, symprec=1e-3, angle_tolerance=5.0).get_point_group_symbol()
        except Exception:
            site_sym = "N/A"
        defect_center = self._infer_defect_center(perfect_structure, defect_structure)

        entry_dict = {
            "@module": "pydefect.input_maker.defect_entry",
            "@class": "DefectEntry",
            "@version": "0.2.4",
            "name": name,
            "charge": int(charge),
            "structure": defect_structure.as_dict(),
            "perturbed_structure": None,
            "site_symmetry": site_sym,
            "defect_center": defect_center,
        }
        out = d / "defect_entry.json"
        with open(out, "w", encoding="utf-8") as f:
            json.dump(entry_dict, f, ensure_ascii=False, indent=2)
        # Also ensure prior_info.yaml exists/updated
        EntriesOps._write_prior_info_yaml(d, int(charge))
        self.ctx.logger.info(f"Wrote {out}")
        return out


class InterstitialOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    def append(self, frac_coords: List[float], info: str = "") -> Path:
        # Standalone: directly edit JSON structure; no pydefect
        paths = self.ctx.cfg["paths"]
        supercell_info_path = Path(paths["supercell_info"]).resolve()
        if not supercell_info_path.is_file():
            raise FileNotFoundError(f"Missing {supercell_info_path}")
        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] Would append interstitial at {frac_coords} with info '{info}'")
            return supercell_info_path
        # Load JSON
        with open(supercell_info_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        inters = data.get("interstitials")
        if inters is None or not isinstance(inters, list):
            inters = []
        inters.append({"frac_coords": [float(frac_coords[0]), float(frac_coords[1]), float(frac_coords[2])],
                       "info": info})
        data["interstitials"] = inters
        with open(supercell_info_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        self.ctx.logger.info(f"Appended interstitial. Updated {supercell_info_path}")
        return supercell_info_path

    def pop(self, index: Optional[int], all_: bool = False, assume_yes: bool = False) -> Path:
        paths = self.ctx.cfg["paths"]
        supercell_info_path = Path(paths["supercell_info"]).resolve()
        if not supercell_info_path.is_file():
            raise FileNotFoundError(f"Missing {supercell_info_path}")
        if self.ctx.cfg.get("dry_run"):
            act = "pop all" if all_ else f"pop index {index}"
            self.ctx.logger.info(f"[dry-run] Would {act} interstitial(s)")
            return supercell_info_path
        # Confirm
        if not assume_yes:
            msg = "Remove ALL interstitials? [y/N]" if all_ else f"Remove interstitial #{index}? [y/N]"
            try:
                ans = input(msg + " ").strip().lower()
            except EOFError:
                ans = "n"
            if ans not in ("y", "yes"):
                self.ctx.logger.info("Cancelled")
                return supercell_info_path
        # Load JSON and modify
        with open(supercell_info_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        inters = data.get("interstitials")
        if inters is None or not isinstance(inters, list):
            inters = []
        if all_:
            inters = []
        else:
            if not index or index <= 0:
                raise ValueError("index must be positive (1-based)")
            try:
                inters.pop(index - 1)
            except IndexError:
                raise IndexError(f"No interstitial at index {index}")
        data["interstitials"] = inters
        with open(supercell_info_path, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        self.ctx.logger.info(f"Updated {supercell_info_path}")
        return supercell_info_path


class LocalExtremaOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    def _load_volumetric(self, files: List[str]):
        from pymatgen.io.vasp.outputs import Chgcar, Locpot
        vol = None
        for fp in files:
            p = Path(fp)
            if not p.is_file():
                raise FileNotFoundError(f"Volumetric data not found: {fp}")
            name = p.name.upper()
            obj = None
            try:
                if ("AECCAR" in name) or ("CHGCAR" in name):
                    obj = Chgcar.from_file(str(p))
                elif "LOCPOT" in name:
                    obj = Locpot.from_file(str(p))
                else:
                    # try CHGCAR first, then LOCPOT
                    try:
                        obj = Chgcar.from_file(str(p))
                    except Exception:
                        obj = Locpot.from_file(str(p))
            except Exception as e:
                raise ValueError(f"Failed to read volumetric data from {fp}: {e}")
            if vol is None:
                vol = obj
            else:
                vol = vol + obj
        return vol

    @staticmethod
    def _grid_and_shape(vol) -> Tuple[Any, Tuple[int, int, int]]:
        grid = vol.data.get("total")
        if grid is None:
            # Some versions use 'data' directly
            grid = getattr(vol, "data", None)
        if grid is None:
            raise ValueError("Volumetric data does not contain 'total' grid")
        return grid, grid.shape

    @staticmethod
    def _find_local_extrema(grid, find_max: bool) -> List[Tuple[int, int, int]]:
        import numpy as np
        try:
            from skimage.feature import peak_local_max
        except Exception:
            # Fallback: 6-neighbor check with PBC
            mask = np.ones(grid.shape, dtype=bool)
            for axis in range(3):
                rolled_p = np.roll(grid, 1, axis=axis)
                rolled_n = np.roll(grid, -1, axis=axis)
                if find_max:
                    mask &= (grid >= rolled_p) & (grid >= rolled_n)
                else:
                    mask &= (grid <= rolled_p) & (grid <= rolled_n)
            idxs = np.argwhere(mask)
            return [tuple(idx) for idx in idxs]
        # Use skimage with PBC via 3x3x3 tiling
        sign = 1.0 if find_max else -1.0
        tiled = np.tile(sign * grid, reps=(3, 3, 3))
        coords = peak_local_max(tiled, min_distance=1)
        # Map back to central cell [1,2) in fractional of the tiled box
        fcoords = [coord / np.array(tiled.shape) * 3.0 for coord in coords]
        fcoords = [f - 1.0 for f in fcoords if np.all(np.array(f) >= 1.0) and np.all(np.array(f) < 2.0)]
        # Convert fractional in central cell to original indices
        nx, ny, nz = grid.shape
        idxs: List[Tuple[int, int, int]] = []
        for fa, fb, fc in fcoords:
            ia = int(np.clip(np.floor(fa * nx), 0, nx - 1))
            ib = int(np.clip(np.floor(fb * ny), 0, ny - 1))
            ic = int(np.clip(np.floor(fc * nz), 0, nz - 1))
            idxs.append((ia, ib, ic))
        return idxs

    @staticmethod
    def _frac_from_index(idx: Tuple[int, int, int], shape: Tuple[int, int, int]) -> List[float]:
        nx, ny, nz = shape
        i, j, k = idx
        # approximate center of voxel
        return [float((i + 0.5) / nx), float((j + 0.5) / ny), float((k + 0.5) / nz)]

    @staticmethod
    def _min_distance_to_atoms(structure, fcoord: List[float]) -> float:
        import numpy as np
        from pymatgen.core import Structure
        if not isinstance(structure, Structure):
            return float("inf")
        f1 = np.array([fcoord])
        f2 = structure.frac_coords
        dmat = structure.lattice.get_all_distances(f1, f2)
        return float(dmat.min()) if dmat.size else float("inf")

    def _cluster_nodes_scipy(self, structure, fcoords: List[List[float]], tol: float) -> List[List[float]]:
        # Try SciPy hierarchical clustering with PBC; fallback to greedy if unavailable
        if not fcoords:
            return []
        try:
            import numpy as np
            from scipy.cluster.hierarchy import linkage, fcluster
            from scipy.spatial.distance import squareform
        except Exception:
            # Fallback: greedy cluster using first-come
            cands = [{"frac_coords": fc, "value": 0.0} for fc in fcoords]
            merged = self._greedy_cluster(structure, cands, tol)
            return [m["frac_coords"] for m in merged]
        import numpy as np
        lat = structure.lattice
        # Distance matrix with PBC
        dist = np.array(lat.get_all_distances(fcoords, fcoords), dtype=float)
        dist = (dist + dist.T) / 2.0
        np.fill_diagonal(dist, 0.0)
        Z = linkage(squareform(dist), method="single")
        cn = fcluster(Z, tol, criterion="distance")
        # Merge by averaging with PBC images to the first point
        merged: List[List[float]] = []
        for n in set(cn):
            idxs = np.where(cn == n)[0].tolist()
            anchor = fcoords[idxs[0]]
            acc = []
            for i in idxs:
                f = np.array(fcoords[i])
                d, image = lat.get_distance_and_image(anchor, f)
                acc.append(f + image)
            ave = np.mean(np.array(acc), axis=0)
            ave = ave - np.floor(ave)
            ave = ave * (np.abs(ave - 1.0) > 1e-15)
            merged.append(ave.tolist())
        return merged

    @staticmethod
    def _greedy_cluster(structure, candidates: List[Dict[str, Any]], tol: float) -> List[Dict[str, Any]]:
        kept: List[Dict[str, Any]] = []
        for cand in candidates:
            keep = True
            for kc in kept:
                d = structure.lattice.get_distance_and_image(kc["frac_coords"], cand["frac_coords"])[0]
                if d < tol:
                    keep = False
                    break
            if keep:
                kept.append(cand)
        return kept

    def run(self, files: List[str], find_max: bool = False,
            threshold_frac: Optional[float] = None,
            threshold_abs: Optional[float] = None,
            min_dist: float = 0.5,
            tol: float = 0.5,
            radius: float = 0.4,
            info: str = "") -> Path:
        import numpy as np
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.core import Element
        vol = self._load_volumetric(files)
        grid, shape = self._grid_and_shape(vol)
        structure = vol.structure
        nx, ny, nz = shape
        ngrid = nx * ny * nz

        self.ctx.logger.info(f"Grid shape: {shape}; find_max={find_max}")
        if ngrid > 2_000_000:
            self.ctx.logger.info("Large grid detected; operations may be slow. Consider coarser grids if needed.")

        # 1) raw extrema (PBC-aware)
        idxs = self._find_local_extrema(grid, find_max=find_max)
        self.ctx.logger.info(f"Raw extrema points: {len(idxs)}")

        # 2) build candidates
        cands: List[Dict[str, Any]] = []
        for idx in idxs:
            val = float(grid[idx])
            fc = self._frac_from_index(idx, shape)
            cands.append({"idx": idx, "frac_coords": fc, "value": val})

        # 3) threshold filters
        if threshold_abs is not None:
            if find_max:
                cands = [c for c in cands if c["value"] >= float(threshold_abs)]
            else:
                cands = [c for c in cands if c["value"] <= float(threshold_abs)]
        if threshold_frac is not None and 0.0 <= threshold_frac <= 1.0 and cands:
            cands.sort(key=lambda x: x["value"], reverse=find_max)
            n_keep = int(len(cands) * float(threshold_frac))
            cands = cands[:n_keep] if n_keep > 0 else []
        self.ctx.logger.info(f"After thresholding: {len(cands)} candidates")

        # 4) remove collisions near atoms
        if min_dist and min_dist > 0:
            before = len(cands)
            cands = [c for c in cands if self._min_distance_to_atoms(structure, c["frac_coords"]) > float(min_dist)]
            self.ctx.logger.info(f"After remove_collisions(min_dist={min_dist}): {len(cands)} (removed {before - len(cands)})")

        # 5) clustering by tol
        if tol and tol > 0 and cands:
            fcoords = [c["frac_coords"] for c in cands]
            merged_fcoords = self._cluster_nodes_scipy(structure, fcoords, float(tol))
            # Re-sample values at nearest grid voxel centers
            recands: List[Dict[str, Any]] = []
            for fc in merged_fcoords:
                ia = int(np.clip(np.floor(fc[0] * nx), 0, nx - 1))
                ib = int(np.clip(np.floor(fc[1] * ny), 0, ny - 1))
                ic = int(np.clip(np.floor(fc[2] * nz), 0, nz - 1))
                recands.append({"idx": (ia, ib, ic), "frac_coords": [float(fc[0]), float(fc[1]), float(fc[2])], "value": float(grid[ia, ib, ic])})
            cands = recands
            self.ctx.logger.info(f"After clustering(tol={tol}): {len(cands)} representatives")

        # 6) local average within radius
        def ave_value_for_fc(fc: List[float]) -> float:
            if radius is None or radius <= 0:
                return float("nan")
            aa = np.linspace(0, 1, nx, endpoint=False)
            bb = np.linspace(0, 1, ny, endpoint=False)
            cc = np.linspace(0, 1, nz, endpoint=False)
            AA, BB, CC = np.meshgrid(aa, bb, cc, indexing="ij")
            pts = np.vstack([AA.ravel(), BB.ravel(), CC.ravel()]).T
            d = structure.lattice.get_all_distances(pts, np.array(fc)).reshape(nx, ny, nz)
            mask = d < float(radius)
            vol_sphere = structure.volume * (mask.sum() / float(ngrid))
            if mask.sum() == 0 or vol_sphere == 0:
                return float(grid[int(fc[0]*nx)%nx, int(fc[1]*ny)%ny, int(fc[2]*nz)%nz])
            chg = float(np.sum(grid * mask) / mask.size / vol_sphere)
            return chg

        for c in cands:
            c["ave_value"] = ave_value_for_fc(c["frac_coords"]) if radius and radius > 0 else c["value"]

        # 7) symmetry grouping and coordination
        added = structure.copy()
        start = len(added)
        for c in cands:
            added.append(Element("Og"), c["frac_coords"])  # dummy species
        end = len(added)
        extrema_points = []
        try:
            sga = SpacegroupAnalyzer(added, symprec=1e-3)
            ds = sga.get_symmetry_dataset()
            eq = ds.get("equivalent_atoms")
            syms = ds.get("site_symmetry_symbols")
            groups: Dict[int, List[int]] = {}
            for idx in range(start, end):
                label = int(eq[idx]) if eq is not None else idx
                groups.setdefault(label, []).append(idx)

            def compute_coordination(fc: List[float]) -> Dict[str, Any]:
                dists = structure.lattice.get_all_distances([fc], structure.frac_coords)[0]
                mind = float(np.min(dists)) if len(dists) else 0.0
                cutoff = 1.3 * mind
                dist_dict: Dict[str, List[float]] = {}
                neigh_idx: List[int] = []
                for i, site in enumerate(structure):
                    dist = float(dists[i])
                    if dist < cutoff and dist > 1e-5:
                        elem = str(site.specie)
                        dist_dict.setdefault(elem, []).append(round(dist, 2))
                        neigh_idx.append(i)
                for k in list(dist_dict.keys()):
                    dist_dict[k] = [float(x) for x in sorted(dist_dict[k])]
                return {"distance_dict": dist_dict, "cutoff": round(cutoff, 3), "neighboring_atom_indices": neigh_idx}

            for label, idx_list in groups.items():
                idx_list = sorted(idx_list)
                repr_idx = min(idx_list)
                # map appended indices to candidate list order
                coords: List[List[float]] = []
                quants: List[float] = []
                for i in idx_list:
                    j = i - start
                    fc = cands[j]["frac_coords"]
                    coords.append(fc)
                    # choose ave_value if finite else value
                    av = cands[j].get("ave_value", float("nan"))
                    quants.append(float(av if av == av else cands[j]["value"]))
                site_sym = (syms[repr_idx] if syms is not None else "N/A").replace(".", "")
                coord = compute_coordination(coords[0])
                extrema_points.append({"site_symmetry": site_sym, "coordination": coord, "frac_coords": coords, "quantities": quants})
        except Exception as e:
            self.ctx.logger.warning(f"Symmetry grouping failed ({e}); falling back to no-group output")
            for c in cands:
                coord = {"distance_dict": {}, "cutoff": 0.0, "neighboring_atom_indices": []}
                extrema_points.append({"site_symmetry": "N/A", "coordination": coord, "frac_coords": [c["frac_coords"]], "quantities": [c["ave_value"]]})

        out = {
            "unit_cell": structure.as_dict(),
            "is_min": (not find_max),
            "extrema_points": extrema_points,
            "info": info or "",
            "params": {
                "threshold_frac": threshold_frac,
                "threshold_abs": threshold_abs,
                "min_dist": min_dist,
                "tol": tol,
                "radius": radius,
            }
        }

        # 确定输出目录：默认写入与输入体积数据文件相同的目录（取第一个输入文件的父目录）
        try:
            input_paths = [Path(fp).resolve() for fp in files]
            parent_dirs = {p.parent for p in input_paths}
            out_dir = next(iter(parent_dirs))
            if len(parent_dirs) > 1:
                self.ctx.logger.warning(f"检测到多个输入文件来自不同目录，输出将保存到首个目录：{out_dir}")
        except Exception as e:
            # 回退策略：若输入路径解析失败，则退回当前工作目录
            out_dir = Path.cwd()
            self.ctx.logger.warning(f"解析输入路径失败，退回当前目录：{out_dir}（原因：{e}）")

        out_path = (out_dir / "volumetric_data_local_extrema.json").resolve()
        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] 将写出本地极值 JSON 到：{out_path}（点数：{len(extrema_points)}）")
            return out_path
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as f:
            json.dump(out, f, ensure_ascii=False, indent=2)
        self.ctx.logger.info(f"已写出本地极值 JSON：{out_path}")
        return out_path


class AddInterstitialsOps:
    def __init__(self, ctx: Context):
        self.ctx = ctx

    def run(self, local_extrema_path: str, indices: Optional[List[int]] = None) -> Path:
        lep = Path(local_extrema_path)
        if not lep.is_file():
            raise FileNotFoundError(f"Missing {local_extrema_path}")
        try:
            with open(lep, "r", encoding="utf-8") as f:
                le = json.load(f)
        except Exception as e:
            # 友好提示：用户可能把 LOCPOT/AECCAR/CHGCAR 当作 --local-extrema 传入
            name = lep.name.upper()
            if name == "LOCPOT" or "AECCAR" in name or name == "CHGCAR":
                raise ValueError(
                    "参数 --local-extrema 需要传入 local-extrema 步骤生成的 JSON 文件（例如 volumetric_data_local_extrema.json）。\n"
                    f"检测到传入的是 {lep.name}，这看起来是体积数据文件。请先运行：\n"
                    "  python generate_defect/generate_defect.py local-extrema --volumetric-data <LOCPOT/AECCAR0 AECCAR2>\n"
                    "然后将生成的 JSON 路径传给 add-interstitials 的 --local-extrema 参数。"
                ) from e
            raise
        points = le.get("extrema_points", [])
        if not points:
            raise ValueError("No extrema_points in local extrema JSON")
        # Determine indices (1-based)
        if not indices:
            indices = list(range(1, len(points) + 1))
        # Load supercell_info
        sc_path = Path(self.ctx.cfg["paths"]["supercell_info"]).resolve()
        if not sc_path.is_file():
            raise FileNotFoundError(f"Missing {sc_path}")
        with open(sc_path, "r", encoding="utf-8") as f:
            sc = json.load(f)
        inters = sc.get("interstitials")
        if inters is None or not isinstance(inters, list):
            inters = []
        base_info = le.get("info", "")
        # Transform unit-cell fractional coords to supercell fractional coords via inv(transformation_matrix)
        import numpy as np
        trans_mat = np.array(sc.get("transformation_matrix") or [[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        try:
            inv_mat = np.linalg.inv(trans_mat)
        except Exception as e:
            raise ValueError(f"Invalid transformation_matrix in supercell_info: {trans_mat}") from e
        for idx in indices:
            if idx <= 0 or idx > len(points):
                raise IndexError(f"Index out of range: {idx}")
            ep = points[idx - 1]
            fc_u = ep.get("frac_coords", [[0, 0, 0]])[0]
            if not isinstance(fc_u, list) or len(fc_u) != 3:
                raise ValueError(f"Invalid frac_coords for extrema index {idx}")
            # Right-multiply row vector by inv_mat: (a_s, b_s, c_s) = (a_u, b_u, c_u) . inv(T)
            fc_s = (np.array(fc_u, dtype=float) @ inv_mat).tolist()
            # Wrap into [0,1)
            fc_s = [float(x % 1.0) for x in fc_s]
            info = (base_info + f" #{idx}").strip()
            inters.append({"frac_coords": [float(fc_s[0]), float(fc_s[1]), float(fc_s[2])], "info": info})
        sc["interstitials"] = inters
        if self.ctx.cfg.get("dry_run"):
            self.ctx.logger.info(f"[dry-run] Would append {len(indices)} interstitial(s) to {sc_path}")
            return sc_path
        with open(sc_path, "w", encoding="utf-8") as f:
            json.dump(sc, f, ensure_ascii=False, indent=2)
        self.ctx.logger.info(f"Updated supercell_info: {sc_path}")
        return sc_path




# ------------- Template generation -------------
class TemplateOps:
    def __init__(self, logger: logging.Logger):
        self.logger = logger

    @staticmethod
    def _build_template_yaml() -> str:
        # 中文注释丰富的 YAML 模板，覆盖所有可配置项与示例值
        return (
            "# generate_defect 配置模板 (YAML)\n"
            "# 说明：\n"
            "# - 配置优先级：CLI 参数 > 此配置文件 > 脚本内置默认\n"
            "# - Windows 路径建议使用双引号并避免空格；如需绝对路径请自行替换\n"
            "# - 带有示例值与中文注释，修改为你的项目实际值后再使用\n"
            "\n"
            "paths:\n"
            "  # 超胞信息文件 (由本脚本 supercell-info 步骤生成)\n"
            "  supercell_info: \"supercell_info.json\"\n"
            "  # 缺陷集合输出文件 (供 defect-entries 使用)\n"
            "  defect_in: \"defect_in.yaml\"\n"
            "  # 存放完美结构 POSCAR 的目录名称\n"
            "  perfect_dir: \"perfect\"\n"
            "  # 日志文件；设为 null 则只输出到控制台\n"
            "  log_file: \"generate_defect.log\"\n"
            "\n"
            "defect_set:\n"
            "  # 手动覆盖价态 (示例：氧化物)；键为元素，值为整数价态\n"
            "  overwritten_oxi_states: {O: -2, Zn: 2}  # 示例：ZnO\n"
            "  # 掺杂元素列表 (示例：Si 中 B 掺杂、P 掺杂)\n"
            "  dopants: [B, P]  # 可留空 []\n"
            "  # 允许的电负性差阈值 (浮点数)\n"
            "  ele_neg_diff: 2.0\n"
            "  # 关键字筛选，常见有 Va(空位), i(间隙), X_on_Y(置换) 等\n"
            "  keywords: [Va, i]  # 示例：只生成空位和间隙\n"
            "\n"
            "entries:\n"
            "  # 批量生成时是否跳过已存在目录\n"
            "  skip_existing: true\n"
            "\n"
            "entry:\n"
            "  # 单条目 make-entry 的电荷；null 表示自动从 INCAR/POTCAR 估算\n"
            "  charge: null\n"
            "  # 完美结构 POSCAR；为 null 时使用 <paths.perfect_dir>/POSCAR\n"
            "  perfect_poscar: null\n"
            "\n"
            "interstitial:\n"
            "  # 原胞/标准胞的 POSCAR 路径；若 supercell_info 中已包含 unitcell_structure 可置 null\n"
            "  base_structure: null\n"
            "  # 追加的间隙位点备注信息\n"
            "  info: \"\"\n"
        )

    def run(self, output: Optional[str], assume_yes: bool = False) -> Path:
        out = Path(output or "generate_defect_template.yaml").resolve()
        if out.exists() and not assume_yes:
            try:
                ans = input(f"目标文件已存在：{out}\n是否覆盖? [y/N] ").strip().lower()
            except EOFError:
                ans = "n"
            if ans not in ("y", "yes"):
                self.logger.info("已取消生成模板。")
                return out
        out.parent.mkdir(parents=True, exist_ok=True)
        content = self._build_template_yaml()
        # 使用 UTF-8 BOM 以便在 Windows 记事本/部分终端中正确显示中文
        with out.open("w", encoding="utf-8-sig") as f:
            f.write(content)
        self.logger.info(f"已生成配置模板：{out}")
        return out

# ------------- CLI -------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="generate_defect", description="Standalone defect-generation CLI (no PyDefect dependency)")
    p.add_argument("--config", type=str, help="YAML config file")

    p.add_argument("--dry-run", action="store_true", help="Preview actions without writing")
    p.add_argument("--yes", "-y", action="store_true", help="Assume yes to prompts")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    p.add_argument("--quiet", "-q", action="store_true", help="Quiet logging (warnings only)")
    p.add_argument("--log-file", type=str, help="Log file path")

    sp = p.add_subparsers(dest="cmd", required=True)

    # template/init-config
    p_tpl = sp.add_parser("template", aliases=["init-config"], help="Generate a commented YAML config template")
    p_tpl.add_argument("--output", type=str, default="generate_defect_template.yaml", help="Output template filename")

    # supercell-info
    p_si = sp.add_parser("supercell-info", aliases=["si"], help="Generate supercell_info.json from POSCAR (primitive or conventional)")
    p_si.add_argument("--poscar", required=True, type=str, help="Path to POSCAR (primitive or conventional cell)")
    p_si.add_argument("--supercell-info", type=str, help="Output path for supercell_info.json")
    p_si.add_argument("--matrix", type=int, nargs=9, help="3x3 integer transformation matrix flattened (row-major), 9 ints")
    p_si.add_argument("--symprec", type=float, help="Symmetry length tolerance (spglib)")
    p_si.add_argument("--angle-tol", type=float, help="Symmetry angle tolerance (degrees)")

    # defect-set
    p_ds = sp.add_parser("defect-set", aliases=["ds"], help="Generate defect_in.yaml from supercell_info.json")
    p_ds.add_argument("--supercell-info", type=str, help="Path to supercell_info.json")
    p_ds.add_argument("--defect-in", type=str, help="Output defect_in.yaml path")
    p_ds.add_argument("--overwritten-oxi", type=str, help="JSON mapping for oxi states, e.g. '{\"O\": -2}'")
    p_ds.add_argument("--dopants", type=str, nargs="*", help="List of dopant elements")
    p_ds.add_argument("--ele-neg-diff", type=float, help="Electronegativity difference threshold")
    p_ds.add_argument("--keywords", type=str, nargs="*", help="Keyword filters")
    p_ds.add_argument("--q", type=int, help="Force all defects to a single charge (e.g., --q 0)")

    # defect-entries
    p_de = sp.add_parser("defect-entries", aliases=["de"], help="Batch create defect directories and files")
    p_de.add_argument("--supercell-info", type=str)
    p_de.add_argument("--defect-in", type=str)
    p_de.add_argument("--perfect-dir", type=str)
    p_de.add_argument("--no-skip-existing", action="store_true", help="Do not skip existing directories")

    # make-entry
    p_me = sp.add_parser("make-entry", aliases=["me"], help="Create defect_entry.json for a single directory")
    p_me.add_argument("--dir", required=True, type=str, help="Defect folder containing POSCAR")
    p_me.add_argument("--name", required=True, type=str, help="Defect name (e.g., Va_O1)")
    p_me.add_argument("--charge", type=int, help="Charge state (override); default: auto from INCAR/POTCAR")
    p_me.add_argument("--perfect-poscar", type=str, help="Path to perfect POSCAR; default: <perfect_dir>/POSCAR")

    # append-interstitial
    p_ai = sp.add_parser("append-interstitial", aliases=["ai"], help="Append an interstitial site to supercell_info.json")
    p_ai.add_argument("--supercell-info", type=str)
    p_ai.add_argument("--base-structure", type=str, help="POSCAR for primitive/unit cell; if omitted, use unitcell in supercell_info")
    p_ai.add_argument("--frac-coords", type=float, nargs=3, required=True, help="Fractional coords (a b c)")
    p_ai.add_argument("--info", type=str, default="", help="Optional label/info")

    # pop-interstitial
    p_pi = sp.add_parser("pop-interstitial", aliases=["pi"], help="Remove interstitial site(s) from supercell_info.json")
    p_pi.add_argument("--supercell-info", type=str)
    g = p_pi.add_mutually_exclusive_group(required=True)
    g.add_argument("--index", type=int, help="1-based index to remove")
    g.add_argument("--all", action="store_true", dest="all_", help="Remove all interstitials")

    # local-extrema (new)
    p_le = sp.add_parser("local-extrema", aliases=["le"], help="Find local extrema from volumetric data (AECCAR/CHGCAR/LOCPOT)")
    p_le.add_argument("--volumetric-data", type=str, nargs="+", required=True, help="Paths to volumetric data files (e.g., AECCAR0 AECCAR2 or LOCPOT)")
    p_le.add_argument("--find-max", action="store_true", help="Search local maxima instead of minima")
    p_le.add_argument("--threshold-frac", type=float, default=None, help="Fraction of extrema to keep (0-1)")
    p_le.add_argument("--threshold-abs", type=float, default=None, help="Absolute threshold on value (min: <=, max: >=)")
    p_le.add_argument("--min-dist", type=float, default=0.5, help="Minimum distance to host atoms in Angstrom")
    p_le.add_argument("--tol", type=float, default=0.5, help="Clustering tolerance in Angstrom")
    p_le.add_argument("--radius", type=float, default=0.4, help="Spherical radius for local averaging (Angstrom)")
    p_le.add_argument("--info", type=str, help="Info string to record in JSON")

    # add-interstitials (new)
    p_ais = sp.add_parser("add-interstitials", aliases=["ais"], help="Append interstitials from local extrema JSON to supercell_info.json")
    p_ais.add_argument("--supercell-info", type=str)
    p_ais.add_argument("--local-extrema", type=str, required=True, help="Path to volumetric_data_local_extrema.json")
    p_ais.add_argument("--indices", type=int, nargs="*", help="1-based indices to add; default: all")

    return p


def _apply_cli_overrides(cfgm: ConfigManager, args: argparse.Namespace) -> None:
    overrides: Dict[str, Any] = {"paths": {}, "defect_set": {}, "entries": {}, "entry": {}, "interstitial": {}}
    if args.log_file:
        overrides["paths"]["log_file"] = args.log_file
    if getattr(args, "supercell_info", None):
        overrides["paths"]["supercell_info"] = args.supercell_info
    if getattr(args, "defect_in", None):
        overrides["paths"]["defect_in"] = args.defect_in
    if getattr(args, "perfect_dir", None):
        overrides["paths"]["perfect_dir"] = args.perfect_dir
    if getattr(args, "no_skip_existing", False):
        overrides["entries"]["skip_existing"] = False
    if getattr(args, "overwritten_oxi", None):
        try:
            overrides["defect_set"]["overwritten_oxi_states"] = json.loads(args.overwritten_oxi)
        except Exception as e:
            raise ConfigError(f"--overwritten-oxi must be JSON: {e}")
    if getattr(args, "dopants", None) is not None:
        overrides["defect_set"]["dopants"] = args.dopants
    if getattr(args, "ele_neg_diff", None) is not None:
        overrides["defect_set"]["ele_neg_diff"] = args.ele_neg_diff
    if getattr(args, "keywords", None) is not None:
        overrides["defect_set"]["keywords"] = args.keywords
    if getattr(args, "q", None) is not None:
        overrides["defect_set"]["force_charge"] = int(args.q)
    if getattr(args, "perfect_poscar", None):
        overrides["entry"]["perfect_poscar"] = args.perfect_poscar
    if getattr(args, "charge", None) is not None:
        overrides["entry"]["charge"] = args.charge
    if getattr(args, "base_structure", None):
        overrides["interstitial"]["base_structure"] = args.base_structure
    if getattr(args, "info", None) is not None:
        overrides["interstitial"]["info"] = args.info
    overrides["dry_run"] = bool(getattr(args, "dry_run", False))
    cfgm.merge_cli_overrides(overrides)


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Prepare config
    cfgm = ConfigManager(args.config)

    _apply_cli_overrides(cfgm, args)
    cfgm.validate()

    # Imports


    # Logger
    logger = LoggerManager.setup(verbose=args.verbose, quiet=args.quiet, log_file=cfgm.cfg["paths"].get("log_file"))
    ctx = Context(logger=logger, cfg=cfgm.cfg)

    try:
        if args.cmd in ("supercell-info", "si"):
            # 可选：覆盖输出路径
            if getattr(args, "supercell_info", None):
                _apply_cli_overrides(cfgm, args)
            out = SupercellInfoOps(ctx).run(
                poscar_path=args.poscar,
                matrix=getattr(args, "matrix", None),
                symprec=getattr(args, "symprec", None),
                angle_tol=getattr(args, "angle_tol", None),
            )
            logger.info(f"supercell_info.json 已生成: {out}")
        elif args.cmd in ("defect-set", "ds"):
            out = DefectSetOps(ctx).run()
            logger.info(f"defect_in.yaml ready at: {out}")
        elif args.cmd in ("defect-entries", "de"):
            EntriesOps(ctx).run()
        elif args.cmd in ("make-entry", "me"):
            charge = cfgm.cfg["entry"].get("charge")
            perfect_poscar = cfgm.cfg["entry"].get("perfect_poscar")
            out = EntryOps(ctx).run(args.dir, args.name, charge=charge, perfect_poscar=perfect_poscar)
            logger.info(f"defect_entry.json ready at: {out}")
        elif args.cmd in ("append-interstitial", "ai"):
            frac = args.frac_coords
            info = cfgm.cfg["interstitial"].get("info", "")
            path = InterstitialOps(ctx).append(frac, info)
            logger.info(f"Updated supercell_info: {path}")
        elif args.cmd in ("pop-interstitial", "pi"):
            path = InterstitialOps(ctx).pop(index=args.index, all_=args.all_, assume_yes=args.yes)
            logger.info(f"Updated supercell_info: {path}")
        elif args.cmd in ("local-extrema", "le"):
            out = LocalExtremaOps(ctx).run(
                files=args.volumetric_data,
                find_max=bool(getattr(args, "find_max", False)),
                threshold_frac=getattr(args, "threshold_frac", None),
                threshold_abs=getattr(args, "threshold_abs", None),
                min_dist=args.min_dist,
                tol=args.tol,
                radius=args.radius,
                info=(args.info if getattr(args, "info", None) is not None else ""),
            )
            logger.info(f"local extrema JSON ready at: {out}")
        elif args.cmd in ("add-interstitials", "ais"):
            if getattr(args, "supercell_info", None):
                _apply_cli_overrides(cfgm, args)
            path = AddInterstitialsOps(ctx).run(args.local_extrema, indices=getattr(args, "indices", None))
            logger.info(f"Updated supercell_info: {path}")
        elif args.cmd in ("template", "init-config"):
            out = TemplateOps(logger).run(output=getattr(args, "output", None), assume_yes=args.yes)
            logger.info(f"Template ready at: {out}")
        else:
            parser.print_help()
            return 2
        return 0
    except Exception as e:
        logger.error(f"Error: {e}")
        if args.verbose:
            logger.exception(e)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

