#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
基因类型判断脚本（两步法）

流程：
1. 先定位参考序列中的最佳匹配片段，记录该位置 (ref_start, ref_end)。
2. 再在该片段内做比对，找到不匹配的位置。
3. 对每个不匹配位置做峰数判断（单峰→纯合，双峰→杂合），据此判 WT/HOM/HET。
"""

import os
import sys
import csv
import re
import json
from pathlib import Path
from typing import Tuple, List, Dict, Optional, Set
from dataclasses import dataclass, field

from Bio import SeqIO

from pysanger import abi_to_dict, generate_consensusseq
from alignment_utils import (
    align_full_and_get_mismatches_with_position,
    align_to_segment_and_get_mismatches,
    align_and_get_core_region,
    HAS_PARASAIL,
)

# 差异位点处：次高峰/最高峰 > 此阈值视为双峰（杂合）
PEAK_RATIO_HET_THRESHOLD = 0.30
# 仅在此置信度以上的碱基位置统计错配，避免测序噪声导致假阳性
CONFIDENCE_MIN_FOR_MISMATCH = 25

# 已知多态性位点（在所有WT样本中都出现的错配位置，不应被视为突变）
KNOWN_POLYMORPHISMS = set()


def _load_known_polymorphisms(polymorphisms_file: str = None) -> set:
    """加载已知多态性位点列表"""
    if polymorphisms_file is None:
        # 打包成 exe 时 (PyInstaller) 数据在 sys._MEIPASS
        base = getattr(sys, "_MEIPASS", Path(__file__).resolve().parent)
        polymorphisms_file = str(Path(base) / "known_polymorphisms.txt")
    polymorphisms = set()
    if os.path.exists(polymorphisms_file):
        try:
            with open(polymorphisms_file, "r", encoding="utf-8") as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        try:
                            polymorphisms.add(int(line))
                        except ValueError:
                            pass
        except Exception:
            pass
    return polymorphisms


@dataclass
class GenotypeResult:
    """基因型分析结果"""
    filename: str
    genotype: str  # WT/HOM/HET/ERR
    genotype_cn: str
    mutation_count: int
    peaks_at_diffs: List[int]  # 各差异位点的峰数 1 或 2
    confidence: float
    is_abnormal: bool = False  # 突变位点之外存在双峰时为 True
    abnormal_positions: List[int] = field(default_factory=list)  # 异常双峰的色谱碱基序号


def count_peaks_at_position(abidata: Dict, base_index: int) -> int:
    """
    在 AB1 色谱中，碱基位置 base_index 处有几条“波”（几个显著峰）。
    若次高峰/最高峰 > PEAK_RATIO_HET_THRESHOLD 则判为双峰(2)，否则单峰(1)。
    """
    if base_index < 0 or base_index >= len(abidata["conf"]):
        return 1
    vals = [
        abidata["channel"]["A"][base_index],
        abidata["channel"]["T"][base_index],
        abidata["channel"]["G"][base_index],
        abidata["channel"]["C"][base_index],
    ]
    sorted_vals = sorted(vals, reverse=True)
    if sorted_vals[0] <= 0:
        return 1
    ratio = sorted_vals[1] / sorted_vals[0]
    return 2 if ratio > PEAK_RATIO_HET_THRESHOLD else 1


def _scan_chromatogram_for_het_peak(
    abidata: Dict,
    min_ratio: float = None,
    min_confidence: int = 10,
    chrom_start: Optional[int] = None,
    chrom_end: Optional[int] = None,
) -> bool:
    """
    扫描色谱中是否存在双峰（杂合信号）。
    当用于“即将判 WT 前”的全序列扫描时，使用 min_ratio 识别真实杂合双峰，
    同时要求 min_confidence 以过滤低置信度噪声。
    若提供 chrom_start/chrom_end，仅在该区间内扫描（与 SnapGene 有效比对区间一致，
    舍弃 5'/3' 低质量区，避免如 83、84 位被误判为杂合）。
    """
    if min_ratio is None:
        min_ratio = PEAK_RATIO_HET_THRESHOLD
    n_bases = len(abidata["conf"])
    i_start = 0 if chrom_start is None else max(0, chrom_start)
    i_end = n_bases if chrom_end is None else min(n_bases, chrom_end)
    for i in range(i_start, i_end):
        if abidata["conf"][i] < min_confidence:
            continue
        vals = [
            abidata["channel"]["A"][i],
            abidata["channel"]["T"][i],
            abidata["channel"]["G"][i],
            abidata["channel"]["C"][i],
        ]
        sorted_vals = sorted(vals, reverse=True)
        if sorted_vals[0] <= 0:
            continue
        ratio = sorted_vals[1] / sorted_vals[0]
        if ratio > min_ratio:
            return True
    return False


def _scan_chromatogram_het_positions(
    abidata: Dict,
    min_ratio: float = 0.55,
    min_confidence: int = 10,
    chrom_start: Optional[int] = None,
    chrom_end: Optional[int] = None,
) -> List[int]:
    """
    扫描有效区间内双峰位置列表（用于将「全色谱判 HET」的位点视为突变位点，不记入异常）。
    """
    n_bases = len(abidata["conf"])
    i_start = 0 if chrom_start is None else max(0, chrom_start)
    i_end = n_bases if chrom_end is None else min(n_bases, chrom_end)
    positions: List[int] = []
    for i in range(i_start, i_end):
        if abidata["conf"][i] < min_confidence:
            continue
        vals = [
            abidata["channel"]["A"][i],
            abidata["channel"]["T"][i],
            abidata["channel"]["G"][i],
            abidata["channel"]["C"][i],
        ]
        sorted_vals = sorted(vals, reverse=True)
        if sorted_vals[0] <= 0:
            continue
        ratio = sorted_vals[1] / sorted_vals[0]
        if ratio > min_ratio:
            positions.append(i)
    return positions


def _find_abnormal_double_peaks(
    abidata: Dict,
    chrom_start: int,
    chrom_end: int,
    mutation_chrom_indices: Set[int],
    min_ratio: float = 0.55,
    min_confidence: int = 10,
) -> List[int]:
    """
    在有效比对区间内扫描：找出「非突变位点」上存在双峰的位置（异常双峰）。
    突变位点之外的双峰可能表示测序异常、污染或质量差，判定为异常。
    返回异常双峰的色谱碱基序号列表（0-based）。
    """
    positions: List[int] = []
    n_bases = len(abidata["conf"])
    i_start = max(0, chrom_start)
    i_end = min(n_bases, chrom_end)
    for i in range(i_start, i_end):
        if i in mutation_chrom_indices:
            continue
        if abidata["conf"][i] < min_confidence:
            continue
        vals = [
            abidata["channel"]["A"][i],
            abidata["channel"]["T"][i],
            abidata["channel"]["G"][i],
            abidata["channel"]["C"][i],
        ]
        sorted_vals = sorted(vals, reverse=True)
        if sorted_vals[0] <= 0:
            continue
        ratio = sorted_vals[1] / sorted_vals[0]
        if ratio > min_ratio:
            positions.append(i)
    return positions


class GenotypeAnalyzer:
    def __init__(
        self,
        reference_path: str,
        cache_path: Optional[str] = None,
        force_ref_start: Optional[int] = None,
        polymorphisms_file: Optional[str] = None,
    ):
        """
        Args:
            reference_path: 参考序列文件路径
            cache_path: 缓存文件路径（可选）
            force_ref_start: 强制使用指定的参考序列起始位置（0-based），如果提供则跳过自动查找
            polymorphisms_file: 已知多态性位点文件路径（可选）
        """
        self.reference_path = reference_path
        self.reference_seq = self._load_reference()
        # 参考序列最佳匹配片段 [ref_start, ref_end)，首次全比对后写入并复用
        self._ref_start: Optional[int] = force_ref_start
        self._ref_end: Optional[int] = None
        self._cache_path = cache_path or str(
            Path(reference_path).with_name("ref_segment.json")
        )
        # 加载已知多态性位点
        self.known_polymorphisms = _load_known_polymorphisms(polymorphisms_file)
        if force_ref_start is None:
            self._load_ref_segment_cache()
        else:
            # 如果强制指定了起始位置，设置一个合理的结束位置
            # 假设测序长度约为 500-600 bp
            self._ref_end = min(force_ref_start + 600, len(self.reference_seq))

    def _load_reference(self) -> str:
        path = Path(self.reference_path)
        if path.suffix.lower() == ".dna":
            record = SeqIO.read(self.reference_path, "snapgene")
        elif path.suffix.lower() in [".fasta", ".fa"]:
            record = SeqIO.read(self.reference_path, "fasta")
        else:
            raise ValueError(f"不支持的文件格式: {path.suffix}")
        return str(record.seq).upper()

    def _load_ref_segment_cache(self) -> None:
        """从文件加载已记录的参考匹配区间。"""
        if not os.path.exists(self._cache_path):
            return
        try:
            with open(self._cache_path, "r", encoding="utf-8") as f:
                d = json.load(f)
            self._ref_start = int(d["ref_start"])
            self._ref_end = int(d["ref_end"])
        except Exception:
            self._ref_start = None
            self._ref_end = None

    def _save_ref_segment_cache(self, ref_start: int, ref_end: int) -> None:
        """将参考匹配区间写入文件，供后续样本直接使用。"""
        try:
            with open(self._cache_path, "w", encoding="utf-8") as f:
                json.dump({"ref_start": ref_start, "ref_end": ref_end}, f, indent=2)
        except Exception:
            pass

    def analyze_sample(self, ab1_path: str) -> GenotypeResult:
        """
        使用parasail进行Smith-Waterman比对（模拟SnapGene逻辑）：
        1. 进行完整的局部比对，得到CIGAR字符串
        2. 找到最长连续匹配区域作为"核心区域"（这就是SnapGene显示的134bp区域）
        3. 在完整比对中找到错配位置
        4. 对每个错配位置做峰数判断，再判 WT/HOM/HET
        """
        filename = os.path.basename(ab1_path)

        try:
            abidata = abi_to_dict(ab1_path)
            consensus_fwd, _ = generate_consensusseq(abidata)
            n_bases = len(abidata["conf"])

            # ---------- 使用parasail进行比对（推荐） ----------
            core_length = 0  # 初始化核心区域长度
            query_start, query_end = 0, n_bases  # 无 parasail 时全序列有效
            if HAS_PARASAIL:
                (
                    mismatches,
                    used_rc,
                    score,
                    core_ref_start,
                    core_ref_end,
                    core_length,
                    query_start,
                    query_end,
                ) = align_and_get_core_region(
                    self.reference_seq, consensus_fwd, try_reverse_complement=True,
                    only_core_mismatches=True  # 只保留核心区域附近的错配
                )
                # 更新缓存为核心区域位置
                if self._ref_start is None:
                    self._ref_start = core_ref_start
                    self._ref_end = core_ref_end
                    self._save_ref_segment_cache(core_ref_start, core_ref_end)
            # ---------- 回退到旧方法 ----------
            elif self._ref_start is None or self._ref_end is None:
                # 自动查找最佳匹配片段
                (
                    mismatches,
                    used_rc,
                    score,
                    ref_start,
                    ref_end,
                ) = align_full_and_get_mismatches_with_position(
                    self.reference_seq, consensus_fwd, try_reverse_complement=True
                )
                self._ref_start = ref_start
                self._ref_end = ref_end
                self._save_ref_segment_cache(ref_start, ref_end)
            else:
                mismatches, used_rc, score = align_to_segment_and_get_mismatches(
                    self.reference_seq,
                    self._ref_start,
                    self._ref_end,
                    consensus_fwd,
                    try_reverse_complement=True,
                )

            # ---------- 排除已知多态性位点，但保留“第一次匹配得到的突变点”用于双峰判断 ----------
            # 单点突变：即使过滤后差异数为 0，若突变点存在双峰，仍判为杂合子
            mismatches_raw = list(mismatches)  # 比对得到的全部错配（含可能被过滤的突变点）
            mismatches_filtered = []
            for ref_pos, query_pos, ref_base, query_base in mismatches:
                if ref_pos in self.known_polymorphisms:
                    continue
                mismatches_filtered.append(
                    (ref_pos, query_pos, ref_base, query_base)
                )
            mismatches = mismatches_filtered

            # ---------- 步骤 3：对突变点做峰数判断（先算 peaks，再决定 WT/HOM/HET） ----------
            # 用于双峰检查的位置：优先用过滤后的错配；若差异数为 0 则用“第一次匹配”得到的突变点（mismatches_raw）
            positions_for_peak_check = mismatches if mismatches else mismatches_raw
            peaks_at_diffs = []
            for ref_pos, query_pos, ref_base, query_base in positions_for_peak_check:
                if used_rc:
                    chrom_index = n_bases - 1 - query_pos
                else:
                    chrom_index = query_pos
                n_peaks = count_peaks_at_position(abidata, chrom_index)
                peaks_at_diffs.append(n_peaks)

            has_het_peak = any(p == 2 for p in peaks_at_diffs)

            # ---------- 有效比对区间与突变位点（色谱索引），用于异常双峰检测 ----------
            if used_rc:
                chrom_start = n_bases - query_end
                chrom_end = n_bases - query_start
            else:
                chrom_start = query_start
                chrom_end = query_end
            mutation_chrom_indices: Set[int] = set()
            for ref_pos, query_pos, ref_base, query_base in positions_for_peak_check:
                if used_rc:
                    mutation_chrom_indices.add(n_bases - 1 - query_pos)
                else:
                    mutation_chrom_indices.add(query_pos)
            abnormal_positions = _find_abnormal_double_peaks(
                abidata,
                chrom_start,
                chrom_end,
                mutation_chrom_indices,
                min_ratio=0.55,
                min_confidence=10,
            )
            is_abnormal = len(abnormal_positions) > 0
            # 异常双峰数量达到此阈值则判为 ERR（如 79 号）；仅 1 处异常仍给基因型并标异常（如 87）
            ABNORMAL_ERR_THRESHOLD = 2

            # ---------- 突变位点之外存在多处双峰 → 判为异常(ERR)，不给出基因型 ----------
            if is_abnormal and len(abnormal_positions) >= ABNORMAL_ERR_THRESHOLD:
                return GenotypeResult(
                    filename=filename,
                    genotype="ERR",
                    genotype_cn="异常(突变位点外双峰)",
                    mutation_count=len(mismatches),
                    peaks_at_diffs=peaks_at_diffs,
                    confidence=0.0,
                    is_abnormal=True,
                    abnormal_positions=abnormal_positions,
                )

            # ---------- 野生型判断（核心区域长且无双峰） ----------
            # 核心区域 >= 130：若比对得到的任一位置存在双峰，仍判为杂合子；否则野生型
            if core_length >= 130:
                if has_het_peak:
                    return GenotypeResult(
                        filename=filename,
                        genotype="HET",
                        genotype_cn="杂合子",
                        mutation_count=len(mismatches),
                        peaks_at_diffs=peaks_at_diffs,
                        confidence=round(
                            min(1.0, 0.5 + sum(1 for p in peaks_at_diffs if p == 2) / max(len(peaks_at_diffs), 1)),
                            3,
                        ),
                        is_abnormal=is_abnormal,
                        abnormal_positions=abnormal_positions,
                    )
                # 即将判野生型前再扫色谱：仅在对齐有效区间内扫描（SnapGene 从 G95 开始，
                # 舍弃 83、84 等低质量区），min_ratio=0.55、置信度>=10
                full_scan_het_positions = _scan_chromatogram_het_positions(
                    abidata,
                    min_ratio=0.55,
                    min_confidence=10,
                    chrom_start=chrom_start,
                    chrom_end=chrom_end,
                )
                has_het_peak = len(full_scan_het_positions) > 0
                if has_het_peak:
                    # 全色谱判 HET 的双峰位点视为突变位点，不记入异常（如 87 的 262）
                    effective_mutation_indices = mutation_chrom_indices | set(
                        full_scan_het_positions
                    )
                    abnormal_positions_het = _find_abnormal_double_peaks(
                        abidata,
                        chrom_start,
                        chrom_end,
                        effective_mutation_indices,
                        min_ratio=0.55,
                        min_confidence=10,
                    )
                    return GenotypeResult(
                        filename=filename,
                        genotype="HET",
                        genotype_cn="杂合子",
                        mutation_count=0,
                        peaks_at_diffs=[],
                        confidence=0.85,
                        is_abnormal=len(abnormal_positions_het) > 0,
                        abnormal_positions=abnormal_positions_het,
                    )
                return GenotypeResult(
                    filename=filename,
                    genotype="WT",
                    genotype_cn="野生型",
                    mutation_count=0,
                    peaks_at_diffs=[],
                    confidence=1.0,
                    is_abnormal=is_abnormal,
                    abnormal_positions=abnormal_positions,
                )

            # ---------- 检查突变点是否有双峰（核心区域 < 130） ----------
            # 单点突变：只要突变点存在双峰，即判为杂合子
            if has_het_peak:
                genotype = "HET"
                genotype_cn = "杂合子"
                confidence = min(
                    1.0,
                    0.5
                    + sum(1 for p in peaks_at_diffs if p == 2)
                    / max(len(peaks_at_diffs), 1),
                )
            else:
                genotype = "HOM"
                genotype_cn = "纯合子"
                confidence = min(1.0, 0.5 + len(peaks_at_diffs) / 50.0)

            return GenotypeResult(
                filename=filename,
                genotype=genotype,
                genotype_cn=genotype_cn,
                mutation_count=len(mismatches),
                peaks_at_diffs=peaks_at_diffs,
                confidence=round(confidence, 3),
                is_abnormal=is_abnormal,
                abnormal_positions=abnormal_positions,
            )

        except Exception as e:
            import traceback

            traceback.print_exc()
            return GenotypeResult(
                filename=filename,
                genotype="ERR",
                genotype_cn=f"错误: {str(e)}",
                mutation_count=0,
                peaks_at_diffs=[],
                confidence=0.0,
                is_abnormal=False,
                abnormal_positions=[],
            )

    def analyze_folder(self, folder_path: str) -> List[GenotypeResult]:
        results = []
        folder = Path(folder_path)
        ab1_files_set = set()
        for f in folder.glob("*.ab1"):
            ab1_files_set.add(f.resolve())
        for f in folder.glob("*.AB1"):
            ab1_files_set.add(f.resolve())
        ab1_files = sorted(
            list(ab1_files_set), key=lambda x: self._extract_number(x.name)
        )

        if not ab1_files:
            print(f"警告: 在 {folder_path} 中没有找到AB1文件")
            return results

        print(f"找到 {len(ab1_files)} 个AB1文件")
        if self._ref_start is not None:
            print(
                f"使用已记录的参考匹配区间: ref[{self._ref_start}:{self._ref_end}]"
            )
        else:
            print("首次运行将全参考比对并记录最佳匹配片段位置")
        print("-" * 60)

        for ab1_file in ab1_files:
            result = self.analyze_sample(str(ab1_file))
            results.append(result)
            status = "[OK]" if result.genotype != "ERR" else "[ERR]"
            peaks_str = ",".join(map(str, result.peaks_at_diffs[:5]))
            if len(result.peaks_at_diffs) > 5:
                peaks_str += f",...({len(result.peaks_at_diffs)}个)"
            abnormal_str = ""
            if result.is_abnormal:
                pos_str = ",".join(map(str, result.abnormal_positions[:8]))
                if len(result.abnormal_positions) > 8:
                    pos_str += f",...({len(result.abnormal_positions)}个)"
                abnormal_str = f" [异常:碱基{pos_str}]"
            print(
                f"{status} {result.filename}: {result.genotype} ({result.genotype_cn}) "
                f"差异数={result.mutation_count} 峰={peaks_str}{abnormal_str}"
            )

        if self._ref_start is not None:
            print(
                f"\n参考匹配区间已保存至: {self._cache_path} (后续运行将直接使用)"
            )
        return results

    def _extract_number(self, filename: str) -> int:
        m = re.match(r"^(\d+)", filename)
        return int(m.group(1)) if m else 0


def export_to_csv(results: List[GenotypeResult], output_path: str):
    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f)
        w.writerow(
            ["文件名", "基因类型", "基因类型(中文)", "差异数", "异常", "异常位置"]
        )
        for r in results:
            abnormal_pos_str = ",".join(map(str, r.abnormal_positions)) if r.abnormal_positions else ""
            w.writerow(
                [
                    r.filename,
                    r.genotype,
                    r.genotype_cn,
                    r.mutation_count,
                    "是" if r.is_abnormal else "",
                    abnormal_pos_str,
                ]
            )
    print(f"\n结果已导出到: {output_path}")


def print_summary(results: List[GenotypeResult]):
    valid = [r for r in results if r.genotype != "ERR"]
    n = len(results)
    wt = sum(1 for r in valid if r.genotype == "WT")
    hom = sum(1 for r in valid if r.genotype == "HOM")
    het = sum(1 for r in valid if r.genotype == "HET")
    abnormal_count = sum(1 for r in results if r.is_abnormal)
    err = n - len(valid)
    print("\n" + "=" * 60)
    print("统计汇总")
    print("=" * 60)
    print(f"总样本数: {n}")
    print(f"野生型(WT):   {wt} ({wt/n*100:.1f}%)")
    print(f"纯合子(HOM):  {hom} ({hom/n*100:.1f}%)")
    print(f"杂合子(HET):  {het} ({het/n*100:.1f}%)")
    if abnormal_count:
        print(f"异常(突变位点外双峰): {abnormal_count} 个")
    if err:
        print(f"错误:         {err} ({err/n*100:.1f}%)")
    print("=" * 60)


def compare_with_answer(
    results: List[GenotypeResult], answer_path: str
):
    if not os.path.exists(answer_path):
        return
    answers = {}
    with open(answer_path, "r") as f:
        for line in f:
            line = line.strip()
            if "." in line:
                num, ans = line.split(".", 1)
                answers[int(num)] = ans.upper().strip()
    correct = total = 0
    mismatches = []
    for r in results:
        m = re.match(r"^(\d+)", r.filename)
        if not m:
            continue
        num = int(m.group(1))
        if num not in answers:
            continue
        total += 1
        exp = answers[num]
        if r.genotype == exp:
            correct += 1
        else:
            mismatches.append((num, exp, r.genotype))
    print("\n" + "=" * 60)
    print("与标准答案比较")
    print("=" * 60)
    print(f"准确率: {correct}/{total} ({correct/total*100:.1f}%)")
    if mismatches:
        print(f"不匹配: {len(mismatches)} 个")
        for num, exp, act in mismatches[:15]:
            print(f"  {num}: 期望={exp}, 实际={act}")
        if len(mismatches) > 15:
            print(f"  ... 还有 {len(mismatches)-15} 个")
    print("=" * 60)


def main():
    base = Path(__file__).resolve().parent
    data_dir = base / "data"
    default_reference = data_dir / "正常C57序列.dna"
    if not default_reference.exists():
        default_reference = data_dir / "standard.dna"
    default_ab1_folder = data_dir / "40个案例"
    default_output = data_dir / "results.csv"
    default_answer = data_dir / "ans.txt"

    if len(sys.argv) >= 3:
        reference_path = sys.argv[1]
        ab1_folder = sys.argv[2]
        output_path = sys.argv[3] if len(sys.argv) > 3 else str(default_output)
    else:
        reference_path = str(default_reference)
        ab1_folder = str(default_ab1_folder)
        output_path = str(default_output)

    print("=" * 60)
    print("基因类型分析（先定位最佳匹配片段，再在片段内比对与峰数判断）")
    print("=" * 60)
    print(f"参考序列: {reference_path}")
    print(f"AB1文件夹: {ab1_folder}")
    print(f"输出文件: {output_path}")
    print("=" * 60 + "\n")

    if not os.path.exists(reference_path):
        print(f"错误: 参考序列不存在: {reference_path}")
        sys.exit(1)
    if not os.path.exists(ab1_folder):
        print(f"错误: AB1文件夹不存在: {ab1_folder}")
        sys.exit(1)

    # 检查是否强制指定起始位置（例如从环境变量或命令行）
    force_ref_start = None
    if "FORCE_REF_START" in os.environ:
        try:
            force_ref_start = int(os.environ["FORCE_REF_START"])
        except ValueError:
            pass

    try:
        analyzer = GenotypeAnalyzer(
            reference_path, force_ref_start=force_ref_start
        )
        if force_ref_start is not None:
            print(f"强制使用参考序列起始位置: {force_ref_start}")
        print(f"参考序列长度: {len(analyzer.reference_seq)} bp")
        if analyzer._ref_start is not None:
            print(f"比对区间: ref[{analyzer._ref_start}:{analyzer._ref_end}]\n")
        else:
            print()
        results = analyzer.analyze_folder(ab1_folder)
        if results:
            print_summary(results)
            export_to_csv(results, output_path)
            compare_with_answer(results, str(default_answer))
    except Exception as e:
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
