#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
序列比对工具：使用parasail进行Smith-Waterman局部比对。

这个实现模拟了SnapGene的比对逻辑：
1. 首先进行完整的Smith-Waterman局部比对
2. 然后找到CIGAR中最长的连续匹配区域作为"核心比对区域"
3. 在核心区域及其附近检测突变

支持的库：
- parasail (首选)：高性能SIMD实现的Smith-Waterman
- pairwise-sequence-alignment (psa)：EMBOSS water
- Biopython PairwiseAligner (回退)
"""

from typing import List, Tuple, Optional

# parasail - 首选库
try:
    import parasail
    HAS_PARASAIL = True
except ImportError:
    HAS_PARASAIL = False

# 可选：pairwise-sequence-alignment (EMBOSS needle/water)
try:
    import psa
    HAS_PSA = True
except ImportError:
    HAS_PSA = False

from Bio import pairwise2
from Bio.Align import PairwiseAligner


def _reverse_complement(seq: str) -> str:
    complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
    return "".join(complement.get(b, "N") for b in reversed(seq))


# ==================== parasail 相关函数 ====================

def _parse_cigar(cigar_str: str) -> List[Tuple[int, str]]:
    """解析CIGAR字符串，返回(长度, 操作)列表"""
    operations = []
    current_len = ""
    for c in cigar_str:
        if c.isdigit():
            current_len += c
        else:
            if current_len:
                operations.append((int(current_len), c))
            current_len = ""
    return operations


def _find_longest_match_region(cigar_str: str, beg_ref: int, beg_query: int) -> dict:
    """
    从CIGAR字符串中找到最长的连续匹配区域。
    这是SnapGene显示的"核心比对区域"。
    """
    operations = _parse_cigar(cigar_str)
    
    best_match_len = 0
    best_ref_start = 0
    best_query_start = 0
    
    ref_pos = beg_ref
    query_pos = beg_query
    
    for length, op in operations:
        if op == '=':  # 匹配
            if length > best_match_len:
                best_match_len = length
                best_ref_start = ref_pos
                best_query_start = query_pos
            ref_pos += length
            query_pos += length
        elif op == 'X':  # 错配
            ref_pos += length
            query_pos += length
        elif op == 'I':  # 插入
            query_pos += length
        elif op == 'D':  # 删除
            ref_pos += length
    
    return {
        "match_length": best_match_len,
        "ref_start": best_ref_start,
        "ref_end": best_ref_start + best_match_len,
        "query_start": best_query_start,
        "query_end": best_query_start + best_match_len,
    }


def _get_mismatches_from_cigar(cigar_str: str, beg_ref: int, beg_query: int,
                                ref_seq: str, query_seq: str) -> List[Tuple[int, int, str, str]]:
    """从CIGAR字符串中提取错配位置"""
    operations = _parse_cigar(cigar_str)
    mismatches = []
    
    ref_pos = beg_ref
    query_pos = beg_query
    
    for length, op in operations:
        if op == '=':
            ref_pos += length
            query_pos += length
        elif op == 'X':
            for i in range(length):
                if ref_pos < len(ref_seq) and query_pos < len(query_seq):
                    mismatches.append((
                        ref_pos,
                        query_pos,
                        ref_seq[ref_pos],
                        query_seq[query_pos]
                    ))
                ref_pos += 1
                query_pos += 1
        elif op == 'I':
            query_pos += length
        elif op == 'D':
            ref_pos += length
    
    return mismatches


def align_with_parasail(ref: str, query: str, 
                        gap_open: int = 10, gap_extend: int = 1) -> dict:
    """
    使用parasail进行Smith-Waterman局部比对。
    返回: {score, cigar, beg_ref, beg_query}
    """
    if not HAS_PARASAIL:
        return {"score": 0, "cigar": None, "beg_ref": 0, "beg_query": 0}
    
    matrix = parasail.matrix_create("ACGT", 2, -1)
    result = parasail.sw_trace_scan_16(query, ref, gap_open, gap_extend, matrix)
    
    cigar_str = None
    beg_query = None
    beg_ref = None
    
    try:
        if result.cigar is not None:
            cigar_str = result.cigar.decode
            beg_query = result.cigar.beg_query
            beg_ref = result.cigar.beg_ref
    except:
        pass
    
    return {
        "score": result.score,
        "cigar": cigar_str,
        "beg_query": beg_query,
        "beg_ref": beg_ref,
    }


def align_and_get_core_region(ref: str, query: str,
                               try_reverse_complement: bool = True,
                               only_core_mismatches: bool = False) -> Tuple[
    List[Tuple[int, int, str, str]], bool, float, int, int, int, int, int]:
    """
    使用parasail进行比对，并找到核心匹配区域（模拟SnapGene）。
    SnapGene 会隐藏低质量末端，实际参与比对的待测序列从“核心区域”对应的 query 位置开始；
    双峰扫描应只在该有效区间内进行，避免舍弃区（如 83、84）的噪声被判为杂合。
    
    Args:
        ref: 参考序列
        query: query序列
        try_reverse_complement: 是否尝试反向互补
        only_core_mismatches: 如果True，只返回核心区域附近的错配
    
    Returns:
        (mismatches, used_rc, score, ref_start, ref_end, core_length, query_start, query_end)
        - mismatches: 错配列表（如果only_core_mismatches=True，只包含核心区域附近的错配）
        - used_rc: 是否使用了反向互补
        - score: 比对分数
        - ref_start: 核心区域的参考序列起始位置(0-based)
        - ref_end: 核心区域的参考序列结束位置(0-based, exclusive)
        - core_length: 核心区域长度
        - query_start: 核心区域在 query 上的起始位置(0-based)，即“有效比对起始”
        - query_end: 核心区域在 query 上的结束位置(0-based, exclusive)
    """
    if not HAS_PARASAIL:
        # 如果没有parasail，回退到其他方法（无 query_start/query_end，用全序列 0, len(query)）
        mm, urc, sc, rs, re = align_full_and_get_mismatches_with_position(
            ref, query, try_reverse_complement
        )
        return (mm, urc, sc, rs, re, re - rs, 0, len(query))
    
    # 正向比对
    result_fwd = align_with_parasail(ref, query)
    best_result = result_fwd
    best_query = query
    used_rc = False
    
    # 反向互补比对
    if try_reverse_complement:
        query_rc = _reverse_complement(query)
        result_rc = align_with_parasail(ref, query_rc)
        
        if result_rc["score"] > result_fwd["score"]:
            best_result = result_rc
            best_query = query_rc
            used_rc = True
    
    # 如果没有有效的CIGAR，回退（无 query_start/query_end，用全序列）
    if not best_result["cigar"]:
        mm, urc, sc, rs, re = align_full_and_get_mismatches_with_position(
            ref, query, try_reverse_complement
        )
        return (mm, urc, sc, rs, re, re - rs, 0, len(query))
    
    # 找到核心匹配区域
    core_region = _find_longest_match_region(
        best_result["cigar"],
        best_result["beg_ref"],
        best_result["beg_query"]
    )
    
    # 获取完整比对中的所有错配
    all_mismatches = _get_mismatches_from_cigar(
        best_result["cigar"],
        best_result["beg_ref"],
        best_result["beg_query"],
        ref,
        best_query
    )
    
    # 如果只需要核心区域附近的错配
    if only_core_mismatches:
        # 核心区域附近：向前扩展20bp，向后扩展20bp
        margin = 20
        core_start = core_region["ref_start"] - margin
        core_end = core_region["ref_end"] + margin
        all_mismatches = [
            m for m in all_mismatches
            if core_start <= m[0] < core_end
        ]
    
    return (
        all_mismatches,
        used_rc,
        float(best_result["score"]),
        core_region["ref_start"],
        core_region["ref_end"],
        core_region["match_length"],
        core_region["query_start"],
        core_region["query_end"],
    )


def align_semiglobal_pairwise2(
    ref: str, query: str
) -> Tuple[str, str, float]:
    """
    半全局比对（两端空位不罚分）。使用 one_alignment_only 避免长序列时超时。
    返回 (aligned_ref, aligned_query, score)。
    """
    alignments = pairwise2.align.globalms(
        ref, query, 2, -1, -2, -0.5,
        penalize_end_gaps=False,
        one_alignment_only=True,
    )
    if not alignments:
        return "", "", 0.0
    a = alignments[0]
    return a.seqA, a.seqB, a.score


def align_local_psa(
    ref: str, query: str
) -> Tuple[str, str, float, int, int]:
    """
    使用 psa (EMBOSS water) 做局部比对，得到最佳匹配片段。
    返回 (aligned_ref, aligned_query, score, ref_start_0based, ref_end_exclusive)。
    psa 使用 1-based：sstart/send 为 subject(ref) 的起止。
    """
    # water: query 与 subject 的局部比对，subject 为参考
    aln = psa.water(moltype="nucl", qseq=query, sseq=ref)
    # saln = subject aligned, qaln = query aligned
    ref_start = aln.sstart - 1  # 转 0-based
    ref_end = aln.send  # 1-based inclusive，片段为 ref[ref_start:ref_end]
    return aln.saln, aln.qaln, float(aln.score), ref_start, ref_end


def align_local_biopython(
    ref: str, query: str
) -> Tuple[str, str, float, int, int]:
    """
    使用 Bio.Align.PairwiseAligner 局部比对（Smith-Waterman），得到最佳匹配片段。
    返回 (aligned_ref, aligned_query, score, ref_start_0based, ref_end_exclusive)。
    """
    aligner = PairwiseAligner(mode="local")
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    # 只取第一个最优比对，避免长序列时最优解过多导致溢出
    alignments = aligner.align(ref, query)
    aln = next(alignments, None)
    if aln is None:
        return "", "", 0.0, 0, 0
    # aligned[0] 为 ref 的 [start,end) 段列表，aligned[1] 为 query
    ref_segments = aln.aligned[0]
    ref_start = min(s[0] for s in ref_segments)
    ref_end = max(s[1] for s in ref_segments)
    ar = str(aln[0])
    aq = str(aln[1])
    return ar, aq, aln.score, ref_start, ref_end


def get_ref_segment_from_alignment(
    aligned_ref: str, aligned_query: str
) -> Tuple[int, int]:
    """
    从比对结果得到参考序列上与 query 实际重叠的区间 [ref_start, ref_end)（0-based）。
    只取 ref 与 query 均非空位的位置（即真正比对的区间）。
    """
    ref_positions = []
    ir = 0
    for ar, aq in zip(aligned_ref, aligned_query):
        if ar != "-" and aq != "-":
            ref_positions.append(ir)
        if ar != "-":
            ir += 1
    if not ref_positions:
        return 0, 0
    return min(ref_positions), max(ref_positions) + 1


def get_mismatches_from_alignment(
    aligned_ref: str, aligned_query: str
) -> List[Tuple[int, int, str, str]]:
    """
    从比对结果提取错配：(ref_0based, query_0based, ref_base, query_base)。
    """
    mismatches = []
    ir, iq = 0, 0
    for ar, aq in zip(aligned_ref, aligned_query):
        if ar != "-" and aq != "-":
            if ar != aq:
                mismatches.append((ir, iq, ar, aq))
            ir += 1
            iq += 1
        elif ar != "-":
            ir += 1
        elif aq != "-":
            iq += 1
    return mismatches


def _mismatches_from_psa_alignment(aln) -> List[Tuple[int, int, str, str]]:
    """从 psa 的 alignment 对象提取错配 (ref_0based, query_0based, ref_base, query_base)。"""
    mismatches = []
    ref_pos = aln.sstart - 1
    query_pos = aln.qstart - 1
    for ar, aq in zip(aln.saln, aln.qaln):
        if ar != "-" and aq != "-":
            if ar != aq:
                mismatches.append((ref_pos, query_pos, ar, aq))
            ref_pos += 1
            query_pos += 1
        elif ar != "-":
            ref_pos += 1
        elif aq != "-":
            query_pos += 1
    return mismatches


def align_full_and_get_mismatches_with_position(
    ref: str, query: str, try_reverse_complement: bool = True
) -> Tuple[List[Tuple[int, int, str, str]], bool, float, int, int]:
    """
    全参考序列比对，返回错配列表与参考匹配区间。
    若已安装 psa 与 EMBOSS，优先用 psa.water 局部比对（更快且直接得到最佳片段）。
    否则用 Biopython pairwise2 半全局比对。
    返回 (mismatches, used_reverse_complement, score, ref_start, ref_end)。
    """
    best_mismatches: List[Tuple[int, int, str, str]] = []
    best_score = float("-inf")
    used_rc = False
    best_ref_start, best_ref_end = 0, 0

    if HAS_PSA:
        try:
            # 正向
            aln = psa.water(moltype="nucl", qseq=query, sseq=ref)
            score = float(aln.score)
            ref_start = aln.sstart - 1
            ref_end = aln.send
            if score > best_score:
                best_score = score
                best_mismatches = _mismatches_from_psa_alignment(aln)
                best_ref_start, best_ref_end = ref_start, ref_end
                used_rc = False
            # 反向互补
            if try_reverse_complement:
                query_rc = _reverse_complement(query)
                aln_rc = psa.water(moltype="nucl", qseq=query_rc, sseq=ref)
                if float(aln_rc.score) > best_score:
                    best_score = float(aln_rc.score)
                    best_mismatches = _mismatches_from_psa_alignment(aln_rc)
                    best_ref_start = aln_rc.sstart - 1
                    best_ref_end = aln_rc.send
                    used_rc = True
            return (
                best_mismatches,
                used_rc,
                best_score,
                best_ref_start,
                best_ref_end,
            )
        except Exception:
            pass

    # 回退到 Biopython 局部比对（Smith-Waterman），精确定位最佳匹配片段
    ar, aq, score, ref_start, ref_end = align_local_biopython(ref, query)
    if score > best_score and ref_end > ref_start:
        best_score = score
        best_mismatches = get_mismatches_from_alignment(ar, aq)
        best_ref_start, best_ref_end = ref_start, ref_end
        used_rc = False

    if try_reverse_complement:
        query_rc = _reverse_complement(query)
        ar_rc, aq_rc, score_rc, ref_start_rc, ref_end_rc = align_local_biopython(
            ref, query_rc
        )
        if score_rc > best_score and ref_end_rc > ref_start_rc:
            best_score = score_rc
            best_mismatches = get_mismatches_from_alignment(ar_rc, aq_rc)
            best_ref_start, best_ref_end = ref_start_rc, ref_end_rc
            used_rc = True

    return (
        best_mismatches,
        used_rc,
        best_score,
        best_ref_start,
        best_ref_end,
    )


def align_to_segment_and_get_mismatches(
    ref: str,
    ref_start: int,
    ref_end: int,
    query: str,
    try_reverse_complement: bool = True,
) -> Tuple[List[Tuple[int, int, str, str]], bool, float]:
    """
    仅用参考序列的片段 ref[ref_start:ref_end] 与 query 做半全局比对。
    错配中的 ref 位置已换算为全参考上的 0-based 坐标。
    返回 (mismatches, used_reverse_complement, score)。
    """
    segment = ref[ref_start:ref_end]
    ar, aq, score = align_semiglobal_pairwise2(segment, query)
    used_rc = False

    if try_reverse_complement:
        query_rc = _reverse_complement(query)
        ar_rc, aq_rc, score_rc = align_semiglobal_pairwise2(segment, query_rc)
        if score_rc > score:
            ar, aq, score = ar_rc, aq_rc, score_rc
            used_rc = True

    raw_mismatches = get_mismatches_from_alignment(ar, aq)
    # 将 segment 内的相对 ref 位置转为全参考位置
    mismatches = [
        (ref_start + r, q, rb, qb) for r, q, rb, qb in raw_mismatches
    ]
    return mismatches, used_rc, score
