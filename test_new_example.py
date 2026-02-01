#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
测试 newExample 文件夹：对其中 .ab1 做基因型分析，并与 ans2.txt 标准答案比对。
答案文件中 "ERROR" 视为 "ERR"（程序异常/错误时的返回值）。
"""

import os
import re
from pathlib import Path

from genotype_analyzer import (
    GenotypeAnalyzer,
    GenotypeResult,
    print_summary,
)


def load_answers_normalized(answer_path: str) -> dict:
    """加载答案文件，将 ERROR 规范为 ERR 以便与程序返回值一致。"""
    answers = {}
    if not os.path.exists(answer_path):
        return answers
    with open(answer_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if "." in line:
                num_str, ans = line.split(".", 1)
                num = int(num_str.strip())
                ans = ans.upper().strip()
                if ans == "ERROR":
                    ans = "ERR"
                answers[num] = ans
    return answers


def compare_with_answers(results, answers):
    """与已加载的答案字典比对，返回 (正确数, 参与比对数, 不匹配列表)。"""
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
    return correct, total, mismatches


def main():
    base = Path(__file__).resolve().parent
    data_dir = base / "data"
    reference = data_dir / "standard.dna"
    ab1_folder = data_dir / "newExample"
    answer_path = data_dir / "ans2.txt"

    if not reference.exists():
        print(f"错误: 参考序列不存在: {reference}")
        return
    if not ab1_folder.is_dir():
        print(f"错误: 文件夹不存在: {ab1_folder}")
        return

    print("=" * 60)
    print("newExample 测试（答案: ans2.txt）")
    print("=" * 60)
    print(f"参考序列: {reference}")
    print(f"AB1 文件夹: {ab1_folder}")
    print(f"答案文件: {answer_path}")
    print("=" * 60 + "\n")

    analyzer = GenotypeAnalyzer(str(reference))
    results = analyzer.analyze_folder(str(ab1_folder))

    if not results:
        print("未得到任何结果。")
        return

    print_summary(results)

    answers = load_answers_normalized(str(answer_path))
    if not answers:
        print("\n未找到答案文件或答案为空，跳过准确率比较。")
        return

    correct, total, mismatches = compare_with_answers(results, answers)
    print("\n" + "=" * 60)
    print("与 ans2.txt 比较")
    print("=" * 60)
    if total == 0:
        print("没有样本在答案文件中有对应条目（请检查文件名是否以数字开头）。")
    else:
        print(f"准确率: {correct}/{total} ({100.0 * correct / total:.1f}%)")
        if mismatches:
            print(f"不匹配: {len(mismatches)} 个")
            for num, exp, act in mismatches:
                print(f"  {num}: 期望={exp}, 实际={act}")
    print("=" * 60)


if __name__ == "__main__":
    main()
