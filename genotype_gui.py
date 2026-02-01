#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
基因型分析可视化界面

选择 .dna 参考序列，以及 .ab1 测序文件（多选）或文件夹（自动查找所有 .ab1），输出基因型结果。
"""

import os
import re
import csv
import sys
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from pathlib import Path
from typing import List

# 确保可导入同目录下的 genotype_analyzer
sys.path.insert(0, str(Path(__file__).resolve().parent))
from genotype_analyzer import GenotypeAnalyzer, GenotypeResult, export_to_csv


def center_window(window: tk.Tk, width: int = 640, height: int = 520):
    """窗口居中"""
    screen_w = window.winfo_screenwidth()
    screen_h = window.winfo_screenheight()
    x = (screen_w - width) // 2
    y = (screen_h - height) // 2
    window.geometry(f"{width}x{height}+{x}+{y}")


def _extract_number(filename: str) -> int:
    m = re.match(r"^(\d+)", filename)
    return int(m.group(1)) if m else 0


def collect_ab1_from_folder(folder_path: str) -> List[str]:
    """从文件夹中收集所有 .ab1 / .AB1 文件路径，按文件名数字排序"""
    folder = Path(folder_path)
    paths = set()
    for f in folder.glob("*.ab1"):
        paths.add(str(f.resolve()))
    for f in folder.glob("*.AB1"):
        paths.add(str(f.resolve()))
    return sorted(list(paths), key=lambda p: _extract_number(Path(p).name))


class GenotypeGUI:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("基因型分析")
        self.root.resizable(True, True)
        center_window(self.root, 640, 520)

        self.dna_path = tk.StringVar(value="")
        self.ab1_paths: List[str] = []  # 当前选中的 .ab1 路径列表
        self.ab1_display = tk.StringVar(value="未选择")  # 用于显示在界面
        self.last_results: List[GenotypeResult] = []

        self._build_ui()

    def _build_ui(self):
        main = ttk.Frame(self.root, padding=16)
        main.pack(fill=tk.BOTH, expand=True)

        # 参考序列 .dna
        ttk.Label(main, text="参考序列 (.dna):").pack(anchor=tk.W)
        row1 = ttk.Frame(main)
        row1.pack(fill=tk.X, pady=(0, 12))
        ttk.Entry(row1, textvariable=self.dna_path, width=55).pack(
            side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 8)
        )
        ttk.Button(row1, text="选择文件", command=self._choose_dna).pack(side=tk.RIGHT)

        # 测序文件：多选 或 文件夹
        ttk.Label(main, text="测序文件 (.ab1):").pack(anchor=tk.W)
        row2 = ttk.Frame(main)
        row2.pack(fill=tk.X, pady=(0, 4))
        ttk.Label(row2, textvariable=self.ab1_display, foreground="gray").pack(
            side=tk.LEFT, fill=tk.X, expand=True
        )
        row2_btns = ttk.Frame(main)
        row2_btns.pack(fill=tk.X, pady=(0, 12))
        ttk.Button(row2_btns, text="选择文件(可多选)", command=self._choose_ab1_multi).pack(
            side=tk.LEFT, padx=(0, 8)
        )
        ttk.Button(row2_btns, text="选择文件夹", command=self._choose_ab1_folder).pack(
            side=tk.LEFT
        )

        # 分析 + 导出
        row3 = ttk.Frame(main)
        row3.pack(fill=tk.X, pady=(0, 12))
        ttk.Button(row3, text="开始分析", command=self._run_analysis).pack(side=tk.LEFT, padx=(0, 8))
        self.export_btn = ttk.Button(row3, text="导出 CSV", command=self._export_csv, state=tk.DISABLED)
        self.export_btn.pack(side=tk.LEFT)

        # 结果表格
        ttk.Label(main, text="分析结果:").pack(anchor=tk.W)
        table_frame = ttk.Frame(main)
        table_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 8))

        columns = ("filename", "genotype", "genotype_cn", "mutation_count", "abnormal")
        self.tree = ttk.Treeview(table_frame, columns=columns, show="headings", height=12)
        self.tree.heading("filename", text="文件名")
        self.tree.heading("genotype", text="基因类型")
        self.tree.heading("genotype_cn", text="中文")
        self.tree.heading("mutation_count", text="差异数")
        self.tree.heading("abnormal", text="异常")
        self.tree.column("filename", width=280)
        self.tree.column("genotype", width=60)
        self.tree.column("genotype_cn", width=60)
        self.tree.column("mutation_count", width=60)
        self.tree.column("abnormal", width=140)
        scroll_y = ttk.Scrollbar(table_frame)
        scroll_x = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL)
        self.tree.configure(yscrollcommand=scroll_y.set, xscrollcommand=scroll_x.set)
        scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # 汇总
        self.summary_label = ttk.Label(main, text="", font=("", 10))
        self.summary_label.pack(anchor=tk.W, pady=(4, 0))

    def _choose_dna(self):
        path = filedialog.askopenfilename(
            title="选择参考序列",
            filetypes=[("SnapGene DNA", "*.dna"), ("所有文件", "*.*")],
        )
        if path:
            self.dna_path.set(path)

    def _choose_ab1_multi(self):
        paths = filedialog.askopenfilenames(
            title="选择测序文件（可多选）",
            filetypes=[("AB1 测序", "*.ab1 *.AB1"), ("所有文件", "*.*")],
        )
        if paths:
            self.ab1_paths = sorted(list(paths), key=lambda p: _extract_number(Path(p).name))
            n = len(self.ab1_paths)
            self.ab1_display.set(f"已选 {n} 个文件" if n > 1 else self.ab1_paths[0])

    def _choose_ab1_folder(self):
        path = filedialog.askdirectory(title="选择包含 .ab1 的文件夹")
        if path:
            self.ab1_paths = collect_ab1_from_folder(path)
            n = len(self.ab1_paths)
            if n == 0:
                messagebox.showinfo("提示", "该文件夹中未找到 .ab1 或 .AB1 文件。")
                self.ab1_display.set("未选择")
            else:
                self.ab1_display.set(f"文件夹: {path}（共 {n} 个 .ab1）")

    def _run_analysis(self):
        dna = self.dna_path.get().strip()
        if not dna:
            messagebox.showwarning("提示", "请先选择参考序列 (.dna) 文件。")
            return
        if not os.path.isfile(dna):
            messagebox.showerror("错误", f"参考序列文件不存在:\n{dna}")
            return
        if not self.ab1_paths:
            messagebox.showwarning("提示", "请先选择测序文件（多选）或选择文件夹。")
            return

        # 清空表格
        for item in self.tree.get_children():
            self.tree.delete(item)
        self.summary_label.config(text="分析中...")
        self.export_btn.config(state=tk.DISABLED)
        self.root.update()

        try:
            analyzer = GenotypeAnalyzer(dna)
            results: List[GenotypeResult] = []
            for i, ab1_path in enumerate(self.ab1_paths):
                self.summary_label.config(text=f"分析中 ({i+1}/{len(self.ab1_paths)})...")
                self.root.update()
                result = analyzer.analyze_sample(ab1_path)
                results.append(result)
                abnormal_str = ""
                if result.is_abnormal:
                    pos_str = ",".join(map(str, result.abnormal_positions[:6]))
                    if len(result.abnormal_positions) > 6:
                        pos_str += f"...({len(result.abnormal_positions)}个)"
                    abnormal_str = f"是(碱基{pos_str})"
                self.tree.insert(
                    "",
                    tk.END,
                    values=(
                        result.filename,
                        result.genotype,
                        result.genotype_cn,
                        result.mutation_count,
                        abnormal_str,
                    ),
                )

            self.last_results = results
            valid = [r for r in results if r.genotype != "ERR"]
            n = len(results)
            wt = sum(1 for r in valid if r.genotype == "WT")
            hom = sum(1 for r in valid if r.genotype == "HOM")
            het = sum(1 for r in valid if r.genotype == "HET")
            err = n - len(valid)
            abnormal_count = sum(1 for r in results if r.is_abnormal)
            summary = f"共 {n} 个 | WT: {wt} | HOM: {hom} | HET: {het}"
            if abnormal_count:
                summary += f" | 异常: {abnormal_count}"
            if err:
                summary += f" | 错误: {err}"
            self.summary_label.config(text=summary)
            self.export_btn.config(state=tk.NORMAL)
        except Exception as e:
            self.summary_label.config(text="分析出错")
            self.export_btn.config(state=tk.DISABLED)
            messagebox.showerror("错误", str(e))

    def _export_csv(self):
        if not self.last_results:
            messagebox.showinfo("提示", "暂无结果可导出，请先完成分析。")
            return
        path = filedialog.asksaveasfilename(
            title="导出 CSV",
            defaultextension=".csv",
            filetypes=[("CSV", "*.csv"), ("所有文件", "*.*")],
        )
        if path:
            try:
                export_to_csv(self.last_results, path)
                messagebox.showinfo("完成", f"已导出到:\n{path}")
            except Exception as e:
                messagebox.showerror("错误", str(e))

    def run(self):
        self.root.mainloop()


if __name__ == "__main__":
    app = GenotypeGUI()
    app.run()
