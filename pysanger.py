"""
PySanger - A Python module to analyze Sanger sequencing results.
From: https://github.com/ponnhide/PySanger
"""

import os 
import sys 
from Bio import SeqIO
from Bio import pairwise2

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}

def abi_to_dict(filename):
    """
    Generate confidence value array and signal intensity arrays of each channel (A, T, G or C) at peak positions.
    """
    record = SeqIO.read(filename, 'abi')
    abi_data = {
        "conf": [],
        "channel": {"A": [], "T": [], "G": [], "C": []},
        "_channel": {"A": [], "T": [], "G": [], "C": []}
    }
    
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        if pos > 4 and pos < len(record.annotations['abif_raw']["DATA9"]) - 5: 
            abi_data["conf"].append(conf)
            abi_data["channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])

            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-5])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+5])
            
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-5])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+5])
            
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-5])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+5])

            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-5])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+5])

    return abi_data 


def generate_consensusseq(abidata):
    """
    Generate the most consensus seq from a sanger sequencing result.
    """
    consensus_seq = "" 
    
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        consensus_seq += _atgc_dict[values.index(max(values))]
    
    return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC", "TACG"))[::-1]) 


def generate_pwm(abidata):
    """
    Generate position weight matrix based on signal intensities of each channel.
    """
    pwm = {"A": [], "T": [], "G": [], "C": []} 
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        v = 100000 / (sum(values) + 1) 
        new_values = (v * values[0], v * values[1], v * values[2], v * values[3])
        new_values = list(map(int, new_values))
        
        while sum(new_values) < 100000:
            for i in range(len(new_values)):
                new_values[i] += 1
                if sum(new_values) == 100000:
                    break 
        
        pwm["A"].append(new_values[0])
        pwm["T"].append(new_values[1])
        pwm["G"].append(new_values[2])
        pwm["C"].append(new_values[3])

    try:
        import pandas as pd
        return pd.DataFrame(pwm)
    except ImportError:
        return pwm  # 无 pandas 时返回 dict，便于打包 exe 不依赖 pandas
