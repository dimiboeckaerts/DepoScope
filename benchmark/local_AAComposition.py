
# -*- coding: utf-8 -*-
"""
The module is used for computing the composition of amino acids, dipetide and
3-mers (tri-peptide) for a given protein sequence.
References
----------
.. [1] Reczko, M. and Bohr, H. (1994) The DEF data base of sequence based protein
   fold class predictions. Nucleic Acids Res, 22, 3616-3619.
.. [2] Hua, S. and Sun, Z. (2001) Support vector machine approach for protein
   subcellular localization prediction. Bioinformatics, 17, 721-728.
.. [3] Grassmann, J., Reczko, M., Suhai, S. and Edler, L. (1999) Protein fold
   class prediction: new methods of statistical classification. Proc Int Conf
   Intell Syst Mol Biol, 106-112.
Authors: Dongsheng Cao and Yizeng Liang.
Date: 2012.3.27
Email: oriental-cds@163.com
"""

# Core Library
import re
from typing import Any, Dict, List

AALetter: List[str] = list("ARNDCEQGHILKMFPSTWYV")

ProteinSequence_docstring = """ProteinSequence: str
        a pure protein sequence"""


def CalculateAAComposition(ProteinSequence: str) -> Dict[str, float]:
    sequence_length = len(ProteinSequence)
    result: Dict[str, float] = {}
    for i in AALetter:
        result[i] = round(float(ProteinSequence.count(i)) / sequence_length * 100, 3)
    return result


def CalculateDipeptideComposition(ProteinSequence: str) -> Dict[str, float]:
    sequence_length = len(ProteinSequence)
    result = {}
    for i in AALetter:
        for j in AALetter:
            dipeptide = i + j
            result[dipeptide] = round(
                float(ProteinSequence.count(dipeptide)) / (sequence_length - 1) * 100, 2
            )
    return result


def Getkmers() -> List[str]:
    kmers = []
    for i in AALetter:
        for j in AALetter:
            for k in AALetter:
                kmers.append(i + j + k)
    return kmers


def GetSpectrumDict(proteinsequence: str) -> Dict[str, int]:
    result = {}
    kmers = Getkmers()
    for i in kmers:
        result[i] = len(re.findall(i, proteinsequence))
    return result


def CalculateAADipeptideComposition(ProteinSequence: str) -> Dict[str, float]:
    result: Dict[Any, Any] = {}
    result.update(CalculateAAComposition(ProteinSequence))
    result.update(CalculateDipeptideComposition(ProteinSequence))
    result.update(GetSpectrumDict(ProteinSequence))
    return result






