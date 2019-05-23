# -*- coding: utf-8 -*-

__doc__ = """
双击便可直接打开运行，
也可以使用pyinstaller -wF hgvs_translate.py打包，就能直接在windows下运行"""

__author__ = "Kai"
__date__ = "2019/05/22"

import PySimpleGUI as sg
import argparse
import logging
import re


def translate(hgvs_name):
    '''
        translate HGVS name into Chinese
        :param: hgvs_name -- MUST be correctly formatted HGVS name
        :return: translated name; however it is possible that this function cannot recognise the provided name,
        return '请报告审核者根据hgvs命名提供' if that's the case.
    '''
    # Those two dictionary will be used for translation...
    _amino_cn = {'A': '丙氨酸', 'C': '半胱氨酸', 'D': '天冬氨酸', 'E': '谷氨酸',
                 'F': '苯丙氨酸', 'G': '甘氨酸', 'H': '组氨酸', 'I': '异亮氨酸',
                 'K': '赖氨酸', 'L': '亮氨酸', 'M': '甲硫氨酸', 'N': '天冬酰胺',
                 'P': '脯氨酸', 'Q': '谷氨酰胺', 'R': '精氨酸', 'S': '丝氨酸',
                 'T': '苏氨酸', 'V': '缬氨酸', 'W': '色氨酸', 'Y': '酪氨酸', 'X': '终止'}

    _amino_d = {'Gly': 'G', 'Cys': 'C', 'Glu': 'E', 'Asp': 'D',
                'Ile': 'I', 'Pro': 'P', 'Tyr': 'Y', 'Lys': 'K', 'Gln': 'Q',
                'Trp': 'W', 'Leu': 'L', 'Phe': 'F', 'Val': 'V', 'Ser': 'S',
                'Met': 'M', 'Ala': 'A', 'His': 'H', 'Ter': 'X', 'Asn': 'N',
                'Thr': 'T', 'Arg': 'R'}

    if hgvs_name.count(":") != 4:
        return "这个HGVS命名本身有问题，看看是否只有4个分号存在，例：BRCA2:NM_000059:exon11:c.2236G>A:p.Val746Met"
    else:
        gene, transcript, exon, aachange, pchange = hgvs_name.split(":")
        ### cDNA part
        # DNA substitution
        m1 = re.match(r'c\.(\d+)([A-Z])>([A-Z])', aachange)
        # DNA deletion
        m2 = re.match(r'c\.(\d+)del', aachange)
        m3 = re.match(r'c\.(\d+)_(\d+)del', aachange)
        # DNA duplication
        m4 = re.match(r'c\.(\d+)dup', aachange)
        m5 = re.match(r'c\.(\d+)_(\d+)dup', aachange)
        # DNA insertion
        m6 = re.match(r'c\.(\d+)_(\d+)ins([A-Z]+)', aachange)
        # DNA delins
        m7 = re.match(r'c\.(\d+)_(\d+)delins([A-Z]+)', aachange)
        flag1 = False
        if m1:
            aachinese = '第{}位碱基{}置换为{}使得'.format(m1.group(1), m1.group(2), m1.group(3))
            flag1 = True
        if m2:
            aachinese = '第{}位碱基缺失使得'.format(m2.group(1))
            flag1 = True
        if m3:
            aachinese = '第{}位碱基缺失使得'.format(m3.group(1) + '_' + m3.group(2))
            flag1 = True
        if m4:
            aachinese = '第{}位碱基重复使得'.format(m4.group(1))
            flag1 = True
        if m5:
            aachinese = '第{}位碱基重复使得'.format(m5.group(1) + '_' + m5.group(2))
            flag1 = True
        if m6:
            aachinese = '第{}位碱基插入{}使得'.format(m6.group(1) + '_' + m6.group(2), m6.group(3))
            flag1 = True
        if m7:
            aachinese = '第{}位碱基缺失插入{}使得'.format(m7.group(1) + '_' + m7.group(2), m7.group(3))
            flag1 = True
        if not flag1:
            return 'hgvs命名碱基部分出错'

        ### protein part
        # substitution - p.Val368Trp
        m1 = re.match(r'p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})', pchange)
        # substitution
        m1_2 = re.match(r'p\.([A-Z][a-z]{2})(\d+)\*', pchange)
        # deletion - p.Val368_Trp378del
        m2 = re.match(r'p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)del', pchange)
        # deletion - p.Val342del
        m2_2 = re.match(r'p\.([A-Z][a-z]{2})(\d+)del', pchange)
        # delins - p.Val368delinsGlyLys
        m3 = re.match(r'p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)delins(([A-Z][a-z]{2})+)', pchange)
        m7 = re.match(r'p\.([A-Z][a-z]{2})(\d+)delins(([A-Z][a-z]{2})+)', pchange)
        # duplication
        m4 = re.match(r'p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)dup', pchange)
        # contain stop codon
        m5 = re.search(r'p\.(\S+)\*(\d+)', pchange)
        # insertion
        m6 = re.match(r'p\.([A-Z][a-z]{2})(\d+)_([A-Z][a-z]{2})(\d+)ins(([A-Z][a-z]{2})+)', pchange)
        flag2 = False
        # frameshift
        m8 = ('fs' in pchange)
        m9 = re.match(r'p\.([A-Z][a-z]{2})(\d+)dup', pchange)
        if m1:
            pchinese = '第{}位氨基酸由{}转变成了{}'.format(m1.group(2), _amino_cn[_amino_d[m1.group(1)]],
                                                 _amino_cn[_amino_d[m1.group(3)]])
            flag2 = True
        if m1_2:
            pchinese = '第{}位氨基酸由{}转变成了终止密码子'.format(m1_2.group(2), _amino_cn[_amino_d[m1_2.group(1)]])
            flag2 = True
        if m2:
            pchinese = '第{}位氨基酸缺失'.format(m2.group(2) + '_' + m2.group(4))
            flag2 = True
        if m2_2:
            pchinese = '第{}的{}缺失'.format(m2_2.group(2), _amino_cn[_amino_d[m2_2.group(1)]])
            flag2 = True
        if m3:
            inserted_aachinse = ''
            inserted_aa = str(m3.group(5))  # e.g. AlaGly
            for i in range(0, len(inserted_aa), 3):
                each_aa = inserted_aa[i:i + 3]
                inserted_aachinse += _amino_cn[_amino_d[each_aa]]
            pchinese = '第{}位氨基酸缺失插入{}'.format(m3.group(2) + '_' + m3.group(4), inserted_aachinse)
            flag2 = True
        if m7:
            inserted_aachinse = ''
            inserted_aa = str(m7.group(3))  # e.g. AlaGly
            for i in range(0, len(inserted_aa), 3):
                each_aa = inserted_aa[i:i + 3]
                inserted_aachinse += _amino_cn[_amino_d[each_aa]]
            pchinese = '第{}位氨基酸缺失插入{}'.format(m7.group(2), inserted_aachinse)
            flag2 = True
        if m4:
            pchinese = '第{}位氨基酸重复'.format(m4.group(2) + '_' + m4.group(4))
            flag2 = True
        if flag2 and m5:
            pchinese += '并于第{}个氨基酸提前终止'.format(m5.group(2))
        if m6:
            inserted_aachinse = ''
            inserted_aa = str(m6.group(5))  # e.g. AlaGly
            for i in range(0, len(inserted_aa), 3):
                each_aa = inserted_aa[i:i + 3]
                inserted_aachinse += _amino_cn[_amino_d[each_aa]]
            pchinese = '第{}位氨基酸插入{}'.format(m6.group(2) + '_' + m6.group(4), inserted_aachinse)
            flag2 = True
        if m9:
            pchinese = '第{}位{}重复'.format(m9.group(2), _amino_cn[_amino_d[m9.group(1)]])
            flag2 = True
        if not flag2:
            return 'hgvs命名氨基酸部分出错'
        # only both flags are true.
        if flag1 and flag2:
            translated_name = gene + '基因' + exon + aachinese + pchinese
        if m8:
            translated_name += '且发生移码'
        return translated_name



if __name__ == '__main__':
    # make gui...
    layout = [[sg.Text("在下方输入HGVS英文命名")],
              [sg.Text("HGVS英文"), sg.InputText('')],
              [sg.Text("HGVS中文"), sg.Output(size=(45, 3), key='_OUTPUT_')],
              [sg.Button("翻译")]]

    window = sg.Window("HGVS翻译～", layout)
    while True:
        event, values = window.Read()
        if event == '翻译':
            window.Element('_OUTPUT_').Update(translate(values[0]))
        elif event is None:
            break
    window.Close()





