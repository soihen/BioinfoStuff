__author__ = 'Kai'
__date__ = '18/03/2019'
__email__ = 'zhentian.kai@outlook.com'

'''
    Two core functions:
    1. correct HGVS name using https://www.mutalyzer.nl
    2. translate HGVS name to Chinese...
'''


from bs4 import BeautifulSoup as bs
import requests
import re



def autocorrection(input_data):
    '''
        Web crawler that correct the HGVS name using https://www.mutalyzer.nl
        N.B. The uncorrected HGVS name maybe so wrong formatted that can even not be recognised by the website :(
        @param: input_data -- uncorrected HGVS name (located on the hgvs命名 column in the annotated result file)
        e.g. CAMTA1:NM_001242701:exon4:c.C238T:p.L80F
             gene : transcript : exon : DNA : protein
    '''
    # reformat the uncorrected name
    searchable = input_data.split(':')[1] + ':' + input_data.split(':')[3]
    # crawling
    r = requests.get('https://www.mutalyzer.nl/name-checker?description=' + searchable)
    html_doc = r.content
    soup = bs(html_doc, 'html.parser')
    findings = soup.find_all('code')
    # get desired texts
    transcript = findings[-3].string.split(":")[1]
    protein = findings[-2].string.split(":")[1]
    protein = re.sub(r'[\(\)]', '', protein)
    outname = input_data.split(":")[:-2]
    outname.append(transcript)
    outname.append(protein)
    return ':'.join(outname)



def translate(hgvs_name):
    '''
        translate HGVS name into Chinese
        @param: hgvs_name -- MUST be correctly formatted HGVS name
        @return: translated name; however it is possible that this function cannot recognise the provided name,
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
        aachinese = '第{}位碱基{}置换为{}使得'.format(m1.group(1), _amino_cn[m1.group(2)], m1.group(3))
        flag1 = True
    if m2:
        aachinese = '第{}位碱基缺失使得'.format(m2.group(1))
        flag1 = True
    if m3:
        aachinese = '第{}位碱基缺失使得'.format(m3.group(1)+'_'+m3.group(2))
        flag1 = True
    if m4:
        aachinese = '第{}位碱基重复使得'.format(m4.group(1))
        flag1 = True
    if m5:
        aachinese = '第{}位碱基重复使得'.format(m5.group(1)+'_'+m5.group(2))
        flag1 = True
    if m6:
        aachinese = '第{}位碱基插入{}使得'.format(m6.group(1)+'_'+m6.group(2), m6.group(3))
        flag1 = True
    if m7:
        aachinese = '第{}位碱基缺失插入{}使得'.format(m7.group(1)+'_'+m7.group(2), m7.group(3))
        flag1 = True
    if not flag1:
        translated_name = '请报告审核者根据hgvs命名提供'

    ### protein part
    # substitution
    m1 = re.match(r'p\.([a-zA-Z]{3})(\d+)(([a-zA-Z]{3}))', pchange)
    # deletion
    m2 = re.match(r'p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)del', pchange)
    # delins
    m3 = re.match(r'p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)delins([a-zA-Z]{3})', pchange)
    # duplication
    m4 = re.match(r'p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)dup', pchange)
    # contain stop codon
    m5 = re.search(r'p\.(\S+)\*(\d+)', pchange)
    # insertion
    m6 = re.match(r'p\.([a-zA-Z]{3})(\d+)_([a-zA-Z]{3})(\d+)ins([a-zA-Z]{3})', pchange)
    flag2 = False
    if m1:
        pchinese = '第{}位氨基酸由{}转变成了{}'.format(m1.group(2), _amino_cn[_amino_d[m1.group(1)]], _amino_cn[_amino_d[m1.group(3)]])
        flag2 = True
    if m2:
        pchinese = '第{}位氨基酸缺失'.format(m2.group(2)+'_'+m2.group(4))
        flag2 = True
    if m3:
        pchinese = '第{}位氨基酸缺失插入{}'.format(m3.group(2)+'_'+m3.group(4), _amino_cn[_amino_d[m3.group(5)]])
        flag2 = True
    if m4:
        pchinese = '第{}位氨基酸重复'.format(m4.group(2)+'_'+m4.group(4))
        flag2 = True
    if flag2 and m5:
        pchinese += '并于第{}个氨基酸提前终止'.format(m5.group(2))
    if m6:
        pchinese = '第{}位氨基酸插入{}'.format(m6.group(2)+'_'+m6.group(4), _amino_cn[_amino_d[m6.group(5)]])
    if not flag2:
        translated_name = '请报告审核者根据hgvs命名提供'
                                                                                 
    # only both flags are true.
    if flag1 and flag2:
        translated_name = gene + '基因' + exon + aachinese + pchinese
    return translated_name



