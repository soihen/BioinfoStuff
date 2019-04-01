## Pysam — VCF

What is VCF format?:   [VCF format details](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40)



pysam is a Python module allowing easy manipulation of SAM/BAM and VCF/BCF files

Manipulation of VCF files is mainly through class `pysam.VariantFile`

Here, I'm gonna demonstrate some basic usages of  `pysam.VariantFile`

---

The `VariantFile` class contain two sub-classes: 

1. `VariantHeader`   —  for VCF header 

2. `VariantRecord`   —  for VCF record



Initiate a `VariantFile` object:

```python
import pysam

# This command opens the VCF file directly
vcf = pysam.VariantFile('MySample.vcf', 'r')
'''
	do some stuff here
'''
vcf.close()
```



---

#### VariantHeader

```python
import pysam

vcf = pysam.VariantFile('MySample.vcf', 'r')
# create header object
header = vcf.header
# get all the header rows
print(header)
# get info, samples, filters, configs
print(list(header.info), list(header.samples), list(header.filters), list(header.contigs))

# get all attributes 
print(dir(pysam.libcbcf.VariantHeader))
# output:
['__bool__', '__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__pyx_vtable__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', 'add_line', 'add_meta', 'add_record', 'add_sample', 'alts', 'contigs', 'copy', 'filters', 'formats', 'info', 'merge', 'new_record', 'records', 'samples', 'version']

# for each row in header
for record in header.records:
  	# print that record
  	print(record)
		# print the name for each record
    print(record.key)
    # print everything under the key
    print(dict(record.items()))

```



A good example of using header:

```python
''' CNVkit '''

meta = collections.defaultdict(list)
for hr in vcf_reader.header.records:
	if hr.key and hr.key not in ('ALT', 'FILTER', 'FORMAT', 'INFO', 'contig'):
		meta[hr.key].append(dict(hr.items()))


if "GATKCommandLine" in meta:
  # GATK 3.0(?) and earlier
  for tag in meta["GATKCommandLine"]:
    if tag.get("ID") == "MuTect":  # any others OK?
      options = dict(kv.split("=", 1) for kv in (tag["CommandLineOptions"].strip('"').split())if '=' in kv)
      
      sample_id = options.get('tumor_sample_name')
      normal_id = options['normal_sample_name']
      yield sample_id, normal_id
```





*******

#### VariantRecord

```python
import pysam

vcf = pysam.VariantFile('MySample.vcf', 'r')
# Each item in VariantFile is a VariantRecord object
for record in vcf:
  	# attributes:
    # ['__class__', '__delattr__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__setstate__', '__sizeof__', '__str__', '__subclasshook__', 'alleles', 'alts', 'chrom', 'contig', 'copy', 'filter', 'format', 'header', 'id', 'info', 'pos', 'qual', 'ref', 'rid', 'rlen', 'samples', 'start', 'stop', 'translate']
    
    # get the whole row
    print(record)
    
    # get positions (That's the position of reference genome !)
    print(record.pos, record.start, record.stop)
    
    # get info for one sample
    SampleInfo = record.samples['MySampleID']

    # for each sample, get their info
    print(dict(SampleInfo))
    # output looks like this:
		#{'GT': (0, 1), 'AD': (71, 54), 'AF': (0.4309999942779541,), 'DP': 125, 'F1R2': (0, 0), 'F2R1': (71, 54), 'MBQ': (30, 35), 'MFRL': (227, 227), 'MMQ': (60,), 'MPOS': (16,), 'ORIGINAL_CONTIG_MISMATCH': 0, 'SA_MAP_AF': (0.4339999854564667, 0.3840000033378601, 0.4320000112056732), 'SA_POST_PROB': (0.02500000037252903, 0.02800000086426735, 0.9470000267028809)}
    
    # get Alleic depth for each called mutation
		print(SampleID['AD'])

vcf.close()
```





