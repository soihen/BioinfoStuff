# Pysam — VCF

What is VCF format?:   [VCF format details](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40)



pysam is a Python module allowing easy manipulation of SAM/BAM and VCF/BCF files

Manipulation of VCF files is mainly through class `pysam.VariantFile`



The `VariantFile` class contain two sub-classes:  **VariantHeader** && **VariantRecord**

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



### VariantHeader

```
Help on VariantHeader object:

class VariantHeader(builtins.object)
 |  VariantHeader()
 |  header information for a :class:`VariantFile` object
 |  
 |  Methods defined here:
 |  
 |  __bool__(self, /)
 |      self != 0
 |  
 |  __init__(self, /, *args, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  __reduce__ = __reduce_cython__(...)
 |      VariantHeader.__reduce_cython__(self)
 |  
 |  __setstate__ = __setstate_cython__(...)
 |      VariantHeader.__setstate_cython__(self, __pyx_state)
 |  
 |  __str__(self, /)
 |      Return str(self).
 |  
 |  add_line(...)
 |      VariantHeader.add_line(self, line)
 |      Add a metadata line to this header
 |  
 |  add_meta(...)
 |      VariantHeader.add_meta(self, key, value=None, items=None)
 |      Add metadata to this header
 |  
 |  add_record(...)
 |      VariantHeader.add_record(self, VariantHeaderRecord record)
 |      Add an existing :class:`VariantHeaderRecord` to this header
 |  
 |  add_sample(...)
 |      VariantHeader.add_sample(self, name)
 |      Add a new sample to this header
 |  
 |  copy(...)
 |      VariantHeader.copy(self)
 |  
 |  merge(...)
 |      VariantHeader.merge(self, VariantHeader header)
 |  
 |  new_record(...)
 |      VariantHeader.new_record(self, contig=None, start=0, stop=0, alleles=None,     id=None, qual=None, filter=None, info=None, samples=None, **kwargs)
 |      Create a new empty VariantRecord.
 |      
 |              Arguments are currently experimental.  Use with caution and expect
 |              changes in upcoming releases.
 |  
 |  ----------------------------------------------------------------------
 |  Static methods defined here:
 |  
 |  __new__(*args, **kwargs) from builtins.type
 |      Create and return a new object.  See help(type) for accurate signature.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  alts
 |      alt metadata (:class:`dict` ID->record).
 |      
 |      The data returned just a snapshot of alt records, is created
 |      every time the property is requested, and modifications will
 |      not be reflected in the header metadata and vice versa.
 |      
 |      i.e. it is just a dict that reflects the state of alt records
 |      at the time it is created.
 |  
 |  contigs
 |      contig information (:class:`VariantHeaderContigs`)
 |  
 |  filters
 |      filter metadata (:class:`VariantHeaderMetadata`)
 |  
 |  formats
 |      format metadata (:class:`VariantHeaderMetadata`)
 |  
 |  info
 |      info metadata (:class:`VariantHeaderMetadata`)
 |  
 |  records
 |      header records (:class:`VariantHeaderRecords`)
 |  
 |  samples
 |      samples (:class:`VariantHeaderSamples`)
 |  
 |  version
 |      VCF version
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __pyx_vtable__ = <capsule object NULL>
```



### VariantRecord

```
Help on VariantRecord object:

class VariantRecord(builtins.object)
 |  VariantRecord(*args, **kwargs)
 |  Variant record
 |  
 |  Methods defined here:
 |  
 |  __eq__(self, value, /)
 |      Return self==value.
 |  
 |  __ge__(self, value, /)
 |      Return self>=value.
 |  
 |  __gt__(self, value, /)
 |      Return self>value.
 |  
 |  __init__(self, /, *args, **kwargs)
 |      Initialize self.  See help(type(self)) for accurate signature.
 |  
 |  __le__(self, value, /)
 |      Return self<=value.
 |  
 |  __lt__(self, value, /)
 |      Return self<value.
 |  
 |  __ne__(self, value, /)
 |      Return self!=value.
 |  
 |  __reduce__ = __reduce_cython__(...)
 |      VariantRecord.__reduce_cython__(self)
 |  
 |  __setstate__ = __setstate_cython__(...)
 |      VariantRecord.__setstate_cython__(self, __pyx_state)
 |  
 |  __str__(self, /)
 |      Return str(self).
 |  
 |  copy(...)
 |      VariantRecord.copy(self)
 |      return a copy of this VariantRecord object
 |  
 |  translate(...)
 |      VariantRecord.translate(self, VariantHeader dst_header)
 |  
 |  ----------------------------------------------------------------------
 |  Static methods defined here:
 |  
 |  __new__(*args, **kwargs) from builtins.type
 |      Create and return a new object.  See help(type) for accurate signature.
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  

 |  alleles
 |      tuple of reference allele followed by alt alleles
 |  
 |  alts
 |      tuple of alt alleles
 |  
 |  chrom
 |      chromosome/contig name
 |  
 |  contig
 |      chromosome/contig name
 |  
 |  filter
 |      filter information (see :class:`VariantRecordFilter`)
 |  
 |  format
 |      sample format metadata (see :class:`VariantRecordFormat`)
 |  
 |  header
 |  
 |  id
 |      record identifier or None if not available
 |  
 |  info
 |      info data (see :class:`VariantRecordInfo`)
 |  
 |  pos
 |      record start position on chrom/contig (1-based inclusive)
 |  
 |  qual
 |      phred scaled quality score or None if not available
 |  
 |  ref
 |      reference allele
 |  
 |  rid
 |      internal reference id number
 |  
 |  rlen
 |      record length on chrom/contig (aka rec.stop - rec.start)
 |  
 |  samples
 |      sample data (see :class:`VariantRecordSamples`)
 |  
 |  start
 |      record start position on chrom/contig (0-based inclusive)
 |  
 |  stop
 |      record stop position on chrom/contig (0-based exclusive)
 |  
 |  ----------------------------------------------------------------------
 |  Data and other attributes defined here:
 |  
 |  __hash__ = None
```







