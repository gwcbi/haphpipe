# haphpipe

HAplotype and PHylodynamics pipeline for viral assembly, population genetics, and phylodynamics.

```bash
conda create -n haphpipe python=2.7

source activate haphpipe
conda install biopython bowtie2 flash freebayes mummer picard samtools trimmomatic trinity
conda install gatk
# # # 
gatk-register GenomeAnalysisTK-3.7.tar.bz2

pip install git+git://github.com/gwcbi/haphpipe.git
```


```bash
source activate haphpipe
hp_assemble -h
```

