# chipquery
Search for public ChIP-seq data

## General info

1. All the data (bed-files) should be downloaded in advance and stored in one directory (see "Data downloading" section).
2. While pre-processing duplicates in bed-files are removed, DB file is generated (see "Pre-processing" section).
3. Various types of search may after that applied.

## Data downloading

Data downloading consists of the following two steps:
1. Generation of report TSV-file. This file contains meta-data about experimental results available via ENCODEproject
- This may be done via Encode site: https://www.encodeproject.org/report/?type=Experiment
- In the left panel data type is to be chosen (assay, organism, assembly etc). Check, that *"ChIP-seq"* is ticked in section **"Assay"** and *"bed narrowPeak"*, *"bigBed narrowPeak"*, *"bed broadPeak"*, *"bigBed broadPeak"* are ticked in section **"Available data"**.
- Through the button columns desired fields may be chosen (it is recommended to tick all of them).
- Pressing the button **"Download TSV"** starts downloading report file.
- Obtained file may check via standard software: Excel, Libre Calc etc.

2. Run bed-files downloader (`encodeproject_files_v2.py`)

usage:
```shell
encodeproject_files_v2.py 
[--assembly [ASSEMBLY]]
[--Assay Type [ASSAY TYPE]]
[--Assay Nickname [ASSAY NICKNAME]]
[--Species [SPECIES]]
[--file_format [FILE_FORMAT]]
[--file_format_type [FILE_FORMAT_TYPE]]
[--output_type [OUTPUT_TYPE]]
```
positional arguments:
* `report`                Report.tsv file downloaded via step 1

optional arguments:
* `--assembly`		Assembly, e.g. hg19, mm10
* `--Assay Type`          Assay type, e.g. *ChIP-seq*, *RNA-seq*, *RIP-seq*, *iCLIP*, *RRBS*, *DNA-seq*
* `--Assay Nickname`      Assay type, e.g. *ChIP-seq*, *'polyA mRNA RNA-seq'*
* `--Species` 		Species, e.g. *'Homo sapiens'*, *'Mus muscus'*
* `--file_format`         Format, e.g. *fastq*, *bam*, *bed*, *bigBed*, *bigWig*
* `--file_format_type`    File format type, e.g. *narrowPeak* or *broadPeak*
* `--output_type`         reads, alignment, *'unfiltered alignments'*, *'signal p-value'*, *'peaks'*, *'replicated peaks'*,*'fold change over control'*, *'conservative idr thresholded peaks'*

Example to process all hg19 peaks in bed format: 
```shell
encodeproject_files_v2.py report.tsv --assembly hg19 --file_format bed
```

## Pre-processing

### Convert everything in sorted bed
1. Convert bigBed to bed (this will automatically generate sorted bed file):
```shell
for filename in *.bigBed; do bigBedToBed ${filename} all_sorted_bed/${filename}.bed; done
```
2. Convert every bed.gz to sorted bed
```shell
for filename in *.bed.gz; do bedtools sort -i ${filename} > all_sorted_bed/${filename}.bed; done
```

### Removing duplicates
The following procedure calculates and compares MD5 hashes for all bed-files in defined directory. Redundant files are then removed.
```shell
python3 remove_duplicates.py -d path_to_directory_with_bed_files
```

### Generating DB file
Now we can gather all the necessary information about resulting files list into DB:

```shell
python3 create_db.py -t DIR -d DIR [-r FILE] [-o OUTPUT] 
```
optional arguments:
* `-t, --tracks`	Directory, containing bed-tracks
* `-d, --diff` 	Directory, containing differnetial tracks
* `-r, --report`	Report summary CSV file
* `-o, --output` 	Name of resulting file containing DB (without extension)

Example: 
```shell
python create_db.py -t ./all_data_in_bed/tracks/ -d ./all_data_in_bed/diff_chipseq/ -r ./report_summary.csv -o summary_DB
```
## Search

```shell
find_track_db.py [-h] -q DIR -d DIR [-i FILE] [-t Int] [-k Int] -o FILE
                        [-p] [-l Str] [-b Str] [-a Str] [-m Int] [-s] [-n]`
```
optional arguments:
*  `-q, --query`   	 Directory, containing query files
*  `-d, --DB`         	 Directory, containing DB files
*  `-i, --info`    	 Database description file
*  `-t, --threads` 	 Number of threads
*  `-k`,            	 Number of the most similar entries
*  `-o, --output`          Results filename (w/o extension)
*  `-p, --pickle_ouput`    If this flag is true - save search results in pickle-file either
*  `-l, --label`   	 Histone modification label (e.g. H3K4me1, H3K27ac, FN etc.). FN means that histone modification is to be found in file name. This filter is not used by default
*  `-b, --biosample`       A string which is to present in biosample field
*  `-a, --assembly`        Assembly type: hg19/GRCh38/ALL
*  `-m, --pubmed` 	 PubMed ID (all by default)
*  `-s, --skip`            Skip files which are not indexed in DB
*  `-n, --nodb`            Search without database (search directory vs. directory)

Example 1 (compare two lists of bed-files everyone vs. everyone withput meta-data): 
```shell
python3 find_track_db.py -q ./query1 -d ./query2 -t 8 -o search_results -k 10000 -n
```
Example 2 (Find top 10 CD14 tracks for every file in directory `./query/`): 
```shell
python3 find_track_db.py -q ./query -d ./all_data_in_bed/tracks/ -i chipseq_db.pickle -t 8 -o search_results_test2 -k 10 -s -b CD14`
```
Example 2 (Find best match within its own histone label, label name is to present within filename directly e.g. `SICER_0.01_YD_H3K4me1.bed`): 
```shell
python3 find_track_db.py -q ./query -d ./all_data_in_bed/tracks/ -i chipseq_db.pickle -t 8 -o search_results_test2 -k 1 -s -l FN
```
