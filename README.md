# Protocol

This code is meant to create a list of consensus driver genes in cancer based on the results of several methods.

## Installation

You can pip install the package:

```bash
$ pip install git+git://github.com/ctokheim/protocol.git
$ driver_protocol --help
```

The required python packages are shown in the [requirements_plotting.txt](https://github.com/ctokheim/protocol/blob/master/requirements_plotting.txt) file. The "driver_protocol" command represents the "protcol.py" file in the source code.

## Data files

### method results

All results should be stored in the same directory (specified by `--input-dir` option). Directories should be named for the method, and the name should match that found in the configuration file (see Configuration file section). Results for each method should be named by cancer type or PANCAN for pan-cancer results. For example, pancancer should be PANCAN.txt and LUAD.txt for lung adenocarcinoma. The files are assumed to be tab-separated. It is **assumed** that you have a pancancer result (PANCAN.txt), but you may or may not include some number of cancer type specific results.

### gene lists

There are four gene lists: Cancer gene census, cancer genome landscapes, Kandoth et al Pancan12 smgs, and the Tamborero et al Pancan12 high confidence drivers (HCD). 

```bash
$ wget https://raw.githubusercontent.com/ctokheim/protocol/master/data/Census_allSat%20Jan%20%207%2018-57-49%202017.tsv
$ wget https://raw.githubusercontent.com/ctokheim/protocol/master/data/cancer_genome_landscapes.txt
$ wget https://raw.githubusercontent.com/ctokheim/protocol/master/data/hcd_pancan12.txt
$ wget https://raw.githubusercontent.com/ctokheim/protocol/master/data/kandoth_pancan12_smgs.txt
```

## Configuration file

Often methods will have different column names for the gene names, p-values, and q-values. The config.yaml file can be tweaked and/or further extended for additional methods. The configuration file has the following format:

```yaml
METHODNAME:
    level: gene
    threshold:
        qvalue: 0.1
        top: low
```

METHODNAME is the name of the method, which should match the method name provided through the command line argument. You can change the column name containing the gene name from "gene" to the column name of your method's result file. The q-value threshold in this case was set at 0.1, and top significant genes are those below the threshold ("top" set to "low"). Note, if your threshold is based on a p-value or score then use "pvalue" or "score", respectively, instead of "qvalue".

## Platform

The protocol package should work on UNIX platforms (mac os x, linux). It likely will work on windows, but has not been tested (installation instructions may differ).

## Support

Please contact the package developer Collin Tokheim (ctokheim at jhu dot edu) with suggestions, questions or bug reports.
