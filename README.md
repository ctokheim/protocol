# Protocol

## Installation

You can pip install the package:

```bash
$ pip install git+git://github.com/ctokheim/protocol.git
$ standard_plots --help
```

The required python packages are shown in the [requirements_plotting.txt](https://github.com/ctokheim/protocol/blob/master/requirements_plotting.txt) file.

## Data files

### method results

All results of a method should be stored in the same directory (specified by `--input-dir` option). Files should be named by cancer type or PANCAN for pan-cancer results. For example, pancancer should be PANCAN.txt and LUAD.txt for lung adenocarcinoma. The files are assumed to be tab-separated.

### gene lists

There are two gene lists: Cancer gene census and cancer genome landscapes. 

```bash
$ wget https://raw.githubusercontent.com/ctokheim/protocol/master/data/Census_allSat%20Jan%20%207%2018-57-49%202017.tsv
$ wget https://raw.githubusercontent.com/ctokheim/protocol/master/data/cancer_genome_landscapes.txt
```

## Configuration file

Often methods will have different column names for the p-value/q-values. By default these will assume the p-value column name is "pvalue", and the q-value column names is "qvalue". This can be changed through a YAML configuration file. We have already provided you with a template configuration file (config.yaml) in the repository. Although, this file can be tweaked and/or further extended for additional methods. The configuration file has the following format:

```yaml
METHODNAME:
    qvalue:
        - qvalCol1
        - qvalCol2
    pvalue:
        - pvalCol1
        - pvalCol2
    threshold:
        score: 0.1
        top: low
```

METHODNAME is the name of the method, which should match the method name provided through the command line argument. In this example, the method reports two p-values (column names: "pvalCol1" and "pvalCol2") and two q-values (column names: "qvalCol1" and "qvalCol2"). In the common case that a method reports one p-value/q-value column, then only one bullet point for each attribute would be needed. The q-value threshold in this case was set at 0.1, and top significant genes are those below the threshold ("top" set to "low").

## Platform

The protocol package should work on UNIX platforms (mac os x, linux). It likely will work on windows, but has not been tested (installation instructions may differ).

## Support

Please contact the package developer Collin Tokheim (ctokheim at jhu dot edu) with suggestions, questions or bug reports.
