# Protocol

## About

A major goal of the huge public investment in large-scale cancer sequencing has been to find the majority of driver genes.  Robust computational prediction of drivers from small somatic variants is critical to this mission, and it is essential that the best methods be identified.  While many such methods have been proposed, it has been difficult to evaluate them because there is no gold standard to use as a benchmark.  Here we developed an evaluation framework for driver gene prediction methods that does not require a gold standard.  The framework includes a large set of small somatic mutations from a wide range of cancer types and five evaluation metrics.  We propose it can be used to systematically evaluate new prediction methods and compare them to existing methods.  

## Links

* [Documentation Wiki](http://github.com/KarchinLab/protocol/wiki/Home)
* [Installation](http://github.com/KarchinLab/protocol/wiki/Installation)
* [Tutorial](http://github.com/KarchinLab/protocol/wiki/Tutorial)

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
    threshold
        score: 0.1
        top: low
```

METHODNAME is the name of the method, which should match the file naming convention described in the data format above. In this example, the method reports two p-values (column names: "pvalCol1" and "pvalCol2") and two q-values (column names: "qvalCol1" and "qvalCol2"). In the common case that a method reports one p-value/q-value column, then only one bullet point for each attribute would be needed. The q-value threshold in this case was set at 0.1, and top significant genes are those below the threshold ("top" set to "low").

## Availability

Releases can be found on github at

* [http://github.com/KarchinLab/protocol/releases](http://github.com/KarchinLab/protocol/releases)

## Platform

The protocol package should work on UNIX platforms (mac os x, linux). It likely will work on windows, but has not been tested (installation instructions may differ).

## Support

Please contact the package developer Collin Tokheim (ctokheim at jhu dot edu) with suggestions, questions or bug reports.
