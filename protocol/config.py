

def fetch_level_names(config, level='gene', exclude=False):
    """Get all method names that predict genes."""
    # get methods that will be excluded
    if exclude:
        exclude_list = config.get('exclude', [])
    else:
        exclude_list = []

    gene_methods = []
    for meth_name in config:
        # skip if excluded
        if meth_name in exclude_list: continue
        # skip if not a method
        if meth_name == 'exclude': continue

        lvl = config[meth_name].get('level', None)
        if str(lvl).lower() == level:
            gene_methods.append(meth_name)
    return gene_methods


def is_valid_config(myconfig, method_name, attribute):
    # check if the above attributes are valid
    has_config = myconfig is not None
    if not has_config:
        return False
    has_method = method_name in myconfig
    if not has_method:
        return False
    has_attribute = attribute in myconfig[method_name]
    if not has_attribute:
        return False

    # no problem, return True
    return True


def fetch_threshold(config, meth_name):
    """Figure out if a score or qvalue is being used."""
    if is_valid_config(config, meth_name, 'threshold'):
        # figure out score/qvalue threshold
        if config[meth_name]['threshold'].has_key('score'):
            thresh_col = 'score'
        elif config[meth_name]['threshold'].has_key('qvalue'):
            thresh_col = 'qvalue'
        elif config[meth_name]['threshold'].has_key('pvalue'):
            thresh_col = 'pvalue'
        elif config[meth_name]['threshold'].has_key('rank'):
            thresh_col = 'rank'
        thresh_val = float(config[meth_name]['threshold'][thresh_col])
        direction = config[meth_name]['threshold']['top']
        return thresh_col, thresh_val, direction
    else:
        logger.error('No threshold attribute specified in config for {0}'.format(meth_name))
        sys.exit(1)
