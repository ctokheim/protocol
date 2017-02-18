
def fetch_gene_level_names(config):
    """Get all method names that predict genes."""
    gene_methods = []
    for meth_name in config:
        lvl = config[meth_name].get('level', None)
        if lvl.lower() == 'gene':
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
        thresh_val = float(config[meth_name]['threshold'][thresh_col])
        direction = config[meth_name]['threshold']['top']
        return thresh_col, thresh_val, direction
    else:
        logger.error('No threshold attribute specified in config for {0}'.format(meth_name))
        sys.exit(1)
