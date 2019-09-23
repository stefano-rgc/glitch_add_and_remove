from inspect import getfullargspec

def get_kwargs(func):
    """Return a dict of (reversed) keyword arguments from a function."""
    spec = getfullargspec(func)
    keys = reversed(spec.args)
    values = reversed(spec.defaults)
    if len(spec.defaults) > 0:
        return { k:v for k,v in zip(keys,values) }
    else:
        return {}