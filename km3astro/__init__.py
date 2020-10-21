try:
    version = get_distribution(__name__).version
except DistributionNotFound:
    version = "unknown version"
